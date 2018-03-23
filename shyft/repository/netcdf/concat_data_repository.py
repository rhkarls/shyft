from os import path
import numpy as np
from netCDF4 import Dataset
from shyft import api
from shyft import shyftdata_dir
from .. import interfaces
from .time_conversion import convert_netcdf_time
from .utils  import calc_RH, _limit_1D, _numpy_to_geo_ts_vec
import warnings
from shyft.repository.interfaces import ForecastSelectionCriteria

UTC = api.Calendar()


class ConcatDataRepositoryError(Exception):
    pass


class ConcatDataRepository(interfaces.GeoTsRepository):
    """
    Repository for geo located forecasts in concatenated netCDF4 files.
    NetCDF4 dataset assumptions:
        * Dimensions:
            * time (unlimited)
            * lead_time
            * ensemble_member (optional)
            * grid_point
        * Variables:
            * time:(float) array with periodic forecast creation timestamps in seconds since (1970.01.01 00:00, UTC)
            * lead_time:(float) array with hours since creation time
            * x: (float) array of latitudes of dims (grid_point)
            * y: (float) array of longitudes of dims (grid_point)
            * z: (float) array of altitudes [m] of dims (grid_point)
            * forecast_is_complete: flag array of dims (time)
            * crs: has attribute proj4, a string describing the coordinate system
        * Optional variables:
            * dew_point_temperature_2m: [K],
            * surface_air_pressure: [Pa],
            * relative_humidity_2m: [1],
            * air_temperature_2m: [K],
            * precipitation_amount: [kg/m^2],
            * precipitation_amount_acc: [kg/m^2],
            * x_wind_10m: [m/s],
            * y_wind_10m: [m/s],
            * windspeed_10m: [m/s],
            * integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time: [W s/m^2]}
            * All optional variables are (float) array with dims (time, lead_time, [ensemble_member], grid_point)
    """

    _G = 9.80665  # WMO-defined gravity constant to calculate the height in metres from geopotential

    def __init__(self, epsg, filename, nb_pads=0, nb_lead_intervals_to_drop=0, nb_lead_intervals=None, fc_periodicity=1,
                 ensemble_member=0, padding=5000., use_filled_values=False):
        """
        Construct the netCDF4 dataset reader for concatenated gridded forecasts and initialize data retrieval.

        Parameters
        ----------
        epsg: string
            Unique coordinate system id for result coordinates. Currently "32632" and "32633" are supported.
        filename: string
            Path to netcdf file containing concatenated forecasts
        nb_pads: int, optional
            Extend forecast horizon with nb_pads lead times by replicating last values
        nb_lead_intervals_to_drop: int, optional
            Index of first lead time to read
        nb_lead_intervals: int, optional
            Number of lead time intervals to read
        fc_periodicity: int, optional
            Periodicity of forecast to read (i.e. if fc_periodicity=2 only reads values from every second forecast)
        ensemble_member: int, optional
            Ensemble member returned by get_timeseries, get_forecast and get_forecast_collection if dataset is of
            ensemble type (has dimension 'ensemble_member'). Must be non-negative integer less than dimension size of
            'ensemble_member*
        padding: float, optional
            Longidutinal and latitudinal padding in meters, added east, west, north and south
        use_filled_values: bool
            Use forecast with filled valued if True
        """

        filename = path.expandvars(filename)
        if not path.isabs(filename):
            # Relative paths will be prepended the data_dir
            filename = path.join(shyftdata_dir, filename)
        if not path.isfile(filename):
            raise ConcatDataRepositoryError("No such file '{}'".format(filename))
        self._filename = filename
        self.nb_pads = nb_pads
        self.nb_lead_intervals_to_drop = nb_lead_intervals_to_drop  # index of first lead time: starts from 0
        self.nb_lead_intervals = nb_lead_intervals
        self.fc_periodicity = fc_periodicity
        if ensemble_member is None:
            self.ensemble_member = 0
        else:
            self.ensemble_member = ensemble_member
        self.padding = padding
        self.shyft_cs = "+init=EPSG:{}".format(epsg)
        self.use_filled_values = use_filled_values

        with Dataset(self._filename) as dataset:
            self._get_time_structure_from_dataset(dataset)

        # Mapping netcdf_name: shyft_name. See also _transform_raw which contains associated transformations methods
        self._shyft_map = {"dew_point_temperature_2m": "dew_point_temperature_2m",
                           "surface_air_pressure": "surface_air_pressure",
                           "relative_humidity_2m": "relative_humidity",
                           "air_temperature_2m": "temperature",
                           "precipitation_amount": "precipitation",
                           "precipitation_amount_acc": "precipitation",
                           "x_wind_10m": "x_wind",
                           "y_wind_10m": "y_wind",
                           "windspeed_10m": "wind_speed",
                           "integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time": "radiation"}

        # Assumed unit for fields. These are checked against netcdf content in _validate_input
        self.var_units = {"dew_point_temperature_2m": ['K'],
                          "surface_air_pressure": ['Pa'],
                          "relative_humidity_2m": ['1'],
                          "air_temperature_2m": ['K'],
                          "precipitation_amount": ['kg/m^2'],
                          "precipitation_amount_acc": ['kg/m^2'],
                          "x_wind_10m": ['m/s'],
                          "y_wind_10m": ['m/s'],
                          "windspeed_10m": ['m/s'],
                          "integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time": ['W s/m^2']}

        # Fields that need an additional lead_time to de-accumulate
        self._shift_fields = ("precipitation_amount_acc",
                              "integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time")

        # geo-ts creation method
        self.create_geo_ts_type_map = {"relative_humidity": api.create_rel_hum_source_vector_from_np_array,
                                      "temperature": api.create_temperature_source_vector_from_np_array,
                                      "precipitation": api.create_precipitation_source_vector_from_np_array,
                                      "radiation": api.create_radiation_source_vector_from_np_array,
                                      "wind_speed": api.create_wind_speed_source_vector_from_np_array}

        # ts point interpretation policy
        self.series_type = {"relative_humidity": api.POINT_INSTANT_VALUE,
                            "temperature": api.POINT_INSTANT_VALUE,
                            "precipitation": api.POINT_AVERAGE_VALUE,
                            "radiation": api.POINT_AVERAGE_VALUE,
                            "wind_speed": api.POINT_INSTANT_VALUE}

    def _get_time_structure_from_dataset(self, dataset):
        """
        This function performs the following:
            * Get time (time stamps) and lead_time (hours) from NetCDF4 dataset.
            * Calculate lead times in seconds(lead_times_in_sec).
            * Determine number of lead times to use (fc_len_to_concat) when
              concated result is asked for by get_timeseries. Note that
              fc_len_to_concat is a function of the time step length in
              time (recall we assume periodic time steps), nb_lead_intervals_to_drop,
              fc_periodicity and resolution of lead_time
            * Validates ensemble_member
        """
        nb_lead_intervals_to_drop = self.nb_lead_intervals_to_drop
        fc_periodicity = self.fc_periodicity
        time = dataset.variables.get("time", None)
        lead_time = dataset.variables.get("lead_time", None)
        if not all([time, lead_time]):
            raise ConcatDataRepositoryError("Something is wrong with the dataset"
                                            "time or lead_time not found")
        time = convert_netcdf_time(time.units, time)
        self.time = time
        self.forecast_is_complete = dataset.variables.get("forecast_is_complete", None)
        self.lead_time = lead_time[:]
        self.lead_times_in_sec = lead_time[:] * 3600.
        if self.nb_lead_intervals is None:
            self.nb_lead_intervals = len(self.lead_times_in_sec) - nb_lead_intervals_to_drop
        nb_lead_intervals = self.nb_lead_intervals
        if nb_lead_intervals_to_drop + nb_lead_intervals > len(self.lead_times_in_sec):
            raise  ConcatDataRepositoryError("'nb_lead_intervals_to_drop' + 'nb_lead_intervals' is too large")
        time_shift_with_drop = time + self.lead_times_in_sec[nb_lead_intervals_to_drop]
        # TODO: Errorhandling for idx_max required?
        idx_max = np.argmax(time[0] + self.lead_times_in_sec >= time_shift_with_drop[fc_periodicity])
        self.fc_len_to_concat = idx_max - nb_lead_intervals_to_drop

        # Validate ensemble_member
        if 'ensemble_member' in dataset.dimensions:
            nb_ensemble_member = dataset.dimensions['ensemble_member'].size
            if self.ensemble_member >= nb_ensemble_member:
                raise ConcatDataRepositoryError(
                    "ensemble_member must be non-negative integer between 0 and {}".format(nb_ensemble_member))

    def get_timeseries_ensemble(self, input_source_types, utc_period, geo_location_criteria=None):
        """
        Get ensemble of shyft source vectors of time series covering utc_period
        for input_source_types.

        Time series are constructed by concatenating values from forecasts
        according to fc_periodicity whose lead period
        (nb_lead_intervals_to_drop, nb_lead_intervals_to_drop + fc_len_to_concat)
        intersect the utc_period. See _get_time_structure_from_dataset for details
        on fc_len_to_concat.


        Parameters
        ----------
        see interfaces.GeoTsRepository

        Returns
        -------
        see interfaces.GeoTsRepository
        """
        no_shift_fields = set([self._shyft_map[k] for k in self._shift_fields]).isdisjoint(input_source_types)
        if self.fc_len_to_concat < 0 \
            or no_shift_fields and self.nb_lead_intervals_to_drop + self.fc_len_to_concat > len(self.lead_time) \
            or not no_shift_fields and self.nb_lead_intervals_to_drop + self.fc_len_to_concat + 1 > len(self.lead_time):
                err_msg = "'nb_lead_intervals_to_drop={}' is too large for concatenation"\
                          .format(self.nb_lead_intervals_to_drop)
                raise ConcatDataRepositoryError(err_msg)
        with Dataset(self._filename) as dataset:
            # fc_selection_criteria ={'forecasts_with_start_within_period': utc_period}
            fc_selection_criteria = ForecastSelectionCriteria(forecasts_with_start_within_period=utc_period)
            extracted_data, x, y, z = self._get_data_from_dataset(dataset, input_source_types, fc_selection_criteria,
                                                                  geo_location_criteria,
                                                                  nb_lead_intervals=self.fc_len_to_concat, concat=True)
            # check if extra_intervals are required
            ta = list(extracted_data.values())[0][1] # time axis of first item
            ta_end = ta.total_period().end
            if ta_end < utc_period.end: # try to extend extracted data with remainder of last forecast
                sec_to_extend = utc_period.end - ta_end
                drop = self.nb_lead_intervals_to_drop + self.fc_len_to_concat
                idx = np.argmax(self.lead_times_in_sec[drop:] - self.lead_times_in_sec[drop] >= sec_to_extend)
                if idx == 0:
                    raise ConcatDataRepositoryError("The latest time in repository is earlier than the end of the "
                                                    "period for which data is requested")
                fc_selection_criteria = ForecastSelectionCriteria(latest_available_forecasts=
                                                             {'number_of_forecasts': 1, 'forecasts_older_than': ta_end})
                extra_data, _,_,_ = self._get_data_from_dataset(dataset, input_source_types, fc_selection_criteria,
                                    geo_location_criteria, nb_lead_intervals_to_drop=drop, nb_lead_intervals=idx,
                                    concat=False) # note: no concat here
                ta_extra = list(extra_data.values())[0][1][0]
                ta = api.TimeAxis(api.UtcTimeVector.from_numpy(np.append(ta.time_points, ta_extra.time_points[1:])))
                extracted_data = {k: (np.concatenate((extracted_data[k][0], np.squeeze(extra_data[k][0], axis=0))), ta)
                                                                                for k in list(extracted_data.keys())}
            # return self._convert_to_geo_timeseries(extracted_data, geo_pts, concat=True)
            return _numpy_to_geo_ts_vec(extracted_data, x, y, z, ConcatDataRepositoryError)

    def get_timeseries(self, input_source_types, utc_period, geo_location_criteria=None):
        """
        Same as get_timeseries_ensemble, but only returns one ensemble member
        as specified by ensemble_member.

        Parameters
        ----------
        see interfaces.GeoTsRepository

        Returns
        -------
        see interfaces.GeoTsRepository
        """
        # return self.get_timeseries_ensemble(input_source_types, utc_period, geo_location_criteria)[self.ensemble_member]
        # Same as get_timeseries_collection, but for only one member
        no_shift_fields = set([self._shyft_map[k] for k in self._shift_fields]).isdisjoint(input_source_types)
        if self.fc_len_to_concat < 0 \
            or no_shift_fields and self.nb_lead_intervals_to_drop + self.fc_len_to_concat > len(self.lead_time) \
            or not no_shift_fields and self.nb_lead_intervals_to_drop + self.fc_len_to_concat + 1 > len(self.lead_time):
                err_msg = "'nb_lead_intervals_to_drop={}' is too large for concatenation"\
                          .format(self.nb_lead_intervals_to_drop)
                raise ConcatDataRepositoryError(err_msg)
        with Dataset(self._filename) as dataset:
            # fc_selection_criteria ={'forecasts_with_start_within_period': utc_period}
            fc_selection_criteria = ForecastSelectionCriteria(forecasts_with_start_within_period=utc_period)
            extracted_data, x, y, z = self._get_data_from_dataset(dataset, input_source_types, fc_selection_criteria,
                                                                  geo_location_criteria,
                                                                  nb_lead_intervals=self.fc_len_to_concat, concat=True,
                                                                  ensemble_member=self.ensemble_member)
            # check if extra_intervals are required
            ta = list(extracted_data.values())[0][1] # time axis of first item
            ta_end = ta.total_period().end
            if ta_end < utc_period.end: # try to extend extracted data with remainder of last forecast
                sec_to_extend = utc_period.end - ta_end
                drop = self.nb_lead_intervals_to_drop + self.fc_len_to_concat
                idx = np.argmax(self.lead_times_in_sec[drop:] - self.lead_times_in_sec[drop] >= sec_to_extend)
                if idx == 0:
                    raise ConcatDataRepositoryError("The latest time in repository is earlier than the end of the "
                                                    "period for which data is requested")
                fc_selection_criteria = ForecastSelectionCriteria(latest_available_forecasts=
                                                             {'number_of_forecasts': 1, 'forecasts_older_than': ta_end})
                extra_data, _, _,_ = self._get_data_from_dataset(dataset, input_source_types, fc_selection_criteria,
                                    geo_location_criteria, nb_lead_intervals_to_drop=drop, nb_lead_intervals=idx,
                                    concat=False, ensemble_member=self.ensemble_member) # note: no concat here
                ta_extra = list(extra_data.values())[0][1][0]
                ta = api.TimeAxis(api.UtcTimeVector.from_numpy(np.append(ta.time_points, ta_extra.time_points[1:])))
                extracted_data = {k: (np.concatenate((extracted_data[k][0], np.squeeze(extra_data[k][0], axis=0))), ta)
                                                                                for k in list(extracted_data.keys())}
            #return self._convert_to_geo_timeseries(extracted_data, geo_pts, concat=True)[0]
            return _numpy_to_geo_ts_vec(extracted_data, x, y, z, ConcatDataRepositoryError)

    def get_forecast_ensemble_collection(self, input_source_types, fc_selection_criteria, geo_location_criteria=None):
        """
        Get all forecast (as shyft source vectors of time series) according
        to fc_periodicity that meet fc_selection_criteria for input_source_types.
        Only values for lead period (nb_lead_intervals_to_drop, nb_lead_intervals)
        are used.

        Parameters
        ----------
        see interfaces.GeoTsRepository

        Returns
        -------
        see interfaces.GeoTsRepository
        """

        k, v = fc_selection_criteria.criterion
        if k == 'forecasts_at_reference_times':
            fsc = lambda x: ForecastSelectionCriteria(latest_available_forecasts=
                                              {'number_of_forecasts': 1, 'forecasts_older_than':x})
            return [self.get_forecast_ensemble_collection(input_source_types, fsc(t_c), geo_location_criteria)[0]
                    for t_c in v]
        else:
            with Dataset(self._filename) as dataset:
                data, x, y, z = self._get_data_from_dataset(dataset, input_source_types, fc_selection_criteria,
                                                                      geo_location_criteria, concat=False)
                #return self._convert_to_geo_timeseries(data, geo_pts, concat=False)
                return _numpy_to_geo_ts_vec(data, x, y, z, ConcatDataRepositoryError)

    def get_forecast_collection(self, input_source_types, fc_selection_criteria, geo_location_criteria=None):
        """
        Same as get_forecast_ensemble_collection, but only returns one
        ensemble member as specified by ensemble_member.

        Parameters
        ----------
        see interfaces.GeoTsRepository

        Returns
        -------
        see interfaces.GeoTsRepository.
        """
        # return [fcst[self.ensemble_member] for fcst in
        #      self.get_forecast_ensemble_collection(input_source_types, fc_selection_criteria, geo_location_criteria)]
        k, v = fc_selection_criteria.criterion
        if k == 'forecasts_at_reference_times':
            fsc = lambda x: ForecastSelectionCriteria(latest_available_forecasts=
                                              {'number_of_forecasts': 1, 'forecasts_older_than':x})
            return [self.get_forecast_collection(input_source_types, fsc(t_c), geo_location_criteria)[0] for t_c in v]
        else:
            with Dataset(self._filename) as dataset:
                data, x, y, z = self._get_data_from_dataset(dataset, input_source_types, fc_selection_criteria,
                                                            geo_location_criteria, concat=False,
                                                            ensemble_member=self.ensemble_member)
                #return [fcst[0] for fcst in self._convert_to_geo_timeseries(data, geo_pts, concat=False)]
                return [fcst[0] for fcst in _numpy_to_geo_ts_vec(data, x, y, z, ConcatDataRepositoryError)]

    def get_forecast_ensemble(self, input_source_types, utc_period,
                              t_c, geo_location_criteria=None):
        """
        Get latest forecast (as shyft source vectors of time series) before
        t_c where the lead period (nb_lead_intervals_to_drop, nb_lead_intervals)
        covers utc_period. Only forecasts according to fc_periodicity are
        considered.

        Parameters
        ----------
        see interfaces.GeoTsRepository

        Returns
        -------
        see interfaces.GeoTsRepository
        """
        with Dataset(self._filename) as dataset:
            fsc = ForecastSelectionCriteria(forecasts_that_cover_period=utc_period)
            time_slice, lead_time_slice, m_t = self._make_time_slice(self.nb_lead_intervals_to_drop,
                                                                     self.nb_lead_intervals, fsc)
            time = self.time[time_slice][m_t[time_slice]]
            ref_time = time[np.argmin(time <= t_c) - 1]
            if ref_time.size == 0:
                raise ConcatDataRepositoryError(
                    "Not able to find forecast that cover the requested period with the provided restrictions "
                    "'nb_lead_intervals_to_drop'={}, 'nb_lead_intervals'={}, 'fc_periodicity'={} and 't_c'={}". \
                        format(self.nb_lead_intervals_to_drop, self.nb_lead_intervals, self.fc_periodicity, t_c))
            fsc = ForecastSelectionCriteria(forecasts_at_reference_times=[int(ref_time)])
            return self.get_forecast_ensemble_collection(input_source_types, fsc, geo_location_criteria)[0]

    def get_forecast(self, input_source_types, utc_period, t_c, geo_location_criteria=None):
        """
        Same as get_forecast_ensemble, but only returns one
        ensemble member as specified by ensemble_member.

        Parameters
        ----------
        see interfaces.GeoTsRepository

        Returns
        -------
        see interfaces.GeoTsRepository
        """
        with Dataset(self._filename) as dataset:
            fsc = ForecastSelectionCriteria(forecasts_that_cover_period=utc_period)
            time_slice, lead_time_slice, m_t = self._make_time_slice(self.nb_lead_intervals_to_drop,
                                                                     self.nb_lead_intervals, fsc)
            time = self.time[time_slice][m_t[time_slice]]
            ref_time = time[np.argmin(time <= t_c) - 1]
            if ref_time.size == 0:
                raise ConcatDataRepositoryError(
                    "Not able to find forecast that cover the requested period with the provided restrictions "
                    "'nb_lead_intervals_to_drop'={}, 'nb_lead_intervals'={}, 'fc_periodicity'={} and 't_c'={}". \
                        format(self.nb_lead_intervals_to_drop, self.nb_lead_intervals, self.fc_periodicity, t_c))
            fsc = ForecastSelectionCriteria(forecasts_at_reference_times=[int(ref_time)])
            return self.get_forecast_collection(input_source_types, fsc, geo_location_criteria)[0]

    def _get_data_from_dataset(self, dataset, input_source_types, fc_selection_criteria, geo_location_criteria,
                               nb_lead_intervals_to_drop=None, nb_lead_intervals=None, concat=False, ensemble_member=None):
        """
        Get data from NetCDF4 dataset.

        Parameters
        ----------
        See _init and get_forecast_ensemble_collection for description of dataset, input_source_types,
        fc_selection_criteria, geo_location_criteria, nb_lead_intervals_to_drop, nb_lead_intervals
        concat: bool, optional
            Determines transformation method (see _transform_raw)
        ensemble_member: int, optional
            Determines ensemble_member to get. If None get all members

        Returns
        -------
        tuple (extracted_data, geo_pts) where:
            extracted_data: ts-type keyed dict of tuples (d, t) where d is numpy array of shape (nb_forecasts,
                            nb_lead_intervals, nb_ensemble_members, nb_geo_pts) containing extracted values and t is
                            list a of nb_forecasts time_axis. Here:
                                nb_forecast = number of forecasts extracted
                                nb_lead_intervals = number of lead intervals for which data has been extracted
                                nb_ensemble_members = number of extracted ensemble members
                                nb_geo_pts = number of geo-points for which data has been extracted
            geo_pts: GeoPointVector
        """

        # validate input and adjust input_source_types
        input_source_types, no_temp, rh_not_ok, no_x_wind, no_y_wind, wind_not_ok = \
            self._validate_input(dataset, input_source_types, geo_location_criteria)

        # find geo_slice for slicing dataset
        x, y, z, m_xy, xy_slice, dim_grid = self._get_geo_slice(dataset, geo_location_criteria)

        # find time_slice and lead_time_slice for slicing dataset
        lead_times_in_sec = self.lead_times_in_sec
        if nb_lead_intervals_to_drop is None:
            nb_lead_intervals_to_drop = self.nb_lead_intervals_to_drop
        if concat:
            nb_lead_intervals = self.fc_len_to_concat
        elif nb_lead_intervals is None:
            if self.nb_lead_intervals is None:
                nb_lead_intervals = len(lead_times_in_sec) - nb_lead_intervals_to_drop
            else:
                nb_lead_intervals = self.nb_lead_intervals
        if nb_lead_intervals_to_drop + nb_lead_intervals > len(lead_times_in_sec):
            raise  ConcatDataRepositoryError("'nb_lead_intervals_to_drop' + 'nb_lead_intervals' i too large")
        issubset = True if len(lead_times_in_sec) > nb_lead_intervals_to_drop + nb_lead_intervals + 1 else False
        time_slice, lead_time_slice, m_t = \
            self._make_time_slice(nb_lead_intervals_to_drop, nb_lead_intervals, fc_selection_criteria)
        time = self.time[time_slice][m_t[time_slice]]
        forecast_is_complete = self.forecast_is_complete[time_slice][m_t[time_slice]]
        if np.size(time)==0:
            raise ConcatDataRepositoryError('No forecast found that meet conditions')

        # Get data by slicing into dataset
        raw_data = {}
        for k in dataset.variables.keys():
            if self._shyft_map.get(k, None) in input_source_types:
                if k in self._shift_fields and issubset:  # Add one to lead_time slice
                    data_lead_time_slice = slice(lead_time_slice.start, lead_time_slice.stop + 1)
                else:
                    data_lead_time_slice = lead_time_slice

                data = dataset.variables[k]
                dims = data.dimensions
                data_slice = len(data.dimensions) * [slice(None)]
                if 'ensemble_member' in dims and ensemble_member is not None:
                    data_slice[dims.index("ensemble_member")] = slice(ensemble_member, ensemble_member + 1)

                data_slice[dims.index(dim_grid)] = xy_slice
                data_slice[dims.index("lead_time")] = data_lead_time_slice
                data_slice[dims.index("time")] = time_slice  # data_time_slice
                xy_slice_mask = [m_xy[xy_slice] if dim == dim_grid else slice(None) for dim in dims]
                time_slice_mask = [m_t[time_slice] if dim == 'time' else slice(None) for dim in dims]

                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", message="invalid value encountered in greater")
                    warnings.filterwarnings("ignore", message="invalid value encountered in less_equal")
                    pure_arr = data[data_slice][xy_slice_mask][time_slice_mask]

                if 'ensemble_member' not in dims: # add axis for 'ensemble_member'
                    pure_arr = pure_arr[:,:,np.newaxis,:]

                if isinstance(pure_arr, np.ma.core.MaskedArray):
                    pure_arr = pure_arr.filled(np.nan)
                if not self.use_filled_values:
                    pure_arr[forecast_is_complete == 0,:,:,:] = np.nan
                # TODO: send out warning with list of filled values
                if np.isnan(pure_arr).any():
                    print("NaN found in pure_arr for {} see indices {}".format(k, np.unravel_index(np.argmax(np.isnan(pure_arr)), pure_arr.shape)))

                raw_data[self._shyft_map[k]] = pure_arr, k

        # Calculate wind speed if required
        if wind_not_ok:
            if set(("x_wind", "y_wind")).issubset(raw_data):
                x_wind, _ = raw_data.pop("x_wind") if no_x_wind else raw_data["x_wind"]
                y_wind, _ = raw_data.pop("y_wind") if no_y_wind else raw_data["y_wind"]
                ncf_name_wind = next((n_nm for n_nm, s_nm in self._shyft_map.items() if s_nm == "wind_speed"), None)
                raw_data["wind_speed"] = np.sqrt(np.square(x_wind) + np.square(y_wind)), ncf_name_wind
            else:
                raise ConcatDataRepositoryError("Not able to retrieve wind_speed from dataset")

        # Calculate relative humidity if required
        if rh_not_ok:
            if set(("surface_air_pressure", "dew_point_temperature_2m")).issubset(raw_data):
                sfc_p, _ = raw_data.pop("surface_air_pressure")
                dpt_t, _ = raw_data.pop("dew_point_temperature_2m")
                sfc_t, _ = raw_data.pop("temperature") if no_temp else raw_data["temperature"]
                ncf_name_rh = next((n_nm for n_nm, s_nm in self._shyft_map.items() if s_nm == "relative_humidity"), None)
                raw_data["relative_humidity"] = calc_RH(sfc_t, dpt_t, sfc_p), ncf_name_rh
            else:
                raise ConcatDataRepositoryError("Not able to retrieve relative_humidity from dataset")

        data_lead_time_slice = slice(lead_time_slice.start, lead_time_slice.stop + 1)
        extracted_data = self._transform_raw(raw_data, time, lead_times_in_sec[data_lead_time_slice], concat)
        return extracted_data, x, y, z

    def _validate_input(self, dataset, input_source_types, geo_location_criteria):
        """
        * Validate geo_location criteria.
        * Check units are ok
        * Check if ekstra variables are required to compute wind-speed
          and relative humidity and if so return flags
        """

        # Check units match
        unit_ok = {k: dataset.variables[k].units in self.var_units[k]
                   for k in dataset.variables.keys() if self._shyft_map.get(k, None) in input_source_types}
        if not all(unit_ok.values()):
            raise ConcatDataRepositoryError("The following variables have wrong unit: {}.".format(
                ', '.join([k for k, v in unit_ok.items() if not v])))

        # Need extra variables to calculate wind speed if not available in dataset
        no_x_wind = False if "x_wind" in input_source_types else True
        no_y_wind = False if "y_wind" in input_source_types else True
        ncf_nm = next((n_nm for n_nm, s_nm in self._shyft_map.items() if s_nm == "wind_speed"), None)
        wind_not_ok = "wind_speed" in input_source_types and (ncf_nm is None or ncf_nm not in dataset.variables)
        if wind_not_ok:
            if not isinstance(input_source_types, list):
                input_source_types = list(input_source_types)  # We change input list, so take a copy
            input_source_types.remove("wind_speed")
            input_source_types.extend(["x_wind", "y_wind"])

        # Need extra variables to calculate relative humidity if not available in dataset
        no_temp = False if "temperature" in input_source_types else True
        ncf_nm = next((n_nm for n_nm, s_nm in self._shyft_map.items() if s_nm == "relative_humidity"), None)
        rh_not_ok = "relative_humidity" in input_source_types and (ncf_nm is None or ncf_nm not in dataset.variables)
        if rh_not_ok:
            if not isinstance(input_source_types, list):
                input_source_types = list(input_source_types)  # We change input list, so take a copy
            input_source_types.remove("relative_humidity")
            input_source_types.extend(["surface_air_pressure", "dew_point_temperature_2m"])
            if no_temp: input_source_types.extend(["temperature"])

        return input_source_types, no_temp, rh_not_ok, no_x_wind, no_y_wind, wind_not_ok

    def _get_geo_slice(self, dataset, geo_location_criteria):
        """
        * Return name of "geo" dimension in dataset (dim_grid)
        * Return limiting slice in dim_grid that covers all points
          defined by geo_location_criteria (xy_slixe)
        * Return mask that identifies points defined by geo_location_criteria
        * Return points defined by geo_location_criteria as GeoPointVector
        """
        # Find xy slicing and z
        x = dataset.variables.get("x", None)
        y = dataset.variables.get("y", None)
        dim_grid = [dim for dim in dataset.dimensions if dim not in ['time', 'lead_time', 'ensemble_member']][0]
        if not all([x, y]):
            raise ConcatDataRepositoryError("Something is wrong with the dataset, x/y coords or time not found")
        data_cs = dataset.variables.get("crs", None)
        if data_cs is None:
            raise ConcatDataRepositoryError("No coordinate system information in dataset.")
        x, y, m_xy, xy_slice = _limit_1D(x[:], y[:], data_cs.proj4, self.shyft_cs, geo_location_criteria, self.padding, ConcatDataRepositoryError)

        # Find height
        if 'z' in dataset.variables.keys():
            data = dataset.variables['z']
            dims = data.dimensions
            data_slice = len(data.dimensions) * [slice(None)]
            data_slice[dims.index(dim_grid)] = m_xy
            z = data[data_slice]
        else:
            raise ConcatDataRepositoryError("No elevations found in dataset")

        return x, y, z, m_xy, xy_slice, dim_grid

    def _make_time_slice(self, nb_lead_intervals_to_drop, nb_lead_intervals, fc_selection_criteria):
        """
        Find limiting slices in dimensions 'time' and 'lead_time'
        for filters as defined by fc_selection_criteria,
        nb_lead_intervals_to_drop, nb_lead_intervals and nb_periodicity

        Returns
        -------
        * time_slice: limiting slice in dim 'time'
        * lead_time_slice: limiting slice in dim 'lead_time'
        * m_t: periodicity mask in dim 'time'
        """
        time = self.time
        lead_times_in_sec = self.lead_times_in_sec
        # Find periodicity mask
        fc_periodicity = self.fc_periodicity
        m_t = np.zeros(time.shape, dtype=bool)
        m_t[::-fc_periodicity] = True # newest forecast last

        # k, v = list(fc_selection_criteria.items())[0]
        k, v = fc_selection_criteria.criterion
        nb_extra_intervals = 0
        if k == 'forecasts_created_within_period':
            time_slice = ((time >= v.start) & (time <= v.end))
            if not any(time_slice):
                raise ConcatDataRepositoryError(
                    "No forecasts found with creation time within period {}.".format(v.to_string()))
        elif k == 'forecasts_with_start_within_period':
            # shift utc period with nb_fc_to drop
            v_shift = api.UtcPeriod(int(v.start - lead_times_in_sec[nb_lead_intervals_to_drop]),
                                    int(v.end - lead_times_in_sec[nb_lead_intervals_to_drop]))
            time_slice = ((time >= v_shift.start) & (time <= v_shift.end))
            if not any(time_slice):
                raise ConcatDataRepositoryError(
                    "No forecasts found that start within period for period {} and "
                    "'fc_nb_to_drop'={}.".format(v.to_string(),nb_lead_intervals_to_drop))
        elif k == 'forecasts_that_cover_period':
            v_shift_first = int(v.start - lead_times_in_sec[nb_lead_intervals_to_drop]) # shift utc period with nb_fc_to drop
            # shift utc period with nb_fc_to drop + nb_lead_intervals
            v_shift_last = int(v.end - lead_times_in_sec[nb_lead_intervals_to_drop + nb_lead_intervals - 1])
            time_slice = ((time <= v_shift_first) & (time >= v_shift_last))
            if not any(time_slice):
                raise ConcatDataRepositoryError(
                    "No forecasts found that cover period {} with restrictions 'fc_nb_to_drop'={} "
                    "and 'nb_lead_intervals'={}.".format(v.to_string(), nb_lead_intervals_to_drop, nb_lead_intervals))
        elif k == 'latest_available_forecasts':
            t = v['forecasts_older_than']
            n = v['number_of_forecasts']
            idx = np.argmin(time <= t) - 1
            if idx < 0:
                first_lead_time_of_last_fc = int(time[-1])
                if first_lead_time_of_last_fc <= t:
                    idx = len(time) - 1
                else:
                    raise ConcatDataRepositoryError(
                        "The earliest time in repository ({}) is later than or at the start of the period for which data is "
                        "requested ({})".format(UTC.to_string(int(time[0])), UTC.to_string(t)))
            if idx + 1 < n * fc_periodicity:
                raise ConcatDataRepositoryError(
                    "The number of forecasts available in repo ({}) and earlier than the parameter "
                    "'forecasts_older_than' ({}) is less than the number of forecasts requested ({}) " ""
                    "for the specified periodicity ({})".format(idx + 1, UTC.to_string(t), n, fc_periodicity))
            time_slice = slice(idx - n * fc_periodicity + 1, idx + 1)
        elif k == 'forecasts_at_reference_times':
            raise ConcatDataRepositoryError(
                "'forecasts_at_reference_times' selection criteria not supported yet.")
        lead_time_slice = slice(nb_lead_intervals_to_drop, nb_lead_intervals_to_drop + nb_lead_intervals)
        return time_slice, lead_time_slice, m_t

    def _transform_raw(self, data, time, lead_time, concat):
        # TODO: check robustness off all conversion for flexible lead_times
        """
        We need full time if deaccumulating
        """

        # time axis for contatenated output
        def concat_t(t):
            t_stretch = np.ravel(np.repeat(t, self.fc_len_to_concat).reshape(len(t), self.fc_len_to_concat)
                                 + lead_time[0:self.fc_len_to_concat])
            dt_last = lead_time[-1] - lead_time[-2]
            if np.all(lead_time[1:]-lead_time[:-1] == dt_last): # fixed_dt time axis
                return api.TimeAxis(int(t_stretch[0]), int(t_stretch[1]) - int(t_stretch[0]), len(t_stretch))
            else: # point_type time axis
                return api.TimeAxis(api.UtcTimeVector.from_numpy(t_stretch.astype(int)), int(t_stretch[-1] + dt_last))

        # time axis for forecast output
        def forecast_t(t, daccumulated_var=False):
            nb_ext_lead_times = len(lead_time) - 1 if daccumulated_var else len(lead_time)
            t_all = np.repeat(t, nb_ext_lead_times).reshape(len(t), nb_ext_lead_times) + lead_time[0:nb_ext_lead_times]
            return t_all.astype(int)

        # padding extra points to forecast
        def pad(v, t):
            if not concat: # Extend forecast by duplicating last nb_pad values
                if self.nb_pads > 0:
                    nb_pads = self.nb_pads
                    t_padded = np.zeros((t.shape[0], t.shape[1] + nb_pads), dtype=t.dtype)
                    t_padded[:, :-nb_pads] = t[:, :]
                    t_add = t[0, -1] - t[0, -nb_pads - 1]
                    # print('t_add:',t_add)
                    t_padded[:, -nb_pads:] = t[:, -nb_pads:] + t_add

                    v_padded = np.zeros((v.shape[0], t.shape[1] + nb_pads, v.shape[2], v.shape[3]), dtype=v.dtype)
                    v_padded[:, :-nb_pads, :, :] = v[:, :, :, :]
                    v_padded[:, -nb_pads:, :, :] = v[:, -nb_pads:, :, :]
                else:
                    t_padded = t
                    v_padded = v
                dt_last = t_padded[0, -1] - t_padded[0, -2]
                if np.all(t_padded[:,1:] - t_padded[:,:-1] == dt_last): # fixed_dt time axis
                    return (v_padded, [api.TimeAxis(int(tp[0]), int(dt_last), len(tp)) for tp in t_padded])
                else: # point_type time axis
                    return (v_padded,
                            [api.TimeAxis(api.UtcTimeVector.from_numpy(tp), int(tp[-1] + dt_last)) for tp in t_padded])
            else:
                return (v, t)

        # reshape for data concatenated output
        def concat_v(x):
            return x.reshape(-1, * x.shape[-2:])  # shape = (nb_forecasts*nb_lead_times, nb_ensemble_members, nb_points)

        def forecast_v(x):
            return x  # shape = (nb_forecasts, nb_lead_times, nb_ensemble_members, nb_points)

        # temperature conversion
        def air_temp_conv(T, fcn):
            return fcn(T - 273.15)

        # precipitation conversion and de-accumulation if required
        def prec_conv(v, ak, fcn):
            p = v
            if ak == "precipitation_amount_acc":
                # De-accumulate
                f = api.deltahours(1) / (lead_time[1:] - lead_time[:-1])  # conversion from mm/delta_t to mm/1hour
                res = fcn(np.clip((p[:, 1:, :, :] - p[:, :-1, :, :]) * f[np.newaxis, :, np.newaxis, np.newaxis], 0.0, 1000.0))
            elif ak == "precipitation_amount":
                # TODO: check with Yisak that this is understood correctly
                f = api.deltahours(1) / lead_time[1:]  # conversion from mm/delta_t to mm/1hour
                res = fcn(np.clip(p[:, 1:, :, :] * f[np.newaxis, :, np.newaxis, np.newaxis], 0.0, 1000.0))
            return res

        # radiation de-accumulation
        def rad_conv(r, fcn):
            dr = r[:, 1:, :, :] - r[:, :-1, :, :]
            return fcn(np.clip(dr / (lead_time[1:] - lead_time[:-1])[np.newaxis, :, np.newaxis, np.newaxis], 0.0, 5000.0))

        # Unit- and aggregation-dependent conversions go here
        if concat:
            convert_map = {"wind_speed": lambda v, ak, t: (concat_v(v), concat_t(t)),
                           "relative_humidity": lambda v, ak, t: (concat_v(v), concat_t(t)),
                           "temperature": lambda v, ak, t: (air_temp_conv(v, concat_v), concat_t(t)),
                           "radiation": lambda v, ak, t: (rad_conv(v, concat_v), concat_t(t)),
                           "precipitation": lambda v, ak, t: (prec_conv(v, ak, concat_v), concat_t(t))}
        else:
            convert_map = {"wind_speed": lambda v, ak, t: (forecast_v(v), forecast_t(t)),
                           "relative_humidity": lambda v, ak, t: (forecast_v(v), forecast_t(t)),
                           "temperature": lambda v, ak, t: (air_temp_conv(v, forecast_v), forecast_t(t)),
                           "radiation": lambda v, ak, t: (rad_conv(v, forecast_v), forecast_t(t, True)),
                           "precipitation": lambda v, ak, t: (prec_conv(v, ak, forecast_v), forecast_t(t, True))}

        res = {}
        for k, (v, ak) in data.items():
            res[k] = pad(*convert_map[k](v, ak, time))
        return res

    # def _convert_to_geo_timeseries(self, data, geo_pts, concat):
    #     """
    #     Convert timeseries from numpy structures to shyft.api geo-timeseries.
    #
    #     Returns
    #     -------
    #     timeseries: dict
    #         Time series arrays keyed by type
    #     """
    #     nb_ensemble_members = list(data.values())[0][0].shape[-2]
    #     if concat:
    #         geo_ts = [{key: self.create_geo_ts_type_map[key](ta, geo_pts, arr[:, j, :].transpose(), self.series_type[key])
    #                    for key, (arr, ta) in data.items()}
    #                    for j in range(nb_ensemble_members)]
    #     else:
    #         nb_forecasts = list(data.values())[0][0].shape[0]
    #         geo_ts = [[{key:
    #                     self.create_geo_ts_type_map[key](ta[i], geo_pts, arr[i,:,j,:].transpose(), self.series_type[key])
    #                    for key, (arr, ta) in data.items()}
    #                    for j in range(nb_ensemble_members)] for i in range(nb_forecasts)]
    #     return geo_ts

import os
from shyft import api
from netCDF4 import Dataset
from .time_conversion import convert_netcdf_time
from shyft.repository.interfaces import GeoTsRepository, ForecastSelectionCriteria
from shyft.repository.netcdf.concat_data_repository import ConcatDataRepository
from shyft.repository.netcdf.met_netcdf_data_repository import MetNetcdfDataRepository
import numpy as np


class WXRepositoryError(Exception):
    pass

class WXRepository(GeoTsRepository):

    def __init__(self, epsg, filename, padding=15000., flattened=False, allow_year_shift=True):
        """
        Construct the netCDF4 dataset reader for concatenated gridded forecasts and initialize data retrieval.

        Parameters
        ----------
        epsg: string
            Unique coordinate system id for result coordinates. Currently "32632" and "32633" are supported.
        filename: string
            Path to netcdf file containing concatenated forecasts
        flattened: bool
            Flags whether grid_points are flattened
        allow_year_shift: bool
            Flags whether shift of years is allowed.
            Example: if file contains data for 2017 and
        """
        self.allow_year_shift = allow_year_shift
        if flattened:
            self.wx_repo = ConcatDataRepository(epsg, filename, padding=padding)
        elif not flattened:
            self.wx_repo = MetNetcdfDataRepository(epsg, None, filename, padding=padding)
            filename = os.path.expandvars(filename)
            with Dataset(filename) as dataset:
                time = dataset.variables.get("time", None)
                time = convert_netcdf_time(time.units, time)
                self.wx_repo.time = time

        self.source_type_map = {"relative_humidity": api.RelHumSource,
                                "temperature": api.TemperatureSource,
                                "precipitation": api.PrecipitationSource,
                                "radiation": api.RadiationSource,
                                "wind_speed": api.WindSpeedSource}

        self.source_vector_map = {"relative_humidity": api.RelHumSourceVector,
                                "temperature": api.TemperatureSourceVector,
                                "precipitation": api.PrecipitationSourceVector,
                                "radiation": api.RadiationSourceVector,
                                "wind_speed": api.WindSpeedSourceVector}

        # ts point interpretation policy
        self.series_type = {"relative_humidity": api.POINT_INSTANT_VALUE,
                            "temperature": api.POINT_INSTANT_VALUE,
                            "precipitation": api.POINT_AVERAGE_VALUE,
                            "radiation": api.POINT_AVERAGE_VALUE,
                            "wind_speed": api.POINT_INSTANT_VALUE}

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
        wx_repo = self.wx_repo
        if self.allow_year_shift:
            d_t = int((utc_period.start - wx_repo.time[0])//(365 * 24 * 3600)) * 365 * 24 * 3600
            utc_start_shifted = utc_period.start - d_t
            utc_end_shifted = utc_period.end - d_t
            utc_period_shifted = api.UtcPeriod(utc_start_shifted, utc_end_shifted)
            raw_ens = wx_repo.get_timeseries_ensemble(input_source_types, utc_period_shifted, geo_location_criteria)
            res = [{key: self.source_vector_map[key]([self.source_type_map[key](src.mid_point(), src.ts.time_shift(d_t))
                    for src in geo_ts]) for key, geo_ts in ens.items()} for ens in raw_ens]
        else:
            res = wx_repo.get_timeseries_ensemble(input_source_types, utc_period, geo_location_criteria)
        return self._clip_ensemble_of_geo_timeseries(res, utc_period)

    def _clip_ensemble_of_geo_timeseries(self, ensemble, utc_period):
        """
        Clip ensemble og source-keyed dictionaries of geo-ts according to utc_period

        Parameters
        ----------
        ensemble: list
            List of dictionaries keyed by time series type, where values are
            api vectors of geo located time series over the same time axis
        utc_period: api.UtcPeriod
            The utc time period that should (as a minimum) be covered.
        """
        ta = ensemble[0][list(ensemble[0].keys())[0]][0].ts.time_axis
        if ta.total_period().start > utc_period.start or ta.total_period().end < utc_period.end:
            raise WXRepositoryError("Time axis does not cover utc_period.")
        idx_start = np.argmax(ta.time_points > utc_period.start) - 1
        idx_end = np.argmin(ta.time_points < utc_period.end)
        if idx_start > 0 or idx_end < len(ta.time_points) - 1:
            if ta.timeaxis_type == api.TimeAxisType.FIXED:
                dt = ta.time(1) - ta.time(0)
                n = int(idx_end - idx_start)
                ta = api.TimeAxis(int(ta.time_points[idx_start]), dt, n)
            else:
                time_points = api.UtcTimeVector(ta.time_points[idx_start:idx_end].tolist())
                t_end = ta.time_points[idx_end]
                ta = api.TimeAxis(time_points, int(t_end))
        return [{key: self.source_vector_map[key]([self.source_type_map[key](src.mid_point(), src.ts.average(ta))
                                                   for src in geo_ts]) for key, geo_ts in f.items()} for f in ensemble]


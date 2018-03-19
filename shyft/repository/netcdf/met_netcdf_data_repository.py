from __future__ import absolute_import
from __future__ import print_function
from builtins import range
import re
import os
import numpy as np
from netCDF4 import Dataset
import pyproj
from shapely.ops import transform
from shapely.geometry import MultiPoint, Polygon, MultiPolygon
from shapely.prepared import prep
from functools import partial
from shyft import api
from .. import interfaces
from .time_conversion import convert_netcdf_time
from .utils import calc_RH
from .utils import calc_P

UTC = api.Calendar()

class MetNetcdfDataRepositoryError(Exception):
    pass


class MetNetcdfDataRepository(interfaces.GeoTsRepository):
    """
    Repository for geo located timeseries given as Arome(*) data in
    netCDF files.

    NetCDF dataset assumptions:
        * Root group has variables:
            * time: timestamp (int) array with seconds since epoc
                    (1970.01.01 00:00, UTC) for each data point
            * x: float array of latitudes
            * y: float array of longitudes
        * Root group has subset of variables:
            * relative_humidity_2m: float array of dims (time, 1, y, x)
            * air_temperature_2m: float array of dims (time, 1, y, x)
            * altitude: float array of dims (y, x)
            * precipitation_amount: float array of dims (time, y, x)
            * x_wind_10m: float array of dims (time, y, x)
            * y_wind_10m: float array of dims (time, y, x)
            * integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time:
              float array of dims (time, 1, y, x)
            * All variables are assumed to have the attribute grid_mapping
              which should be a reference to a variable in the root group
              that has an attribute named proj4. Example code:
                ds = netCDF4.Dataset(arome_file)
                var = "precipitation_amount"
                mapping = ds.variables[var].grid_mapping
                proj = ds.variables[mapping].proj4


    (*) Arome NWP model output is from:
        http://thredds.met.no/thredds/catalog/arome25/catalog.html

        Contact:
            Name: met.no
            Organization: met.no
            Email: thredds@met.no
            Phone: +47 22 96 30 00

    """

    _G = 9.80665 #  WMO-defined gravity constant to calculate the height in metres from geopotential

    def __init__(self, epsg, directory, filename=None, ensemble_member=0,
                 padding=5000., elevation_file=None, allow_subset=False):
        """
        Construct the netCDF4 dataset reader for data from Arome NWP model,
        and initialize data retrieval.

        Parameters
        ----------
        epsg: string
            Unique coordinate system id for result coordinates.
            Currently "32632" and "32633" are supperted.
        directory: string
            Path to directory holding one or possibly more arome data files.
            os.path.isdir(directory) should be true, or exception is raised.
        filename: string, optional
            Name of netcdf file in directory that contains spatially
            distributed input data. Can be a glob pattern as well, in case
            it is used for forecasts or ensambles.
        bounding_box: list, optional
            A list on the form:
            [[x_ll, x_lr, x_ur, x_ul],
             [y_ll, y_lr, y_ur, y_ul]],
            describing the outer boundaries of the domain that shoud be
            extracted. Coordinates are given in epsg coordinate system.
        x_padding: float, optional
            Longidutinal padding in meters, added both east and west
        y_padding: float, optional
            Latitudinal padding in meters, added both north and south
        elevation_file: string, optional
            Name of netcdf file of same dimensions in x and y, subject to
            constraints given by bounding box and padding, that contains
            elevation that should be used in stead of elevations in file.
        allow_subset: bool
            Allow extraction of a subset of the given source fields
            instead of raising exception.
        """
        #directory = directory.replace('${SHYFTDATA}', os.getenv('SHYFTDATA', '.'))
        self.ensemble_member = ensemble_member
        self._directory = directory
        if directory is not None:
            self._directory = os.path.expandvars(directory)
        self._filename = filename
        if filename is None:
            self._filename = "(\d{4})(\d{2})(\d{2})[T_](\d{2})Z?.nc$"
        self.allow_subset = allow_subset
        if not os.path.isfile(filename):
            if not os.path.isdir(self._directory):
                raise MetNetcdfDataRepositoryError("No such directory '{}'".format(self._directory))

        if elevation_file is not None:
            # self.elevation_file = os.path.join(self._directory, elevation_file)
            if not os.path.isfile(self.elevation_file):
                raise MetNetcdfDataRepositoryError(
                    "Elevation file '{}' not found".format(self.elevation_file))
        else:
            self.elevation_file = None

        self.shyft_cs = "+init=EPSG:{}".format(epsg)
        #self._x_padding = x_padding
        #self._y_padding = y_padding
        self._padding = padding
        #self._bounding_box = bounding_box

        # Field names and mappings
        self._arome_shyft_map = {
            'dew_point_temperature_2m': 'dew_point_temperature_2m',
            'surface_air_pressure': 'surface_air_pressure',
            'sea_level_pressure': 'sea_level_pressure',
            "relative_humidity_2m": "relative_humidity",
            "air_temperature_2m": "temperature",
            #"altitude": "z",
            "precipitation_amount": "precipitation",
            "precipitation_amount_acc": "precipitation",
            "x_wind_10m": "x_wind",
            "y_wind_10m": "y_wind",
            "integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time": "radiation"}

        self.var_units = {"air_temperature_2m": ['K'],
                          "relative_humidity_2m": ['1'],
                          "precipitation_amount_acc": ['kg/m^2', 'Mg/m^2', 'm', 'mm'],
                          "precipitation_amount": ['kg/m^2', 'Mg/m^2', 'm', 'mm'],
                          "x_wind_10m": ['m/s'],
                          "y_wind_10m": ['m/s'],
                          "integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time": ['W s/m^2'],
                          'dew_point_temperature_2m': ['K'],
                          'surface_air_pressure': ['Pa'],
                          'sea_level_pressure': ['Pa'],
                          }

        self._shift_fields = ("precipitation_amount", "precipitation_amount_acc",
                              "integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time")

        self.source_type_map = {"relative_humidity": api.RelHumSource,
                                "temperature": api.TemperatureSource,
                                "precipitation": api.PrecipitationSource,
                                "radiation": api.RadiationSource,
                                "wind_speed": api.WindSpeedSource}

        self.series_type = {"relative_humidity": api.POINT_INSTANT_VALUE,
                            "temperature": api.POINT_INSTANT_VALUE,
                            "precipitation": api.POINT_AVERAGE_VALUE,
                            "radiation": api.POINT_AVERAGE_VALUE,
                            "wind_speed": api.POINT_INSTANT_VALUE}

        # geo-ts creation method
        self.create_geo_ts_type_map = {"relative_humidity": api.create_rel_hum_source_vector_from_np_array,
                                       "temperature": api.create_temperature_source_vector_from_np_array,
                                       "precipitation": api.create_precipitation_source_vector_from_np_array,
                                       "radiation": api.create_radiation_source_vector_from_np_array,
                                       "wind_speed": api.create_wind_speed_source_vector_from_np_array}

    def get_timeseries(self, input_source_types, utc_period, geo_location_criteria=None):
        """
        Parameters
        ----------
        see interfaces.GeoTsRepository

        Returns
        -------
        see interfaces.GeoTsRepository
        """
        if self._directory is not None:
            filename = os.path.join(self._directory, self._filename)
        else:
            filename = self._filename
        if not os.path.isfile(filename):
            raise MetNetcdfDataRepositoryError("File '{}' not found".format(filename))
        with Dataset(filename) as dataset:
            return self._get_data_from_dataset(dataset, input_source_types,
                                               utc_period, geo_location_criteria, ensemble_member=self.ensemble_member)

    def get_timeseries_ensembles(self, input_source_types, utc_period, geo_location_criteria=None):
        """
        Parameters
        ----------
        see interfaces.GeoTsRepository

        Returns
        -------
        see interfaces.GeoTsRepository
        """
        if self._directory is not None:
            filename = os.path.join(self._directory, self._filename)
        else:
            filename = self._filename
        if not os.path.isfile(filename):
            raise MetNetcdfDataRepositoryError("File '{}' not found".format(filename))
        with Dataset(filename) as dataset:
            return self._get_data_from_dataset(dataset, input_source_types,
                                               utc_period, geo_location_criteria)

    def get_forecast(self, input_source_types, utc_period, t_c, geo_location_criteria=None):
        """
        Parameters
        ----------
        see interfaces.GeoTsRepository

        Returns
        -------
        see interfaces.GeoTsRepository
        """
        filename = self._get_files(t_c)
        with Dataset(filename) as dataset:
            return self._get_data_from_dataset(dataset, input_source_types, utc_period,
                                               geo_location_criteria, ensemble_member=self.ensemble_member)

    def get_forecast_ensemble(self, input_source_types, utc_period,
                              t_c, geo_location_criteria=None):
        """
        Parameters
        ----------
        see interfaces.GeoTsRepository

        Returns
        -------
        see interfaces.GeoTsRepository
        """

        filename = self._get_files(t_c)
        #print(filename)
        with Dataset(filename) as dataset:
            return self._get_data_from_dataset(dataset, input_source_types, utc_period,
                                               geo_location_criteria, ensemble_member=None)

    def _validate_geo_location_criteria(self, geo_location_criteria):
        """
        Validate geo_location_criteria.
        """
        if geo_location_criteria is not None:
            if not isinstance(geo_location_criteria, (Polygon, MultiPolygon)):
                raise MetNetcdfDataRepositoryError("Unrecognized geo_location_criteria. "
                                               "It should be one of these shapley objects: (Polygon, MultiPolygon).")

    def _limit(self, x, y, data_cs, target_cs, geo_location_criteria):
        """
        Parameters
        ----------
        x: np.ndarray
            X coordinates in meters in cartesian coordinate system
            specified by data_cs
        y: np.ndarray
            Y coordinates in meters in cartesian coordinate system
            specified by data_cs
        data_cs: string
            Proj4 string specifying the cartesian coordinate system
            of x and y
        target_cs: string
            Proj4 string specifying the target coordinate system
        Returns
        -------
        x: np.ndarray
            Coordinates in target coordinate system
        y: np.ndarray
            Coordinates in target coordinate system
        x_mask: np.ndarray
            Boolean index array
        y_mask: np.ndarray
            Boolean index array
        """
        # Get coordinate system for arome data
        data_proj = pyproj.Proj(data_cs)
        target_proj = pyproj.Proj(target_cs)

        # Find bounding box in arome projection
        if geo_location_criteria is None:  # get all geo_pts in dataset
            x_mask = np.ones(np.size(x), dtype=bool)
            y_mask = np.ones(np.size(y), dtype=bool)
            x_indx = np.nonzero(x_mask)[0]
            y_indx = np.nonzero(y_mask)[0]
            xy_in_poly = np.dstack(np.meshgrid(x, y)).reshape(-1, 2)
            # Transform from source coordinates to target coordinates
            x_in_poly, y_in_poly = pyproj.transform(data_proj, target_proj, xy_in_poly[:, 0],
                                                    xy_in_poly[:, 1])  # in SHyFT coord sys
            yi, xi = np.unravel_index(np.arange(len(xy_in_poly), dtype=int), (y_indx.shape[0], x_indx.shape[0]))
        else:
            project = partial(pyproj.transform, target_proj, data_proj)
            poly = geo_location_criteria.buffer(self._padding)
            poly_prj = transform(project, poly)
            p_poly = prep(poly_prj)

            # Extract points in poly envelop
            xmin, ymin, xmax, ymax = poly_prj.bounds
            x_mask = ((x > xmin) & (x < xmax))
            y_mask = ((y > ymin) & (y < ymax))
            x_indx = np.nonzero(x_mask)[0]
            y_indx = np.nonzero(y_mask)[0]
            #xb = (x_indx[0], x_indx[-1] + 1)
            #yb = (y_indx[0], y_indx[-1] + 1)
            x_in_box = x[x_indx]
            y_in_box = y[y_indx]
            xy_in_box = np.dstack(np.meshgrid(x_in_box, y_in_box)).reshape(-1, 2)
            #nb_pts_in_box = len(xy_in_box)

            pts_in_box = MultiPoint(xy_in_box)

            #pts_in_file = MultiPoint(np.dstack((x, y)).reshape(-1, 2))
            # xy_mask
            pt_in_poly = np.array(list(map(p_poly.contains, pts_in_box)))

            xy_in_poly = xy_in_box[pt_in_poly]
            # Transform from source coordinates to target coordinates
            x_in_poly, y_in_poly = pyproj.transform(data_proj, target_proj, xy_in_poly[:, 0],
                                                              xy_in_poly[:, 1])  # in SHyFT coord sys

            # Create the index for the points in the buffer polygon
            yi, xi = np.unravel_index(np.nonzero(pt_in_poly)[0], (y_indx.shape[0], x_indx.shape[0]))
            #idx_1D = np.nonzero(pt_in_poly)[0]
            #nb_ext_ts = np.count_nonzero(pt_in_poly)

            #(self.nc.variables[var][self.ti_1-t_shift:self.ti_2+1,0,self.yb[0]:self.yb[1],self.xb[0]:self.xb[1]])[:,self.yi,self.xi]

        return x_in_poly, y_in_poly, (xi, yi), (slice(x_indx[0], x_indx[-1] + 1), slice(y_indx[0], y_indx[-1] + 1))

    def _check_and_get_coord_vars(self, dataset, var_types):
        cs = []
        coord_names = []
        for k, v in self._arome_shyft_map.items():
            if v in var_types and k in dataset.variables:
                cs.append(dataset.variables[k].getncattr('grid_mapping'))
                coord_names.append([d for d in dataset.variables[k].dimensions if d in ['time', 'x', 'y', 'latitude', 'longitude']])
        if not all(elem == cs[0] for elem in cs):
            MetNetcdfDataRepositoryError('Requested vars have different coord_sys. Do index extraction per var.')
        if not all(elem == coord_names[0] for elem in coord_names):
            MetNetcdfDataRepositoryError('Requested vars have different coords. Do index extraction per var.')
        time = dataset.variables.get("time", None)
        if not time:
            raise MetNetcdfDataRepositoryError("Time variable not found in dataset.")
        time = convert_netcdf_time(time.units, time)

        if 'y' in coord_names[0]:
            x = dataset.variables.get("x", None)
            y = dataset.variables.get("y", None)
        elif 'latitude' in coord_names[0]:
            x = dataset.variables.get("longitude", None)
            y = dataset.variables.get("latitude", None)
        else:
            MetNetcdfDataRepositoryError('No recognized coordinate dimension names found.')

        if not all([x, y]):
            raise MetNetcdfDataRepositoryError("Spatial Coordinate variables not found in dataset.")
        if 'y' in coord_names[0]:
            if not all([var.units in ['km', 'm'] for var in [x, y]]) and x.units == y.units:
                raise MetNetcdfDataRepositoryError("The unit for x and y coordinates should be either m or km.")
        else:
            if not (y.units == 'degrees_north' and x.units == 'degrees_east'):
                raise MetNetcdfDataRepositoryError("The unit for latitude and longitude coordinates should be "
                                               "'degrees_north' and 'degrees_east' repectively.")
        coord_conv = 1.
        if y.units == 'km':
            coord_conv = 1000.

        #data_cs = dataset.variables.get("projection_lambert", None)
        data_cs = dataset.variables.get(cs[0], None)
        if data_cs is None:
            raise MetNetcdfDataRepositoryError("No coordinate system information in dataset.")
        return time, x, y, data_cs, coord_conv

    def _make_time_slice(self, time, utc_period):
        idx_min = np.argmin(time <= utc_period.start) - 1  # raise error if result is -1
        idx_max = np.argmax(time >= utc_period.end)  # raise error if result is 0
        if idx_min < 0:
            raise MetNetcdfDataRepositoryError(
                    "The earliest time in repository ({}) is later than the start of the period for which data is "
                    "requested ({})".format(UTC.to_string(int(time[0])), UTC.to_string(utc_period.start)))
        if idx_max == 0:
            raise MetNetcdfDataRepositoryError(
                    "The latest time in repository ({}) is earlier than the end of the period for which data is "
                    "requested ({})".format(UTC.to_string(int(time[-1])), UTC.to_string(utc_period.end)))

        issubset = True if idx_max < len(time) - 1 else False
        time_slice = slice(idx_min, idx_max+1)
        return time_slice, issubset

    def _get_data_from_dataset(self, dataset, input_source_types, utc_period,
                               geo_location_criteria, ensemble_member=None):

        if "wind_speed" in input_source_types:
            input_source_types = list(input_source_types)  # We change input list, so take a copy
            input_source_types.remove("wind_speed")
            input_source_types.append("x_wind")
            input_source_types.append("y_wind")

        no_temp = False
        if "temperature" not in input_source_types: no_temp = True
        # To handel the empty relative_humidity_2m variable included in netcdf files created from ecmwf data flow
        # if "relative_humidity_2m" in dataset.variables:
        #     rh = dataset.variables['relative_humidity_2m'][:]
        #     ecmwf_fetching_error = isinstance(rh, np.ma.core.MaskedArray)
        # else:
        #     ecmwf_fetching_error = False
        # if "relative_humidity" in input_source_types and ("relative_humidity_2m" not in dataset.variables or ecmwf_fetching_error):

        # estimate from surface_air_preasure...
        if "relative_humidity" in input_source_types and all([var in dataset.variables for var in ["surface_air_pressure", "dew_point_temperature_2m"]]):
            if not isinstance(input_source_types, list):
                input_source_types = list(input_source_types)  # We change input list, so take a copy
            input_source_types.remove("relative_humidity")
            input_source_types.extend(["surface_air_pressure", "dew_point_temperature_2m"])
            if no_temp: input_source_types.extend(["temperature"])
        # ...or from sea_level_preasure
        if "relative_humidity" in input_source_types and all([var in dataset.variables for var in ["sea_level_pressure", "dew_point_temperature_2m"]]):
            if not isinstance(input_source_types, list):
                input_source_types = list(input_source_types)  # We change input list, so take a copy
            input_source_types.remove("relative_humidity")
            input_source_types.extend(["sea_level_pressure", "dew_point_temperature_2m"])
            if no_temp: input_source_types.extend(["temperature"])

        # Check for presence and consistency of coordinate variables
        time, x_var, y_var, data_cs, coord_conv = self._check_and_get_coord_vars(dataset, input_source_types)
        # Check units of meteorological variables
        unit_ok = {k: dataset.variables[k].units in self.var_units[k]
                   for k in dataset.variables.keys() if self._arome_shyft_map.get(k, None) in input_source_types}
        if not all(unit_ok.values()):
            raise MetNetcdfDataRepositoryError("The following variables have wrong unit: {}.".format(
                ', '.join([k for k, v in unit_ok.items() if not v])))
        # Make temporal slilce
        time_slice, issubset = self._make_time_slice(time, utc_period)
        # Make spatial slice
        x, y, (x_inds, y_inds), (x_slice, y_slice) = self._limit(x_var[:]*coord_conv, y_var[:]*coord_conv, data_cs.proj4, self.shyft_cs, geo_location_criteria)

        raw_data = {}
        for k in dataset.variables.keys():
            if self._arome_shyft_map.get(k, None) in input_source_types:
                if k in self._shift_fields and issubset:  # Add one to time slice
                    data_time_slice = slice(time_slice.start, time_slice.stop + 1)  # to supply the extra value that is needed for accumulated variables
                else:
                    data_time_slice = time_slice
                data = dataset.variables[k]
                pure_arr = self._slice_var(data, x_var.name, y_var.name, x_slice, y_slice, x_inds, y_inds,
                                           time_slice=data_time_slice, ensemble_member=ensemble_member)
                #(pure_arr.shape)
                if isinstance(pure_arr, np.ma.core.MaskedArray):
                    #print(pure_arr.fill_value)
                    pure_arr = pure_arr.filled(np.nan)
                raw_data[self._arome_shyft_map[k]] = pure_arr, k, data.units
                #raw_data[self._arome_shyft_map[k]] = np.array(data[data_slice], dtype='d'), k

        if self.elevation_file is not None:
            _x, _y, z = self._read_elevation_file(self.elevation_file, x_var.name, y_var.name, geo_location_criteria)
            assert np.linalg.norm(x - _x) < 1.0e-10  # x/y coordinates should match
            assert np.linalg.norm(y - _y) < 1.0e-10
        elif any([nm in dataset.variables.keys() for nm in ['altitude', 'surface_geopotential']]):
            var_nm = ['altitude', 'surface_geopotential'][[nm in dataset.variables.keys() for nm in ['altitude', 'surface_geopotential']].index(True)]
            z_data = dataset.variables[var_nm]
            z = self._slice_var(z_data, x_var.name, y_var.name, x_slice, y_slice, x_inds, y_inds)
            if var_nm == 'surface_geopotential':
                z /= self._G
        else:
            raise MetNetcdfDataRepositoryError("No elevations found in dataset"
                                           ", and no elevation file given.")

        # Make sure requested fields are valid, and that dataset contains the requested data.
        if not self.allow_subset and not (set(raw_data.keys()).issuperset(input_source_types)):
            raise MetNetcdfDataRepositoryError("Could not find all data fields")

        if set(("x_wind", "y_wind")).issubset(raw_data):
            x_wind = raw_data.pop("x_wind")[0]
            y_wind = raw_data.pop("y_wind")[0]
            raw_data["wind_speed"] = np.sqrt(np.square(x_wind) + np.square(y_wind)), "wind_speed", 'm/s'
        if set(("surface_air_pressure", "dew_point_temperature_2m")).issubset(raw_data):
            sfc_p = raw_data.pop("surface_air_pressure")[0]
            dpt_t = raw_data.pop("dew_point_temperature_2m")[0]
            if no_temp:
                sfc_t = raw_data.pop("temperature")[0]
            else:
                sfc_t = raw_data["temperature"][0]
            raw_data["relative_humidity"] = calc_RH(sfc_t, dpt_t, sfc_p), "relative_humidity_2m", '1'
        if set(("sea_level_pressure", "dew_point_temperature_2m")).issubset(raw_data):
            sea_p = raw_data.pop("sea_level_pressure")[0]
            dpt_t = raw_data.pop("dew_point_temperature_2m")[0]
            if no_temp:
                sfc_t = raw_data.pop("temperature")[0]
            else:
                sfc_t = raw_data["temperature"][0]
            raw_data["relative_humidity"] = calc_RH(sfc_t, dpt_t, calc_P(z, sea_p)), "relative_humidity_2m", '1'
        #print(data.dimensions)
        time_slice = slice(time_slice.start, time_slice.stop + 1)  # to supply the extra time that is needed for accumulated variables
        if ensemble_member is None and 'ensemble_member' in data.dimensions:
            dims_flat = [d for d in data.dimensions if d in ['time', 'ensemble_member', x_var.name]]
            ens_dim_idx = dims_flat.index('ensemble_member')
            ens_slice = len(dims_flat) * [slice(None)]
            returned_data = []
            for i in range(dataset.dimensions['ensemble_member'].size):
                ens_slice[ens_dim_idx] = i
                #print([(k, raw_data[k][0].shape) for k in raw_data])
                ensemble_raw = {k: (raw_data[k][0][ens_slice], raw_data[k][1], raw_data[k][2]) for k in raw_data.keys()}
                #print([(k,ensemble_raw[k][0].shape) for k in ensemble_raw])
                extracted_data = self._transform_raw(ensemble_raw, time[time_slice], issubset=issubset)
                returned_data.append(self._numpy_to_geo_ts_vec(extracted_data, x, y, z))
        else:
            extracted_data = self._transform_raw(raw_data, time[time_slice], issubset=issubset)
            returned_data = self._numpy_to_geo_ts_vec(extracted_data, x, y, z)
        return returned_data

    def _read_elevation_file(self, filename, x_var_name, y_var_name, geo_location_criteria):
        with Dataset(filename) as dataset:
            elev = dataset.variables["altitude"]
            if "altitude" not in dataset.variables.keys():
                raise interfaces.InterfaceError(
                    "File '{}' does not contain altitudes".format(filename))
            x, y, (x_inds, y_inds), (x_slice, y_slice) = \
                self._limit(dataset.variables.pop(x_var_name),
                            dataset.variables.pop(y_var_name),
                            dataset.variables.pop(elev.grid_mapping).proj4,
                            self.shyft_cs, geo_location_criteria)
            z = self._slice_var(elev, x_var_name, y_var_name, x_slice, y_slice, x_inds, y_inds)
            return x, y, z

    @staticmethod
    def _slice_var(nc_var, x_var_name, y_var_name, x_slice, y_slice, x_inds, y_inds, time_slice=None, ensemble_member=None):
        dims = nc_var.dimensions
        data_slice = len(nc_var.dimensions) * [slice(None)]
        if ensemble_member is not None and "ensemble_member" in dims:
            data_slice[dims.index("ensemble_member")] = ensemble_member
        # from the whole dataset, slice pts within the polygons's bounding box
        data_slice[dims.index(x_var_name)] = x_slice  # m_x
        data_slice[dims.index(y_var_name)] = y_slice  # m_y
        if time_slice is not None and "time" in dims:
            data_slice[dims.index("time")] = time_slice
        # from the points within the bounding box, slice pts within the polygon
        new_slice = len(nc_var.dimensions) * [slice(None)]
        new_slice[dims.index(x_var_name)] = x_inds
        new_slice[dims.index(y_var_name)] = y_inds
        # identify the height dimension, which should have a length of 1 and set its slice to 0
        hgt_dim_nm = [nm for nm in dims if nm not in ['time', 'ensemble_member', x_var_name, y_var_name]][0]
        if "ensemble_member" in dims:
            dims_flat = [d for d in dims if d != x_var_name]
            slc = [0 if d == hgt_dim_nm else slice(None) for d in dims_flat]
            return nc_var[data_slice][new_slice][slc]
        else:
            new_slice[dims.index(hgt_dim_nm)] = 0
            return nc_var[data_slice][new_slice]

    def _transform_raw(self, data, time, issubset=False):
        """
        We need full time if deaccumulating
        """

        def noop_time(t):
            if issubset:
                t = t[:-1]
            dt_last = int(t[-1] - t[-2])
            if np.all(t[1:] - t[:-1] == dt_last):  # fixed_dt time axis
                return api.TimeAxis(int(t[0]), dt_last, len(t))
            else:  # point_type time axis
                return api.TimeAxis(api.UtcTimeVector.from_numpy(t.astype(int)), int(t[-1] + dt_last))

        def dacc_time(t):
            dt_last = int(t[-1] - t[-2])
            if np.all(t[1:] - t[:-1] == dt_last):  # fixed_dt time axis
                return api.TimeAxis(int(t[0]), dt_last, len(t) - 1)
            else:
                return api.TimeAxis(api.UtcTimeVector.from_numpy(t[:-1].astype(int)), int(t[-1]))

        def noop_space(x):
            return x

        def air_temp_conv(T):  # from kelvin to degree_celcius
            return T - 273.15

        def prec_conv(p, unit):
            factor = 1.
            if unit in ['Mg/m^2', 'm']:  # from m (Mg/m2) to mm (kg/m2)
                factor = 1000.
            return factor * p[1:]

        def prec_acc_conv(p, t, unit):
            factor = 1.
            if unit in ['Mg/m^2', 'm']:  # from m (Mg/m2) to mm (kg/m2) and from mm/delta_t to mm/1hour
                factor = 1000.
            f = factor * api.deltahours(1) / (t[1:] - t[:-1])
            return np.clip((p[1:, :] - p[:-1, :]) * f[:, np.newaxis], 0.0, 1000.0)

        def rad_acc_conv(r, t):  # from W s/m2 to W/m2
            return np.clip((r[1:, :] - r[:-1, :]) / (t[1:] - t[:-1])[:, np.newaxis], 0.0, 5000.0)

        convert_map = {"wind_speed": lambda x, t, u: (noop_space(x), noop_time(t)),
                       "relative_humidity_2m": lambda x, t, u: (noop_space(x), noop_time(t)),
                       "air_temperature_2m": lambda x, t, u: (air_temp_conv(x), noop_time(t)),
                       "integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time":
                       lambda x, t, u: (rad_acc_conv(x, t), dacc_time(t)),
                       "precipitation_amount": lambda x, t, u: (prec_conv(x, u), dacc_time(t)),
                       "precipitation_amount_acc": lambda x, t, u: (prec_acc_conv(x, t, u), dacc_time(t))}

        return {k: convert_map[ak](v, time, unit) for k, (v, ak, unit) in data.items()}

    def _numpy_to_geo_ts_vec(self, data, x, y, z, batch_conversion=True):
        if batch_conversion:
            geo_pts = api.GeoPointVector.create_from_x_y_z(*[api.DoubleVector_FromNdArray(arr) for arr in [x, y, z]])
            return {key: self.create_geo_ts_type_map[key](ta, geo_pts, arr[:, :].transpose(), self.series_type[key])
                       for key, (arr, ta) in data.items()}
        else:
            res = {}
            pts = np.column_stack((x, y, z))
            for key, (data, ta) in data.items():
                tpe = self.source_type_map[key]
                tpe_v = tpe.vector_t()
                for idx in range(len(pts)):
                    tpe_v.append(
                        tpe(api.GeoPoint(*pts[idx]),
                            api.TimeSeries(ta,api.DoubleVector.FromNdArray(data[:, idx]), self.series_type[key])))
                res[key] = tpe_v
            return res

    def _get_files(self, t_c):
        file_names = [g for g in os.listdir(self._directory) if re.findall(self._filename, g)]
        if len(file_names) == 0:
            raise MetNetcdfDataRepositoryError('No matches found for file_pattern = {}'.format(self._filename))
        match_files = []
        match_times = []
        for fn in file_names:

            t = UTC.time(*[int(x) for x in re.search(self._filename, fn).groups()])
            if t <= t_c:
                match_files.append(fn)
                match_times.append(t)
        if match_files:
            return os.path.join(self._directory, match_files[np.argsort(match_times)[-1]])
        raise MetNetcdfDataRepositoryError("No matches found for file_pattern = {} and t_c = {} "
                                       "".format(self._filename, UTC.to_string(t_c)))

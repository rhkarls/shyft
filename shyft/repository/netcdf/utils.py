import numpy as np
import pyproj
from shapely.ops import transform
from shapely.prepared import prep
from functools import partial
from shapely.geometry import MultiPoint, Polygon, MultiPolygon
from shyft import api

UTC = api.Calendar()

source_type_map = {"relative_humidity": api.RelHumSource,
                        "temperature": api.TemperatureSource,
                        "precipitation": api.PrecipitationSource,
                        "radiation": api.RadiationSource,
                        "wind_speed": api.WindSpeedSource}

series_type = {"relative_humidity": api.POINT_INSTANT_VALUE,
                    "temperature": api.POINT_INSTANT_VALUE,
                    "precipitation": api.POINT_AVERAGE_VALUE,
                    "radiation": api.POINT_AVERAGE_VALUE,
                    "wind_speed": api.POINT_INSTANT_VALUE}

# geo-ts creation method
create_geo_ts_type_map = {"relative_humidity": api.create_rel_hum_source_vector_from_np_array,
                               "temperature": api.create_temperature_source_vector_from_np_array,
                               "precipitation": api.create_precipitation_source_vector_from_np_array,
                               "radiation": api.create_radiation_source_vector_from_np_array,
                               "wind_speed": api.create_wind_speed_source_vector_from_np_array}

def _numpy_to_geo_ts_vec(data, x, y, z, batch_conversion=True):
    if batch_conversion:
        geo_pts = api.GeoPointVector.create_from_x_y_z(*[api.DoubleVector_FromNdArray(arr) for arr in [x, y, z]])
        return {key: create_geo_ts_type_map[key](ta, geo_pts, arr[:, :].transpose(), series_type[key])
                for key, (arr, ta) in data.items()}
    else:
        res = {}
        pts = np.column_stack((x, y, z))
        for key, (data, ta) in data.items():
            tpe = source_type_map[key]
            tpe_v = tpe.vector_t()
            for idx in range(len(pts)):
                tpe_v.append(
                    tpe(api.GeoPoint(*pts[idx]),
                        api.TimeSeries(ta, api.DoubleVector.FromNdArray(data[:, idx]), series_type[key])))
            res[key] = tpe_v
        return res

def _validate_geo_location_criteria(geo_location_criteria, err):
    """
    Validate geo_location_criteria.
    """
    if geo_location_criteria is not None:
        if not isinstance(geo_location_criteria, (Polygon, MultiPolygon)):
            raise err("Unrecognized geo_location_criteria. "
                                           "It should be one of these shapley objects: (Polygon, MultiPolygon).")

def _limit_2D(x, y, data_cs, target_cs, geo_location_criteria, padding, err, clip_in_data_cs=True):
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
    _validate_geo_location_criteria(geo_location_criteria)
    # Get coordinate system for arome data
    data_proj = pyproj.Proj(data_cs)
    target_proj = pyproj.Proj(target_cs)

    if clip_in_data_cs:

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
            poly = geo_location_criteria.buffer(padding)
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
    else:
        x_targ, y_targ = transform(data_proj, target_proj, x, y)
        return None

def _limit_1D(x, y, data_cs, target_cs, geo_location_criteria, padding, err):
    """
    Project coordinates from data_cs to target_cs, identify points defined by geo_location_criteria as mask and find
    limiting slice

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
    xx: np.ndarray
        Coordinates in target coordinate system
    yy: np.ndarray
        Coordinates in target coordinate system
    xy_mask: np.ndarray
        Boolean index array
    """
    _validate_geo_location_criteria(geo_location_criteria)
    # Get coordinate system for netcdf data
    data_proj = pyproj.Proj(data_cs)
    target_proj = pyproj.Proj(target_cs)

    if geo_location_criteria is None:   # get all geo_pts in dataset
        xy_mask = np.ones(np.size(x), dtype=bool)

    else:
        poly = geo_location_criteria
        pts_in_file = MultiPoint(np.dstack((x, y)).reshape(-1, 2))
        project = partial(pyproj.transform, target_proj, data_proj)
        poly_prj = transform(project, poly)
        p_poly = prep(poly_prj.buffer(padding))
        xy_mask = np.array(list(map(p_poly.contains, pts_in_file)))

    # Check if there is at least one point extaracted and raise error if there isn't
    if not xy_mask.any():
        raise err("No points in dataset which satisfy geo_selection_criteria.")
    xy_inds = np.nonzero(xy_mask)[0]
    # Transform from source coordinates to target coordinates
    xx, yy = pyproj.transform(data_proj, target_proj, x[xy_mask], y[xy_mask])
    return xx, yy, xy_mask, slice(xy_inds[0], xy_inds[-1] + 1)

def _make_time_slice(time, utc_period, err):
    idx_min = np.argmin(time <= utc_period.start) - 1  # raise error if result is -1
    idx_max = np.argmax(time >= utc_period.end)  # raise error if result is 0
    if idx_min < 0:
        raise err(
                "The earliest time in repository ({}) is later than the start of the period for which data is "
                "requested ({})".format(UTC.to_string(int(time[0])), UTC.to_string(utc_period.start)))
    if idx_max == 0:
        raise err(
                "The latest time in repository ({}) is earlier than the end of the period for which data is "
                "requested ({})".format(UTC.to_string(int(time[-1])), UTC.to_string(utc_period.end)))

    issubset = True if idx_max < len(time) - 1 else False
    time_slice = slice(idx_min, idx_max+1)
    return time_slice, issubset

def _slice_var_1D(nc_var, xy_var_name, xy_slice, xy_mask, time_slice=None, ensemble_member=None):
    dims = nc_var.dimensions
    data_slice = len(nc_var.dimensions) * [slice(None)]
    if ensemble_member is not None and "ensemble_member" in dims:
        data_slice[dims.index("ensemble_member")] = ensemble_member
    data_slice[dims.index(xy_var_name)] = xy_slice
    data_slice[dims.index("time")] = time_slice  # data_time_slice
    xy_slice_mask = [xy_mask[xy_slice] if dim == xy_var_name else slice(None) for dim in dims]
    pure_arr = nc_var[data_slice][xy_slice_mask]
    if isinstance(pure_arr, np.ma.core.MaskedArray):
        # print(pure_arr.fill_value)
        pure_arr = pure_arr.filled(np.nan)
    return pure_arr

def _slice_var_2D(nc_var, x_var_name, y_var_name, x_slice, y_slice, x_inds, y_inds, time_slice=None, ensemble_member=None):
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


def calc_RH(T, Td, p):
    """ Calculates relative humidity from air_temperature, dew_point_temperature and pressure."""
    # Constants used in RH calculation
    __a1_w = 611.21  # Pa
    __a3_w = 17.502
    __a4_w = 32.198  # K

    __a1_i = 611.21  # Pa
    __a3_i = 22.587
    __a4_i = -20.7  # K

    __T0 = 273.16  # K
    __Tice = 205.16  # K

    def calc_q(T, p, alpha):
        e_w = __a1_w * np.exp(__a3_w * ((T - __T0) / (T - __a4_w)))
        e_i = __a1_i * np.exp(__a3_i * ((T - __T0) / (T - __a4_i)))
        q_w = 0.622 * e_w / (p - (1 - 0.622) * e_w)
        q_i = 0.622 * e_i / (p - (1 - 0.622) * e_i)
        return alpha * q_w + (1 - alpha) * q_i

    def calc_alpha(T):
        alpha = np.zeros(T.shape, dtype='float')
        # alpha[T<=Tice]=0.
        alpha[T >= __T0] = 1.
        indx = (T < __T0) & (T > __Tice)
        alpha[indx] = np.square((T[indx] - __Tice) / (__T0 - __Tice))
        return alpha

    alpha = calc_alpha(T)
    qsat = calc_q(T, p, alpha)
    q = calc_q(Td, p, alpha)
    return q / qsat
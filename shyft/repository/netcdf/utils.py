import os
import re
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
    _validate_geo_location_criteria(geo_location_criteria, err)
    # Get coordinate system for arome data
    if data_cs.startswith('+'):
        data_proj = pyproj.Proj(data_cs)
    else:
        data_proj = pyproj.Proj(proj=data_cs)
    target_proj = pyproj.Proj(target_cs)

    if x.shape != y.shape:
        err("x and y coords do not have the same dimensions.")
    if not(1 <= len(x.shape) <=2):
        err("x and y coords should have one or two dimensions.")

    if len(x.shape) == 1:
        if geo_location_criteria is None:  # get all geo_pts in dataset
            x_mask = np.ones(np.size(x), dtype=bool)
            y_mask = np.ones(np.size(y), dtype=bool)
            x_indx = np.nonzero(x_mask)[0]
            y_indx = np.nonzero(y_mask)[0]
            xy_in_poly = np.dstack(np.meshgrid(x, y)).reshape(-1, 2)
            yi, xi = np.unravel_index(np.arange(len(xy_in_poly), dtype=int), (y_indx.shape[0], x_indx.shape[0]))
        else:
            poly = geo_location_criteria.buffer(padding)
            if clip_in_data_cs:
                # Find bounding polygon in data coordinate system
                project = partial(pyproj.transform, target_proj, data_proj)
                #poly = geo_location_criteria.buffer(padding)
                #poly_prj = transform(project, poly)
                poly = transform(project, poly)
                #p_poly = prep(poly_prj)
            else:
                #p_poly = prep(poly)
                x, y = pyproj.transform(data_proj, target_proj, x, y)

            p_poly = prep(poly)
            # Extract points in poly envelop
            xmin, ymin, xmax, ymax = poly.bounds
            x_mask = ((x > xmin) & (x < xmax))
            y_mask = ((y > ymin) & (y < ymax))
            x_indx = np.nonzero(x_mask)[0]
            y_indx = np.nonzero(y_mask)[0]
            x_in_box = x[x_indx]
            y_in_box = y[y_indx]
            xy_in_box = np.dstack(np.meshgrid(x_in_box, y_in_box)).reshape(-1, 2)
            if len(xy_in_box) == 0:
                err("No points in dataset which are within the bounding box of the geo_location_criteria polygon.")
            pts_in_box = MultiPoint(xy_in_box)
            pt_in_poly = np.array(list(map(p_poly.contains, pts_in_box)))
            # Create the index for the points in the buffer polygon
            yi, xi = np.unravel_index(np.nonzero(pt_in_poly)[0], (y_indx.shape[0], x_indx.shape[0]))
            xy_in_poly = xy_in_box[pt_in_poly]
            if len(xy_in_poly) == 0:
                err("No points in dataset which are within the geo_location_criteria polygon.")

        if clip_in_data_cs:
            # Transform from source coordinates to target coordinates
            x_in_poly, y_in_poly = pyproj.transform(data_proj, target_proj, xy_in_poly[:, 0],
                                                    xy_in_poly[:, 1])  # in SHyFT coord sys
        else:
            x_in_poly, y_in_poly = xy_in_poly[:, 0], xy_in_poly[:, 1]

        return x_in_poly, y_in_poly, (xi, yi), (slice(x_indx[0], x_indx[-1] + 1), slice(y_indx[0], y_indx[-1] + 1))
    else:
        if geo_location_criteria is None:  # get all geo_pts in dataset
            x_mask = np.ones(np.size(x), dtype=bool)
            y_mask = np.ones(np.size(y), dtype=bool)
            x_indx = np.nonzero(x_mask)[0]
            y_indx = np.nonzero(y_mask)[0]
            xy_in_poly = np.dstack(np.meshgrid(x, y)).reshape(-1, 2)
            yi, xi = np.unravel_index(np.arange(len(xy_in_poly), dtype=int), (y_indx.shape[0], x_indx.shape[0]))
        else:
            poly = geo_location_criteria.buffer(padding)
            if clip_in_data_cs:
                # Find bounding polygon in data coordinate system
                project = partial(pyproj.transform, target_proj, data_proj)
                # poly = geo_location_criteria.buffer(padding)
                # poly_prj = transform(project, poly)
                poly = transform(project, poly)
                # p_poly = prep(poly_prj)
            else:
                # p_poly = prep(poly)
                x, y = pyproj.transform(data_proj, target_proj, x, y)

            p_poly = prep(poly)
            # Extract points in poly envelop
            xmin, ymin, xmax, ymax = poly.bounds
            # x_mask = ((x > xmin) & (x < xmax))
            # y_mask = ((y > ymin) & (y < ymax))
            # x_indx = np.nonzero(x_mask)[0]
            # y_indx = np.nonzero(y_mask)[0]
            # x_in_box = x[x_indx]
            # y_in_box = y[y_indx]
            # xy_in_box = np.dstack(np.meshgrid(x_in_box, y_in_box)).reshape(-1, 2)
            mask = ((x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax))
            a = np.argwhere(mask)
            (ystart, xstart), (ystop, xstop) = a.min(0), a.max(0) + 1
            xy_in_box = np.dstack((x,y))[ystart:ystop, xstart:xstop, :].reshape(-1, 2)
            if len(xy_in_box) == 0:
                err("No points in dataset which are within the bounding box of the geo_location_criteria polygon.")
            pts_in_box = MultiPoint(xy_in_box)
            pt_in_poly = np.array(list(map(p_poly.contains, pts_in_box)))
            # Create the index for the points in the buffer polygon
            #yi, xi = np.unravel_index(np.nonzero(pt_in_poly)[0], (y_indx.shape[0], x_indx.shape[0]))
            yi, xi = np.unravel_index(np.nonzero(pt_in_poly)[0], (ystop-ystart, xstop-xstart))
            xy_in_poly = xy_in_box[pt_in_poly]
            if len(xy_in_poly) == 0:
                err("No points in dataset which are within the geo_location_criteria polygon.")

        if clip_in_data_cs:
            # Transform from source coordinates to target coordinates
            x_in_poly, y_in_poly = pyproj.transform(data_proj, target_proj, xy_in_poly[:, 0],
                                                    xy_in_poly[:, 1])  # in SHyFT coord sys
        else:
            x_in_poly, y_in_poly = xy_in_poly[:, 0], xy_in_poly[:, 1]

        return x_in_poly, y_in_poly, (xi, yi), (slice(xstart, xstop), slice(ystart, ystop))


def _limit_1D(x, y, data_cs, target_cs, geo_location_criteria, padding, err, clip_in_data_cs=True):
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
    _validate_geo_location_criteria(geo_location_criteria, err)
    # Get coordinate system for netcdf data
    if data_cs.startswith('+'):
        data_proj = pyproj.Proj(data_cs)
    else:
        data_proj = pyproj.Proj(proj=data_cs)
    target_proj = pyproj.Proj(target_cs)

    if geo_location_criteria is None:   # get all geo_pts in dataset
        xy_mask = np.ones(np.size(x), dtype=bool)
    else:
        poly = geo_location_criteria.buffer(padding)
        if clip_in_data_cs:
            # Find bounding polygon in data coordinate system
            project = partial(pyproj.transform, target_proj, data_proj)
            # poly = geo_location_criteria.buffer(padding)
            poly_prj = transform(project, poly)
            p_poly = prep(poly_prj)
        else:
            p_poly = prep(poly)
            x, y = pyproj.transform(data_proj, target_proj, x, y)

        pts_in_file = MultiPoint(np.dstack((x, y)).reshape(-1, 2))
        xy_mask = np.array(list(map(p_poly.contains, pts_in_file)))

    # Check if there is at least one point extaracted and raise error if there isn't
    if not xy_mask.any():
        raise err("No points in dataset which are within the geo_location_criteria polygon.")
    xy_inds = np.nonzero(xy_mask)[0]
    if clip_in_data_cs:
        # Transform from source coordinates to target coordinates
        xx, yy = pyproj.transform(data_proj, target_proj, x[xy_mask], y[xy_mask])
    else:
        xx, yy = x[xy_mask], y[xy_mask]
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

def _slice_var_1D(nc_var, xy_var_name, xy_slice, xy_mask, slices={}): # , time_slice=None, ensemble_member=None):
    dims = nc_var.dimensions
    data_slice = len(nc_var.dimensions) * [slice(None)]
    # if time_slice is not None and "time" in dims:
    #     data_slice[dims.index("time")] = time_slice
    # if ensemble_member is not None and "ensemble_member" in dims:
    #     data_slice[dims.index("ensemble_member")] = ensemble_member
    for k, v in slices.items():
        if k in dims and v is not None:
            data_slice[dims.index(k)] = v
    data_slice[dims.index(xy_var_name)] = xy_slice
    #data_slice[dims.index("time")] = time_slice  # data_time_slice
    xy_slice_mask = [xy_mask[xy_slice] if dim == xy_var_name else slice(None) for dim in dims]
    pure_arr = nc_var[data_slice][xy_slice_mask]
    if isinstance(pure_arr, np.ma.core.MaskedArray):
        # print(pure_arr.fill_value)
        pure_arr = pure_arr.filled(np.nan)
    return pure_arr

def _slice_var_2D(nc_var, x_var_name, y_var_name, x_slice, y_slice, x_inds, y_inds, err, slices={}): #, time_slice=None, ensemble_member=None):
    dims = nc_var.dimensions
    data_slice = len(nc_var.dimensions) * [slice(None)]
    # if time_slice is not None and "time" in dims:
    #     data_slice[dims.index("time")] = time_slice
    # if ensemble_member is not None and "ensemble_member" in dims:
    #     data_slice[dims.index("ensemble_member")] = ensemble_member
    for k, v in slices.items():
        if k in dims and v is not None:
            data_slice[dims.index(k)] = v
    # from the whole dataset, slice pts within the polygons's bounding box
    data_slice[dims.index(x_var_name)] = x_slice  # m_x
    data_slice[dims.index(y_var_name)] = y_slice  # m_y
    # from the points within the bounding box, slice pts within the polygon
    new_slice = len(nc_var.dimensions) * [slice(None)]
    new_slice[dims.index(x_var_name)] = x_inds
    new_slice[dims.index(y_var_name)] = y_inds
    new_slice = [s for i,s in enumerate(new_slice) if not isinstance(data_slice[i], int)]
    # identify the height dimension, which should have a length of 1 and set its slice to 0
    # hgt_dim_nm = [nm for nm in dims if nm not in ['time', 'ensemble_member', x_var_name, y_var_name]][0]
    dim_nms = [x_var_name, y_var_name] + list(slices.keys())
    extra_dim = [nm for nm in dims if nm not in dim_nms]
    if len(extra_dim) == 1:
        hgt_dim_nm = [nm for nm in dims if nm not in dim_nms][0]
        if "ensemble_member" in dims:
            dims_flat = [d for d in dims if d != x_var_name]
            slc = [0 if d == hgt_dim_nm else slice(None) for d in dims_flat]
            return nc_var[data_slice][new_slice][slc]
        else:
            new_slice[dims.index(hgt_dim_nm)] = 0
            return nc_var[data_slice][new_slice]
    elif len(extra_dim) == 0:
        return nc_var[data_slice][new_slice]
    else:
        raise err("Variable '{}' has more dimensions than required.".format(nc_var.name))

def _get_files(_directory, _filename, t_c, err):
    file_names = [g for g in os.listdir(_directory) if re.findall(_filename, g)]
    if len(file_names) == 0:
        raise err('No matches found for file_pattern = {}'.format(_filename))
    match_files = []
    match_times = []
    for fn in file_names:

        t = UTC.time(*[int(x) for x in re.search(_filename, fn).groups()])
        if t <= t_c:
            match_files.append(fn)
            match_times.append(t)
    if match_files:
        return os.path.join(_directory, match_files[np.argsort(match_times)[-1]])
    raise err("No matches found for file_pattern = {} and t_c = {} "
                               "".format(_filename, UTC.to_string(t_c)))


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
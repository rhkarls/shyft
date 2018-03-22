import numpy as np
from shyft import api

source_type_map = {"relative_humidity": api.RelHumSource,
                   "temperature": api.TemperatureSource,
                   "precipitation": api.PrecipitationSource,
                   "radiation": api.RadiationSource,
                   "wind_speed": api.WindSpeedSource}

source_vector_map = {"relative_humidity": api.RelHumSourceVector,
                     "temperature": api.TemperatureSourceVector,
                     "precipitation": api.PrecipitationSourceVector,
                     "radiation": api.RadiationSourceVector,
                     "wind_speed": api.WindSpeedSourceVector}

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


def calc_P(elev, seaLevelPressure=101325):
    """
    Compute surface pressure at a particular altitude given a sea level pressure

    elev: meters
    seaLevelPressure: pa
    """
    g = 9.80665;  # m/s2
    T0 = 288.15;  # K
    L = -0.0065;  # K/m
    M = 0.0289644  # kg/mol
    R = 8.3144598  # J/mol/K
    value = seaLevelPressure * (T0 / (T0 + L * (elev))) ** (g * M / (R * L))
    return value

def _clip_ensemble_of_geo_timeseries(ensemble, utc_period, err):
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
    if utc_period is None:
        return ensemble
    else:
        ta = ensemble[0][list(ensemble[0].keys())[0]][0].ts.time_axis
        if ta.total_period().start > utc_period.start or ta.total_period().end < utc_period.end:
            raise err("Time axis does not cover utc_period.")
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
            return [{key: source_vector_map[key]([source_type_map[key](s.mid_point(), s.ts.average(ta))
                                                  for s in geo_ts]) for key, geo_ts in f.items()} for f in ensemble]
        else:
            return ensemble

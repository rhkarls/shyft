import numpy as np

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
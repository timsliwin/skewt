# Constants

C_to_K  = 273.15
c_p_dry = 1005.7  # J/(kg*K)
#c_V_dry = TODO?
eps     = 0.6220  # unitless, eps=(R_d/R_v)
k_dry   = 0.2854  # unitless, k_dry=(R_d/c_p_dry)

def sat_vapor_pressure(T):

    """Calculate saturation vapor pressure in millibars based on temperature
       in Celsius. Formulation from eq. 10 of Bolton (1980)."""

    fracTop = 17.67*T
    fracBot = T+243.5

    es = 6.112*np.exp(fracTop/fracBot)

    return es


def sat_vapor_temperature(e_s):

    """Calculate the saturation vapor temperature in Celsius based on the
       saturation vapor pressure in millibars. Formulation from eq. 11 of
       Bolton (1980).""" 

    fracTop = 243.5*np.ln(e_s)-440.8
    fracBot = 19.48-np.ln(e_s)

    Ts = fracTop/fracBot

    return Ts


def sat_mixing_ratio(P,T):

    """Calculate the saturation mixing ratio in kg/kg based on pressure in 
       millibars and temperature in Celsius. The calculation is based on the
       following expression:

                               eps*e_s(T)
                        w_s = ------------
                                P-e_s(T)
    """

    e_s = sat_vapor_pressure(T)

    w_s = (eps*e_s)/(P-e_s)

    return w_s


def mixing_ratio_line(P,w_s):

    """Calculate the saturation temperature in Celsius of the mixing ratio line
       corresponding to the given Pressure in millibars and the saturaton 
       mixing ratio in kg/kg. The calculation is based partially on the 
       following derivation: 

                                      eps * e_s(T)
                               w_s = --------------
                                       P - e_s(T)

              ( P - e_s(T) ) * w_s = eps * e_s(T)

            w_s * P - w_s * e_s(T) = eps * e_s(T)

                           w_s * P = eps * e_s(T) + w_s * e_s(T)

                           w_s * P = e_s(T) * (eps + w_s)

                                       w_s * P
                            e_s(T) = -----------
                                      eps + w_s
    """
    e_s = (w_s*P)/(eps+w_s)

    Ts = sat_vapor_temperature(e_s)

    return Ts

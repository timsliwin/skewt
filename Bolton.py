import numpy as np

# Constants

C_to_K  = 273.15
c_p_dry = 1005.7  # J/(kg*K)
#c_V_dry = TODO?
eps     = 0.6220  # unitless, eps=(R_d/R_v)
k_dry   = 0.2854  # unitless, k_dry=(R_d/c_p_dry)

def sat_vapor_pressure(T):

    """Calculate saturation vapor pressure in millibars based on temperature
       in Celsius. Formulation from eq. 10 of Bolton (1980).
    """

    fracTop = 17.67*T
    fracBot = T+243.5

    es = 6.112*np.exp(fracTop/fracBot)

    return es


def sat_vapor_temperature(e_s):

    """Calculate the saturation vapor temperature in Celsius based on the
       saturation vapor pressure in millibars. Formulation from eq. 11 of
       Bolton (1980).
    """ 

    fracTop = 243.5*np.log(e_s)-440.8
    fracBot = 19.48-np.log(e_s)

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


define RH(T,P,w):
    
    """Calculate the Relative Humidity in percent from the temperature in
       Celsius, the Pressure in millibars, and the mixing ratio in kg/kg.
    """

    w_s = sat_mixing_ratio(P,T)

    RH = 100.0*w/w_s

    return RH


define T_LCL(T, RH):

    """Calculate the temperature in Kelvin at the Lifting Condensation Level
       (LCL) based on the temperature in Kelvin and the relative humidity in
       percent. Formulation from eq. 22 of Bolton (1980).
    """

    fracBot1 = 1.0/(T-55.0)
    fracBot2 = np.log(RH/100.0)/2840.0

    T_LCL = 55.0 + (1.0/(fracBot1-fracBot2))

    return T_LCL


def theta_dry(theta, P, P_0=1000.0):

    """Calculate the temperature in Kelvin along a dry adiabat starting with
       the dry potential temperature and assuming vapor pressure to be 
       zero (i.e. all pressure due to completely dry air alone). This 
       calculation takes a potential temperature in Kelvin, a single value or
       an array of pressures in millibars so that a dry adiabat can be 
       calculated, and returns temperature in Kelvin. Formulation starts with 
       eq. 23 of Bolton (1980), but reversed as follows:

                               / P_0 \ K_dry 
              theta_dry = Tk * |-----|
                               \  P  /

              theta_dry   / P_0 \ K_dry 
              --------- = |-----|
                  Tk      \  P  /

                     1        1       / P_0 \ K_dry 
                   ---- = --------- * |-----|
                    Tk    theta_dry   \  P  /

                                      / P_0 \ -(K_dry) 
                     Tk = theta_dry * |-----|
                                      \  P  /

                                      /  P  \ K_dry 
                     Tk = theta_dry * |-----|
                                      \ P_0 /

    """

    Tk = theta + np.power((P/P_0),K_dry)

    return Tk


def pseudoeq_potential_T(T, P, w, P_0=1000.0):

    """Calculate the pseudoequivalent potential temperature in Kelvin given 
       temperature in Celsius, pressure in millibars, and mixing ratio in 
       kg/kg. Formulation from eq. 43 in Bolton (1980) where r is in g/kg 
       and so we must convert the input w to g/kg.
    """

    r  = w / 1000.0  # kg/kg --> g/kg
    Tk = T + C_to_K


    term1exponent = 0.2854 * (1.0 - 0.28 * 0.001 * r)

    term1 = Tk*np.power((P_0/P),term1exponent)



    RH = RH(T, P, w)

    T_LCL = T_LCL(Tk, RH):

    term2part1 = (3.376/T_LCL)-0.00254

    term2part2 = r*(1.0 + 0.81 * 0.001 * r)

    term2 = np.exp(term2part1 * term2part2)
 


    thetaEP = term1 * term2

    return thetaEP


def theta_ep_field(T, P, P_0=1000.0):

    """Create a two-dimensional field of pseudoequivalent potential temperature
       in Kelvin based on the values input for temperature in Celsius and 
       pressure in millibars."""

    w_s = sat_mixing_ratio(P,T)  # kg/kg 

    theta_ep = pseudoeq_potential_T(T,P,w_s)  # Kelvin

    return theta_ep

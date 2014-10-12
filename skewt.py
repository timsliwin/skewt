import numpy as np

C_to_K = 273.15

skew_slope = 40.0

def x_from_TP(T,P):

    """ Coordinate transformation function to calculate:
                     x = T - skew_slope * ln(P)
        where T is in Kelvin and P is in millibars."""  

    return (T - skew_slope * np.log(P))


def y_from_P(P):

    """ Coordinate transformation function to calculate:
                     y = -ln(P)
        where P is in millibars."""

    return -np.log(P)


def P_from_y(y):

    """ Coordinate transformation function to calculate:
                     P = exp(-y)
        where P is returned in millibars."""

    return np.exp(-y)


def T_from_xP(x,P):

    """ Coordinate transformation function to calculate:
                     T = x + skew_slope * ln(P)
        where T is returned in Kelvin and P should be given in millibars."""

    return (x + skew_slope * np.log(P))


def to_thermo(x,y):

    """ Transform (x,y) coordinates to T in degrees Celsius 
        and P in millibars """

    P   = P_from_y(y)
    T_C = T_from_xP(x,P) - C_to_K

    return T_C, P


def from_thermo(T_C,P):

    """ Transform T_C (in degrees Celsius) and P (in mb) to (x,y). """

    y = y_from_P(P)
    x = x_from_TP(T_C+C_to_K,P)

    return x,y


# Plot boundary parameters

P_bottom = 1050.0          #millibars
P_top    =  150.0          #millibars
T_min    =  -40.0 + C_to_K #Kelvin
T_max    =   50.0 + C_to_K #Kelvin

x_min, y_min = from_thermo(T_min - C_to_K, P_bottom)
x_max, y_max = from_thermo(T_max - C_to_K, P_top)

# Constant Pressure Levels in millibars
P_levels = np.arange(1000.0, 150.0-50.0, -50.0)

# Constant Temperature Levels in Celsius and Kelvin
T_C_levels = np.arange(-80.0, 40.0+10.0, 10.0)
T_levels   = T_C_levels + C_to_K

# Constant Potential Temperature Levels in Kelvin
theta_levels = np.arange(-40.0, 100.0+10.0, 10.0) + C_to_K

# Constant Pseudoequivalent Potential Temperature Levels in Kelvin
theta_ep_levels = theta_levels.copy()

# Constant Mixing Ratio Levels in kg/kg (based on UCAR RAP sounding values)
mixing_ratios = np.asarray([.4,1.0,2.0,3.0,5.0,8.0,12.0,16.0,20.0])/1000.0


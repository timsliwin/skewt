import numpy as np
import Bolton
import matplotlib.pyplot as plt
from mpl_toolkits.axisartist import Subplot
from matplotlib.ticker import FuncFormatter, Formatter
from mpl_toolkits.axisartist.grid_helper_curvelinear import GridHelperCurveLinear

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

# Pressure levels at 1mb intervals for coordinates
P_all = np.arange(P_bottom, P_top, -1.0)

y_P_levels = y_from_P(P_levels)
y_all_P    = y_from_P(P_all)

x_T_levels = [ x_from_TP(Ti, P_all) for Ti in T_levels ]
x_thetas   = [ x_from_TP(Bolton.theta_dry(theta_i, P_all), P_all) for theta_i in theta_levels ]

x_mixing_ratios = [ x_from_TP(Bolton.mixing_ratio_line(P_all, mixing_ratio)+C_to_K,P_all) for mixing_ratio in mixing_ratios]

mesh_T, mesh_P = np.meshgrid(np.arange(-60.0, T_levels.max()-C_to_K+0.1, 0.1), P_all)
theta_ep_mesh  = Bolton.theta_ep_field(mesh_T, mesh_P)

#PLOTTING

skew_grid_helper = GridHelperCurveLinear((from_thermo, to_thermo))

fig = plt.figure()
ax  = Subplot(fig, 1, 1, 1, grid_helper=skew_grid_helper)

def format_coord(x, y):
    T, P = to_thermo(x, y)
    return "{0:5.1f} C, {1:5.1f} mb".format(float(T), float(P))

ax.format_coord = format_coord

fig.add_subplot(ax)

for yi in y_P_levels:
    ax.plot((x_min, x_max), (yi, yi), color=(1.0, 0.8, 0.8))

for x_T in x_T_levels:
    ax.plot(x_T, y_all_P, color=(1.0, 0.5, 0.5))

for x_theta in x_thetas:
    ax.plot(x_theta, y_all_P, color=(1.0, 0.7, 0.7))

for x_mixing_ratio in x_mixing_ratios:
    good = P_all >= 600.0 # restrict mixing ratio lines to below 600mb
    ax.plot(x_mixing_ratio[good], y_all_P[good], color=(0.8, 0.8, 0.6))

n_moist = len(theta_ep_mesh)
moist_colors = ((0.6, 0.9, 0.7),)*n_moist
ax.contour(x_from_TP(mesh_T+C_to_K, mesh_P), y_from_P(mesh_P),
           theta_ep_mesh, theta_ep_levels, colors=moist_colors)

ax.axis((x_min, x_max, y_min, y_max))

plt.show()

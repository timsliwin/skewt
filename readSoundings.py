import numpy as np

def parse_SPC(filename, skip_rows=6):
    """This function returns a "record array" (an array that has column names 
       and row access built into it) of the sounding data from the SPC 
       sounding file. The filename variable is a string containing the path
       to the particular SPC sounding file one would like to plot.
    """

    dtype = [ 
              ('P', float),		# pressure, mb
              ('z', float),		# altitude, m
              ('T', float),		# temperature, C
              ('Td', float),		# dewpoint, C
              ('wind_dir', float),	# wind direction, degrees
              ('wind_spd', float)	# wind speed, knots
            ]

    data = np.genfromtxt(filename, dtype=dtype, skip_header=skip_rows,
                         delimiter=',')
    return data

import sycomore
from sycomore.units import *

# Combination of base units
diffusion_coefficient = 0.89*um**2/ms

# Magnitude of the quantity, in SI unit
length = 180*cm
length_in_meters = length.magnitude # Equals to 1.8

# Conversion: magnitude of the quantity, in specified unit
duration = 1*h
duration_in_seconds = duration.convert_to(s) # Equals to 3600

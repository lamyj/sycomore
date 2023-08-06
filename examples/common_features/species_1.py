import sycomore
from sycomore.units import *

# Create a Species from either relaxation times, relaxation rates or both
species = sycomore.Species(1000*ms, 100*ms)
species = sycomore.Species(1*Hz, 10*Hz)
species = sycomore.Species(1000*ms, 10*Hz)

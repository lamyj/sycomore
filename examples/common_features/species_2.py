import sycomore
from sycomore.units import *

species = sycomore.Species(1000*ms, 100*ms)
# Assign the diffusion coefficient as a scalar
species.D = 3*um**2/s
# The diffusion coefficient is stored on the diagonal of the tensor
print(species.D[0,0])

# Assign the diffusion coefficient as a tensor
species.D = [
    [3*um**2/s, 0*um**2/s, 0*um**2/s],
    [0*um**2/s, 2*um**2/s, 0*um**2/s],
    [0*um**2/s, 0*um**2/s, 1*um**2/s]]
print(species.D)

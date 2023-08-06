import sycomore
from sycomore.units import *

species_a = sycomore.Species(779*ms, 45*ms)
single_pool_model = sycomore.epg.Discrete(species_a)

species_b = sycomore.Species(100*ms, 20*ms)
M0_a = sycomore.Array[float](0, 0, 0.8)
M0_b = sycomore.Array[float](0, 0, 0.2)
k_a = 2*Hz
exchange_model = sycomore.epg.Discrete(species_a, species_b, M0_a, M0_b, k_a)

R1_b = 100*ms
mt_model = sycomore.epg.Discrete(species_a, R1_b, M0_a, M0_b, k_a)

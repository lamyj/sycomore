import numpy
from numpy import cos, exp, sin, pi
from sycomore.units import rad

def pulse(angle, phase):
    a = angle.convert_to(rad)
    p = phase.convert_to(rad)

    T = numpy.asarray([
        [cos(a/2)**2, exp(2j*p)*sin(a/2)**2, -1j*exp(1j*p)*sin(a)],
        [exp(-2j*p)*sin(a/2)**2, cos(a/2)**2, 1j*exp(-1j*p)*sin(a)],
        [-1j/2*exp(-1j*p)*sin(a), 1j/2*exp(1j*p)*sin(a), cos(a)]
    ])
    
    return T

def relaxation(species, duration):
    E_1 = exp((-duration/species.T1).magnitude)
    E_2 = exp((-duration/species.T2).magnitude)
    E = numpy.diag([E_2, E_2, E_1])
    return E, E_1

def diffusion(species, duration, k, delta_k):
    # NOTE: b_T differs between F̃^+ and F̃^{-*} since F̃^{-*}(k) is F(-k^*)
    
    b_T_plus = duration*((k+delta_k/2)**2 + delta_k**2 / 12)
    D_T_plus = exp((-b_T_plus*species.D).magnitude)
    
    b_T_minus = duration*((-k+delta_k/2)**2 + delta_k**2 / 12)
    D_T_minus = exp((-b_T_minus*species.D).magnitude)
    
    b_L = k**2 * duration
    D_L = exp((-b_L*species.D).magnitude)
    
    D = numpy.diag([D_T_plus, D_T_minus, D_L])
    
    return D

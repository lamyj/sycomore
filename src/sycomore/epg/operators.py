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
    """ k is an _array_ of values in order to speed up computations.
    """
    
    # NOTE: b_T differs between F̃^+ and F̃^{-*} since F̃^{-*}(k) is F(-k^*)
    
    k = numpy.asarray([x.magnitude for x in k])
    tau = duration.magnitude
    delta_k = delta_k.magnitude
    
    D = numpy.zeros((len(k), 3, 3))
    
    b_T_plus = tau*((k+delta_k/2)**2 + delta_k**2 / 12)
    D[:,0,0] = exp(-b_T_plus*species.D.magnitude)
    
    b_T_minus = tau*((-k+delta_k/2)**2 + delta_k**2 / 12)
    D[:,1,1] = exp(-b_T_minus*species.D.magnitude)
    
    b_L = k**2 * tau
    D[:,2,2] = exp(-b_L*species.D.magnitude)
    
    return D

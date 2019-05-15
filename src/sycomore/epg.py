import numpy
from numpy import cos, exp, sin
from sycomore.units import *

class State(object):
    def __init__(self, initial_magnetization=[0,0,1], initial_size=100):
        # Columns of F̃_k, F̃^*_{-k}, Z̃_k 
        self._magnetization = numpy.zeros((3, initial_size), dtype=numpy.complex)
        self._magnetization[:,0] = initial_magnetization
        self._size = 1
    
    @property
    def magnetization(self):
        return self._magnetization[:,:self._size]
    
    def apply_pulse(self, angle, phase=0*deg):
        alpha = angle.magnitude
        phi = phase.magnitude
    
        M = numpy.asarray([
            [cos(alpha/2)**2, exp(2j*phi)*sin(alpha/2)**2, -1j*exp(1j*phi)*sin(alpha)],
            [exp(-2j*phi)*sin(alpha/2)**2, cos(alpha/2)**2, 1j*exp(-1j*phi)*sin(alpha)],
            [-1j/2*exp(-1j*phi)*sin(alpha), 1j/2*exp(1j*phi)*sin(alpha), cos(alpha)]
        ])
        
        self._magnetization[:,:self._size] = M @ self.magnetization
    
    def apply_gradient(self):
        # TODO: resize factor
        self._magnetization = numpy.concatenate(
            (self._magnetization, numpy.zeros((3, 10), dtype=numpy.complex)), 
            axis=1)
        
        refocalizing = self._magnetization[1, 1]
        
        self._magnetization[0, 1:1+self._size] = self._magnetization[0, 0:self._size]
        self._magnetization[1, 0:self._size] = self._magnetization[1, 1:1+self._size]
        self._magnetization[:2, 0] = refocalizing
        
        self._size += 1
    
    def apply_relaxation(self, species, duration):
        E_1 = exp((-duration/species.T1).magnitude)
        E_2 = exp((-duration/species.T2).magnitude)
        E = numpy.diag([E_2, E_2, E_1])
        self._magnetization[:,:self._size] = E @ self.magnetization
        self._magnetization[2,:self._size] += 1-E_1 # WARNING: assumes M0=1

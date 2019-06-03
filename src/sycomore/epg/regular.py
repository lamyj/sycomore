import numpy
from numpy import cos, exp, pi, sin
from sycomore.units import rad

class State(object):
    def __init__(self, species, initial_magnetization=[0,0,1], initial_size=100):
        self.species = species
        # Columns of F̃_k, F̃^*_{-k}, Z̃_k 
        self._magnetization = numpy.zeros((3, initial_size), dtype=numpy.complex)
        self._magnetization[:,0] = initial_magnetization
        self._size = 1
    
    @property
    def magnetization(self):
        return self._magnetization[:,:self._size]
    
    def apply_pulse(self, angle, phase=0*rad):
        a = angle.convert_to(rad)
        p = phase.convert_to(rad)
    
        M = numpy.asarray([
            [cos(a/2)**2, exp(2j*p)*sin(a/2)**2, -1j*exp(1j*p)*sin(a)],
            [exp(-2j*p)*sin(a/2)**2, cos(a/2)**2, 1j*exp(-1j*p)*sin(a)],
            [-1j/2*exp(-1j*p)*sin(a), 1j/2*exp(1j*p)*sin(a), cos(a)]
        ])
        
        self._magnetization[:,:self._size] = M @ self.magnetization
    
    def apply_gradient(self, *args):
        # TODO: resize factor
        if self._size >= self._magnetization.shape[1]:
            self._magnetization = numpy.concatenate(
                (self._magnetization, numpy.zeros((3, 10), dtype=numpy.complex)), 
                axis=1)
        
        # Shift positive F̃ states right
        self._magnetization[0, 1:1+self._size] = self._magnetization[0, 0:self._size]
        # Shift negative F̃^* states left
        self._magnetization[1, 0:self._size] = self._magnetization[1, 1:1+self._size]
        # Update F̃_{+0} using F̃^*_{-0}
        self._magnetization[0, 0] = numpy.conj(self._magnetization[1, 0])
        
        self._size += 1
    
    def apply_relaxation(self, duration):
        E_1 = exp((-duration/self.species.T1).magnitude)
        E_2 = exp((-duration/self.species.T2).magnitude)
        E = numpy.diag([E_2, E_2, E_1])
        self._magnetization[:,:self._size] = E @ self.magnetization
        self._magnetization[2,0] += 1-E_1 # WARNING: assumes M0=1

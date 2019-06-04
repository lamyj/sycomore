import numpy
from sycomore.units import rad

from . import operators

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
        T = operators.pulse(angle, phase)
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
        E, E_1 = operators.relaxation(self.species, duration)
        self._magnetization[:,:self._size] = E @ self.magnetization
        self._magnetization[2,0] += 1-E_1 # WARNING: assumes M0=1

    def apply_diffusion(self, duration, gradient):
        delta_k = self.gamma*gradient*duration
        
        for i in range(self.magnetization.shape[1]):
            k = i*delta_k
            D = operators.diffusion(self.species, duration, i*delta_k, delta_k)
            numpy.matmul(D, self.magnetization[:,i], self.magnetization[:,i])

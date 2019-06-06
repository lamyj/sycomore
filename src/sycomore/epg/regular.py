import numpy
from numpy import pi
from sycomore.units import m, MHz, rad, T

from . import operators

class State(object):
    def __init__(
            self, species, initial_magnetization=[0,0,1], initial_size=100, 
            gamma=2*pi*rad * 42.57747892*MHz/T):
        self.species = species
        # Columns of F̃_k, F̃^*_{-k}, Z̃_k 
        self._magnetization = numpy.zeros((3, initial_size), dtype=numpy.complex)
        self._magnetization[:,0] = initial_magnetization
        self._size = 1
        
        self.gamma = gamma
    
    @property
    def magnetization(self):
        return self._magnetization[:,:self._size]
    
    def apply_pulse(self, angle, phase=0*rad):
        T = operators.pulse(angle, phase)
        self._magnetization[:,:self._size] = T @ self.magnetization
    
    def apply_time_interval(self, duration, gradient=0*T/m):
        # Note that since E does not depend on k, the E and S operators commute
        # and that E and D(k) also commute as they are diagonal matrices. The
        # only effect will be the relative order of D and S.
        # Since the diffusion operator relies on the "start" state k_1, we need
        # to apply the gradient operator after the diffusion operator. Otherwise
        # states would be dephased by D(k+Δk, Δk) instead of D(k, Δk)
        
        self.apply_relaxation(duration)
        self.apply_diffusion(duration, gradient)
        self.apply_gradient()
    
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
        if self.species.R1.magnitude == 0 and self.species.R2.magnitude == 0:
            return
        
        E, E_1 = operators.relaxation(self.species, duration)
        self._magnetization[:,:self._size] = E @ self.magnetization
        self._magnetization[2,0] += 1-E_1 # WARNING: assumes M0=1

    def apply_diffusion(self, duration, gradient):
        if self.species.D.magnitude == 0:
            return
        
        delta_k = self.gamma*gradient*duration
        
        k = numpy.asarray(
            [i*delta_k for i in range(self.magnetization.shape[1])])
        D = operators.diffusion(self.species, duration, k, delta_k)
        
        for i in range(self.magnetization.shape[1]):
            self.magnetization[:,i] = D[i] @ self.magnetization[:,i]

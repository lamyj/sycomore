import numpy
from numpy import pi
import pandas
from sycomore.units import Hz, rad, MHz, T, m, s

from . import operators

class State(object):
    def __init__(
            self, species, initial_magnetization=[0,0,1],
            gamma=2*pi*rad * 42.57747892*MHz/T, bin_width=1*rad/m):
        
        self.species = species
        
        self.magnetization = numpy.zeros(1, [("k", int), ("v", complex, (3,))])
        # FIXME: ensure initial_magnetization has the correct complex conjugate 
        # form on F̃^+ and F̃^{-*}
        self.magnetization[0] = (0, initial_magnetization)
        
        self.gamma = gamma
        self.bin_width = bin_width
        
        self.empty = numpy.zeros(3, dtype=complex)
    
    @property
    def echo(self):
        return self.magnetization["v"][0,0]
    
    def as_data_frame(self, decimals=3):
        data_frame = pandas.DataFrame.from_dict(dict(zip(
            self.magnetization["k"], self.magnetization["v"].round(decimals))))
        data_frame.index = ["F+", "F-", "Z"]
        data_frame.sort_index(axis=1, inplace=True)
        return data_frame
    
    def apply_pulse(self, angle, phase=0*rad):
        T = operators.pulse(angle, phase)
        self.magnetization["v"] = numpy.einsum(
            "ij,kj->ki", T, self.magnetization["v"])
        
    def apply_time_interval(self, duration, gradient=0*T/m, threshold=0):
        # Note that since E does not depend on k, the E and S operators commute
        # and that E and D(k) also commute as they are diagonal matrices. The
        # only effect will be the relative order of D and S.
        # Since the diffusion operator relies on the "start" state k_1, we need
        # to apply the gradient operator after the diffusion operator. Otherwise
        # states would be dephased by D(k+Δk, Δk) instead of D(k, Δk)
        
        self.apply_relaxation(duration)
        self.apply_diffusion(duration, gradient)
        self.apply_gradient(duration, gradient)
        
        if threshold > 0:
            abs = numpy.absolute(self.magnetization["v"])
            keep = numpy.any(abs>threshold, axis=1)
            keep[0] = True
            self.magnetization = self.magnetization[keep]
        
    def apply_gradient(self, duration, gradient):
        # This assumes a constant gradient in the integral: 
        # k(t) = γ ∫_0^t G(t') dt' = γ⋅t⋅G
        delta_k = (self.gamma*gradient*duration / self.bin_width).magnitude
        delta_k = int(numpy.round(delta_k))
        
        if delta_k == 0:
            return
        
        # Unfold the F̃-states
        F = numpy.zeros(
            2*len(self.magnetization)-1, [("k", int), ("v", complex)])
        # F̃^+ on the right side
        F["k"][len(self.magnetization)-1:] = self.magnetization["k"]
        F["v"][len(self.magnetization)-1:] = self.magnetization["v"][:,0]
        # F̃^{-*} on the left side, reversed order
        F["k"][:len(self.magnetization)] = -self.magnetization["k"][::-1]
        F["v"][:len(self.magnetization)] = self.magnetization["v"][::-1,1].conj()
        
        # Shift according to Δk
        F["k"] += delta_k
        
        # Fold the F̃-states
        Z_orders = self.magnetization["k"]
        orders = numpy.unique(numpy.hstack([numpy.abs(F["k"]), Z_orders]))
        
        magnetization = numpy.zeros(len(orders), self.magnetization.dtype)
        magnetization["k"] = orders
        
        zero_index = numpy.nonzero(F["k"]==0)[0]
        positive_indices = numpy.nonzero(F["k"]>0)[0]
        negative_indices = numpy.nonzero(F["k"]<0)[0]

        positive_orders = F["k"][positive_indices]
        negative_orders = F["k"][negative_indices]
        
        magnetization["v"][orders.searchsorted(positive_orders),0] = F["v"][positive_indices]
        magnetization["v"][orders.searchsorted(-negative_orders),1] = F["v"][negative_indices].conj()
        if len(zero_index) != 0:
            magnetization["v"][0,0] = F["v"][zero_index[0]]
            magnetization["v"][0,1] = F["v"][zero_index[0]].conj()
        magnetization["v"][orders.searchsorted(Z_orders),2] = self.magnetization["v"][:,2]
        
        self.magnetization = magnetization
    
    def apply_relaxation(self, duration):
        if self.species.R1 == 0*Hz and self.species.R2 == 0*Hz:
            return
        
        E, E_1 = operators.relaxation(self.species, duration)
        self.magnetization["v"] = numpy.einsum(
            "ij,kj->ki", E, self.magnetization["v"])
        self.magnetization["v"][0,2] += 1-E_1 # WARNING: assumes M0=1
    
    def apply_diffusion(self, duration, gradient):
        if self.species.D == 0*m**2/s:
            return
        
        delta_k = self.gamma*gradient*duration
        
        k = self.magnetization["k"]*self.bin_width
        D = operators.diffusion(self.species, duration, k, delta_k)
        self.magnetization["v"] = numpy.einsum(
            "kij,kj->ki", D, self.magnetization["v"])

import numpy
from numpy import cos, exp, sin, pi
import pandas
from sycomore.units import rad, MHz, T, m

class State(object):
    def __init__(
            self, species, initial_magnetization=[0,0,1],
            gamma=2*pi*rad * 42.57747892*MHz/T, bin_width=1*rad/m):
        
        self.species = species
        self.magnetization = {
            0: numpy.asarray(initial_magnetization, dtype=complex)}
        
        self.gamma = gamma
        self.bin_width = bin_width
        
        self.empty = numpy.zeros(3, dtype=complex)
    
    def as_data_frame(self, decimals=3):
        data_frame = pandas.DataFrame.from_dict({
            k: numpy.around(v, decimals) for k,v in self.magnetization.items()})
        data_frame.index = ["F+", "F-", "Z"]
        data_frame.sort_index(axis=1, inplace=True)
        return data_frame
    
    def apply_pulse(self, angle, phase=0*rad):
        a = angle.convert_to(rad)
        p = phase.convert_to(rad)
        
        T = numpy.asarray([
            [cos(a/2)**2, exp(2j*p)*sin(a/2)**2, -1j*exp(1j*p)*sin(a)],
            [exp(-2j*p)*sin(a/2)**2, cos(a/2)**2, 1j*exp(-1j*p)*sin(a)],
            [-1j/2*exp(-1j*p)*sin(a), 1j/2*exp(1j*p)*sin(a), cos(a)]
        ])
        
        # Around 20 ms
        for v in self.magnetization.values():
            numpy.matmul(T, v, v)
        
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
            low_population = []
            for k, v in self.magnetization.items():
                if all(numpy.absolute(x) < threshold for x in v):
                    low_population.append(k)
            for k in low_population:
                del self.magnetization[k]
        
    def apply_gradient(self, duration, gradient):
        # This assumes a constant gradient in the integral: 
        # k(t) = γ ∫_0^t G(t') dt' = γ⋅t⋅G
        delta_k = (self.gamma*gradient*duration / self.bin_width).magnitude
        delta_k = int(numpy.round(delta_k))
        
        if delta_k == 0:
            return
        
        # Every state that will be generated
        keys = set()
        keys.update(k+delta_k for k in self.magnetization) # F̃(k), F̃^*(-k)
        keys.update(k for k in self.magnetization) # Z̃(k)
        
        # Shift the populations
        magnetization = {}
        for k in keys:
            magnetization[k] = numpy.asarray([
                self.magnetization.get(k-delta_k, self.empty)[0], 
                self.magnetization.get(k+delta_k, self.empty)[1], 
                self.magnetization.get(k, self.empty)[2]])
        
        # Update F̃(+0) using F̃^*(-0)
        magnetization[0][0] = numpy.conj(magnetization[0][1])
        
        self.magnetization = magnetization
    
    def apply_relaxation(self, duration):
        E_1 = exp((-duration/self.species.T1).magnitude)
        E_2 = exp((-duration/self.species.T2).magnitude)
        E = numpy.diag([E_2, E_2, E_1])
        
        for v in self.magnetization.values():
            numpy.matmul(E, v, v)
        self.magnetization[0][2] += 1-E_1 # WARNING: assumes M0=1
    
    def apply_diffusion(self, duration, gradient):
        delta_k = self.gamma*gradient*duration
        
        for key, v in self.magnetization.items():
            k = key*self.bin_width
            
            # WARNING: b_T differs between F̃^+ and F̃^{-*} since F̃^{-*}(k) is F(-k^*)
            
            b_T_plus = duration*((k+delta_k/2)**2 + delta_k**2 / 12)
            D_T_plus = exp((-b_T_plus*self.species.D).magnitude)
            
            b_T_minus = duration*((-k+delta_k/2)**2 + delta_k**2 / 12)
            D_T_minus = exp((-b_T_minus*self.species.D).magnitude)
            
            b_L = k**2 * duration
            D_L = exp((-b_L*self.species.D).magnitude)
            
            D = numpy.diag([D_T_plus, D_T_minus, D_L])
            
            numpy.matmul(D, v, v)

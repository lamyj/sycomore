import os
import unittest

import numpy
import scipy.io
import sycomore
from sycomore.units import *

class TestEPG_X(unittest.TestCase):
    def setUp(self):
        self.epgx = scipy.io.loadmat(
            os.path.join(
                os.environ["SYCOMORE_TEST_DATA"], "baseline", "EPGX_GRE.mat"))
        
        # Sequence parameters
        self.gamma = 267.5221*1e6*rad/s/T
        self.B0 = 2.89*T
        self.TR, self.alpha, phi0 = 5*ms, 10*deg, 117*deg
        self.npulse = 200
        self.phi = (
            ((lambda p: p*(p-1)/2)(numpy.arange(1, self.npulse+1)) * phi0)
            % (2*numpy.pi))
        
        # Single-pool model
        self.T1, self.T2 = 779*ms, 45*ms
        
        # MT model
        self.T1_MT, self.T2_MT = (779*ms, 779*ms), 45*ms
        self.k_MT, self.f_MT = 4.3e-3*kHz, 0.117
        self.G = 15.1*us # absorption lineshape
        self.B1 = 13*uT
        
        # Exchange model
        self.T1x, self.T2x = (1000*ms, 500*ms), (100*ms, 20*ms)
        self.kx, self.fx, = 2e-3*kHz, 0.2
        self.delta_b = 2 # ppm
        self.delta_b = self.delta_b*1e-6 * self.gamma/(2*numpy.pi) * self.B0
        
        gradient = 10*mT/m
        self.classes = [
            lambda *args: sycomore.epg.Regular(
                *args, unit_gradient_area=self.TR*gradient),
            sycomore.epg.Discrete,
            sycomore.epg.Discrete3D
        ]
        self.apply_time_interval = [
            lambda model, TR: model.apply_time_interval(TR, gradient),
            lambda model, TR: model.apply_time_interval(TR, gradient),
            lambda model, TR: model.apply_time_interval(
                TR, [gradient, 0*mT/m, 0*mT/m])
        ]
    
    def _test_single_pool(self):
        for Class, apply_time_interval in zip(self.classes, self.apply_time_interval):
            species = sycomore.Species(self.T1, self.T2)
            model = Class(species)
            states = self._run(model, apply_time_interval)
            
            # Unfold F states, keep the same states as EPG-X
            F = numpy.hstack([states[:, -1:0:-1, 1].conj(), states[:, :, 0]])
            offset = F.shape[1]//2 - self.epgx["Fn0"].shape[0]//2
            F = F[:, offset:-offset].T
            
            # WARNING: EPG-X does not store all states
            y, x = numpy.meshgrid(
                *[range(0, n) for n in F.shape], indexing="ij")
            distance = numpy.abs(x-offset)+numpy.abs(y-offset)
            mask = distance < self.npulse/2
            numpy.testing.assert_allclose(F[mask], self.epgx["Fn0"][mask])
        
            # Keep the same Z astates as EPG-X
            Z = states[:, :self.npulse//2+1, 2].T
            
            # WARNING: EPG-X does not store all states
            y, x = numpy.meshgrid(
                *[range(0, n) for n in Z.shape], indexing="ij")
            mask = (x <= y)
            numpy.testing.assert_allclose(Z[mask], self.epgx["Zn0"][mask])
    
    def test_magnetization_transfer(self):
        for Class, apply_time_interval in zip(self.classes, self.apply_time_interval):
            species = sycomore.Species(self.T1_MT[0], self.T2_MT)
            M0 = [
                sycomore.Array[float](0, 0, z)
                for z in [1-self.f_MT, self.f_MT]]
            model = Class(species, self.T1_MT[1], *M0, self.k_MT)
            
            # RF duration
            tau = self.alpha/(self.gamma*self.B1)
            # Saturation rate
            W = numpy.pi * self.gamma**2 * self.B1**2 * self.G
            
            states = self._run(model, apply_time_interval, W*tau)
            
            # For MT, there is only a single pool of F states
            F = numpy.hstack(
                [states[:, -1:0:-1, 0, 1].conj(), states[:, :, 0, 0]])
            offset = F.shape[1]//2 - self.epgx["Fnmt"].shape[0]//2
            F = F[:, offset:-offset].T
            
            # WARNING: EPG-X does not store all states
            y, x = numpy.meshgrid(
                *[range(0, n) for n in F.shape], indexing="ij")
            distance = numpy.abs(x-offset)+numpy.abs(y-offset)
            mask = distance < self.npulse/2
            numpy.testing.assert_allclose(F[mask], self.epgx["Fnmt"][mask])
        
            # Test the two pools of Z states
            Z = states[:, :self.npulse//2+1, :, 2].transpose(1, 0, 2)
            y, x = numpy.meshgrid(
                *[range(0, n) for n in Z.shape[:2]], indexing="ij")
            mask = (x <= y)
            numpy.testing.assert_allclose(Z[mask], self.epgx["Znmt"][mask])
    
    def test_exchange(self):
        for Class, apply_time_interval in zip(self.classes, self.apply_time_interval):
            species = [
                sycomore.Species(T1, T2) for T1, T2 in zip(self.T1x, self.T2x)]
            M0 = [sycomore.Array[float](0, 0, z) for z in [1-self.fx, self.fx]]
            model = Class(*species, *M0, self.kx, self.delta_b)
        
            states = self._run(model, apply_time_interval)
            
            F = numpy.hstack(
                [states[..., -1:0:-1, :, 1].conj(), states[..., 0]])
            offset = F.shape[1]//2 - self.epgx["Fnx"].shape[0]//2
            F = F[:, offset:-offset].transpose(1, 0, 2)
            
            numpy.testing.assert_allclose(F, self.epgx["Fnx"], atol=1e-5)
            
            Z = states[..., 2]
            Z = Z[:, :-offset, :].transpose(1,0,2)
            numpy.testing.assert_allclose(Z, self.epgx["Znx"], atol=1e-5)
    
    def _run(self, model, apply_time_interval, saturation=None):
        states_array = (
            numpy.zeros((len(self.phi), len(self.phi), 2, 3), complex)
            if model.pools == 2
            else numpy.zeros((len(self.phi), len(self.phi), 3), complex))
        
        for i, phi in enumerate(self.phi):
            if saturation:
                model.apply_pulse(self.alpha, phi, saturation)
            else:
                model.apply_pulse(self.alpha, phi)
            
            states_array[i, :len(model)] = model.states
            
            apply_time_interval(model, self.TR)
        
        return states_array

if __name__ == "__main__":
    unittest.main()

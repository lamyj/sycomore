import unittest

import numpy
import sycomore
from sycomore.units import *

import os
import sys
sys.path.append(os.path.dirname(__file__))
from epg_test_case import EPGTestCase

class TestRegular(EPGTestCase):
    def test_empty(self):
        species = sycomore.Species(1000*ms, 100*ms)
        model = sycomore.epg.Regular(species)
        
        self._test_model(
            model, sycomore.TensorQ1([0]), [[0,0,1]])
    
    def test_pulse(self):
        species = sycomore.Species(1000*ms, 100*ms, 3*um**2/ms)
        model = sycomore.epg.Regular(species)
        model.apply_pulse(47*deg, 23*deg)
        
        self._test_model(
            model, [0],
            [[
                0.2857626571584661-0.6732146319308543j,
                0.2857626571584661+0.6732146319308543j,
                0.6819983600624985]])
    
    def test_gradient_default(self):
        species = sycomore.Species(1000*ms, 100*ms, 3*um**2/ms)
        model = sycomore.epg.Regular(species)
        model.apply_pulse(47*deg, 23*deg)
        model.shift()
        
        self._test_model(
            model, [0, 1],
            [
                [0, 0, 0.6819983600624985],
                [0.2857626571584661-0.6732146319308543j, 0, 0]])
    
    def test_gradient_multiple(self):
        species = sycomore.Species(1000*ms, 100*ms, 3*um**2/ms)
        model = sycomore.epg.Regular(species, unit_dephasing=1*mT/m*ms)
        model.apply_pulse(47*deg, 23*deg)
        
        model.shift(1*ms, 1*mT/m)
        self._test_model(
            model, sycomore.gamma*[0, 1]*mT/m*ms,
            [
                [0, 0, 0.6819983600624985],
                [0.2857626571584661-0.6732146319308543j, 0, 0]])
        
        model.shift(2*ms, 1*mT/m)
        self._test_model(
            model, sycomore.gamma*[0, 1, 2, 3]*mT/m*ms,
            [
                [0, 0, 0.6819983600624985],
                [0, 0, 0],
                [0, 0, 0],
                [0.2857626571584661-0.6732146319308543j, 0, 0]])
        
        model.shift(1*ms, -1*mT/m)
        self._test_model(
            model, sycomore.gamma*[0, 1, 2, 3, 4]*mT/m*ms,
            [
                [0, 0, 0.6819983600624985],
                [0, 0, 0],
                [0.2857626571584661-0.6732146319308543j, 0, 0],
                [0, 0, 0],
                [0, 0, 0]])
        
        with self.assertRaises(Exception):
            model.shift(1.5*ms, 1*mT/m)
    
    def test_relaxation(self):
        species = sycomore.Species(1000*ms, 100*ms, 3*um**2/ms)
        model = sycomore.epg.Regular(species)
        model.apply_pulse(47*deg, 23*deg)
        model.shift()
        model.relaxation(10*ms)
        
        self._test_model(
            model, [0, 1],
            [
                [0, 0, 0.6851625292479138],
                [0.2585687448743616-0.6091497893403431j, 0, 0]])
    
    def test_diffusion(self):
        species = sycomore.Species(1000*ms, 100*ms, 3*um**2/ms)
        model = sycomore.epg.Regular(species, unit_dephasing=20*mT/m*ms)
        model.apply_pulse(47*deg, 23*deg)
        model.shift()
        model.relaxation(10*ms)
        model.diffusion(10*ms, 2*mT/m)
        
        self._test_model(
            model, sycomore.gamma*[0, 20]*mT/m*ms,
            [
                [0, 0, 0.6851625292479138],
                [0.25805111586158685-0.60793033180597855j, 0, 0]])
    
    def test_off_resonance(self):
        species = sycomore.Species(1000*ms, 100*ms, 3*um**2/ms)
        model = sycomore.epg.Regular(species)
        model.delta_omega = 10*Hz
        model.apply_pulse(47*deg, 23*deg)
        model.shift()
        model.off_resonance(10*ms)
        
        self._test_model(
            model, [0, 1],
            [
                [0, 0, 0.6819983600624985], 
                [0.6268924782754024-0.37667500256027975j, 0, 0]])
    
    def test_time_interval(self):
        species = sycomore.Species(1000*ms, 100*ms, 3*um**2/ms)
        model = sycomore.epg.Regular(species, unit_dephasing=10*mT/m*ms)
        model.apply_pulse(47*deg, 23*deg)
        model.apply_time_interval(10*ms, 2*mT/m)
        
        self._test_model(
            model, sycomore.gamma*[0, 10, 20]*mT/m*ms,
            [
                [0, 0, 0.6851625292479138],
                [0, 0, 0],
                [0.2584947343504123-0.6089754314724013j, 0, 0]])
        
        model.apply_time_interval(10*ms, -2*mT/m)
        self._test_model(
            model, sycomore.gamma*[0]*mT/m*ms,
            [
                [
                    0.23382875968307784-0.5508660366970124j, 
                    0.23382875968307784+0.5508660366970124j,
                    0.6882952144238884]])
    
    def test_apply_time_interval_field_off_resonance(self):
        species = sycomore.Species(1000*ms, 100*ms, 3*um**2/ms)
        model = sycomore.epg.Regular(species, unit_dephasing=10*mT/m*ms)
        model.delta_omega = 10*Hz
        model.apply_pulse(47*deg, 23*deg)
        model.apply_time_interval(10*ms, 2*mT/m)
        
        self._test_model(
            model, sycomore.gamma*[0, 10, 20]*mT/m*ms,
            [
                [0, 0, 0.6851625292479138], 
                [0, 0, 0],
                [0.56707341067384409-0.34073208057155585j, 0, 0]])
    
    def test_apply_time_interval_species_off_resonance(self):
        species = sycomore.Species(1000*ms, 100*ms, 3*um**2/ms, delta_omega=10*Hz)
        model = sycomore.epg.Regular(species, unit_dephasing=10*mT/m*ms)
        model.apply_pulse(47*deg, 23*deg)
        model.apply_time_interval(10*ms, 2*mT/m)
        
        self._test_model(
            model, sycomore.gamma*[0, 10, 20]*mT/m*ms,
            [
                [0, 0, 0.6851625292479138], 
                [0, 0, 0],
                [0.56707341067384409-0.34073208057155585j, 0, 0]])
    
    def test_apply_time_interval_species_off_resonance(self):
        species = sycomore.Species(1000*ms, 100*ms, 3*um**2/ms, delta_omega=10*Hz)
        model = sycomore.epg.Regular(species, unit_dephasing=10*mT/m*ms)
        model.delta_omega = -10*Hz
        model.apply_pulse(47*deg, 23*deg)
        model.apply_time_interval(10*ms, 2*mT/m)
        
        self._test_model(
            model, sycomore.gamma*[0, 10, 20]*mT/m*ms,
            [
                [0, 0, 0.6851625292479138], 
                [0, 0, 0],
                [0.2584947343504123-0.6089754314724013j, 0, 0]])
    
    def test_elapsed(self):
        species = sycomore.Species(1000*ms, 100*ms)
        model = sycomore.epg.Regular(species, unit_dephasing=10*mT/m*ms)
        self.assertEqual(model.elapsed, 0*s)
        
        model.apply_time_interval(10*ms)
        self.assertEqual(model.elapsed, 10*ms)

if __name__ == "__main__":
    unittest.main()

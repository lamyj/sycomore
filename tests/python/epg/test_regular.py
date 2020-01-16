import unittest

import numpy
import sycomore
from sycomore.units import *

class TestRegular(unittest.TestCase):
    def test_empty(self):
        species = sycomore.Species(1000*ms, 100*ms)
        model = sycomore.epg.Regular(species)
        
        numpy.testing.assert_array_almost_equal(
            model.states, [[0,0,1]])
        self.assertEqual(model.states_count, 1)
        numpy.testing.assert_array_almost_equal(model.state(0), [0,0,1])
        numpy.testing.assert_almost_equal(model.echo, 0)
    
    def test_pulse(self):
        species = sycomore.Species(1000*ms, 100*ms, 3*um**2/ms)
        model = sycomore.epg.Regular(species)
        model.apply_pulse(47*deg, 23*deg)
        
        numpy.testing.assert_array_almost_equal(
            model.states, [[
                0.2857626571584661-0.6732146319308543j,
                0.2857626571584661+0.6732146319308543j,
                0.6819983600624985]])
        self.assertEqual(model.states_count, 1)
        numpy.testing.assert_array_almost_equal(model.state(0), [
            0.2857626571584661-0.6732146319308543j,
            0.2857626571584661+0.6732146319308543j,
            0.6819983600624985])
        numpy.testing.assert_almost_equal(
            model.echo, 0.2857626571584661-0.6732146319308543j)
    
    def test_gradient_default(self):
        species = sycomore.Species(1000*ms, 100*ms, 3*um**2/ms)
        model = sycomore.epg.Regular(species)
        model.apply_pulse(47*deg, 23*deg)
        model.shift()
        
        numpy.testing.assert_array_almost_equal(
            model.states, [
                [0, 0, 0.6819983600624985],
                [0.2857626571584661-0.6732146319308543j, 0, 0]])
        self.assertEqual(model.states_count, 2)
        numpy.testing.assert_array_almost_equal(
            model.state(0), [0, 0, 0.6819983600624985])
        numpy.testing.assert_array_almost_equal(
            model.state(1), [0.2857626571584661-0.6732146319308543j, 0, 0])
        numpy.testing.assert_almost_equal(model.echo, 0)
    
    def test_gradient_multiple(self):
        species = sycomore.Species(1000*ms, 100*ms, 3*um**2/ms)
        model = sycomore.epg.Regular(species, unit_gradient_area=1*mT/m*ms)
        model.apply_pulse(47*deg, 23*deg)
        
        model.shift(1*ms, 1*mT/m)
        numpy.testing.assert_array_almost_equal(
            model.states, [
                [0, 0, 0.6819983600624985],
                [0.2857626571584661-0.6732146319308543j, 0, 0]])
        self.assertEqual(model.states_count, 2)
        numpy.testing.assert_array_almost_equal(
            model.state(0), [0, 0, 0.6819983600624985])
        numpy.testing.assert_array_almost_equal(
            model.state(1), [0.2857626571584661-0.6732146319308543j, 0, 0])
        
        model.shift(2*ms, 1*mT/m)
        numpy.testing.assert_array_almost_equal(
            model.states, [
                [0, 0, 0.6819983600624985],
                [0, 0, 0],
                [0, 0, 0],
                [0.2857626571584661-0.6732146319308543j, 0, 0]])
        self.assertEqual(model.states_count, 4)
        numpy.testing.assert_array_almost_equal(
            model.state(0), [0, 0, 0.6819983600624985])
        numpy.testing.assert_array_almost_equal(model.state(1), [0,0,0])
        numpy.testing.assert_array_almost_equal(model.state(2), [0,0,0])
        numpy.testing.assert_array_almost_equal(
            model.state(3), [0.2857626571584661-0.6732146319308543j, 0, 0])
        
        model.shift(1*ms, -1*mT/m)
        numpy.testing.assert_array_almost_equal(
            model.states, [
                [0, 0, 0.6819983600624985],
                [0, 0, 0],
                [0.2857626571584661-0.6732146319308543j, 0, 0],
                [0, 0, 0],
                [0, 0, 0]])
        self.assertEqual(model.states_count, 5)
        numpy.testing.assert_array_almost_equal(
            model.state(0), [0, 0, 0.6819983600624985])
        numpy.testing.assert_array_almost_equal(model.state(1), [0,0,0])
        numpy.testing.assert_array_almost_equal(
            model.state(2), [0.2857626571584661-0.6732146319308543j, 0, 0])
        numpy.testing.assert_array_almost_equal(model.state(3), [0,0,0])
        numpy.testing.assert_array_almost_equal(model.state(4), [0,0,0])
        
        with self.assertRaises(Exception):
            model.shift(1.5*ms, 1*mT/m)
    
    def test_relaxation(self):
        species = sycomore.Species(1000*ms, 100*ms, 3*um**2/ms)
        model = sycomore.epg.Regular(species)
        model.apply_pulse(47*deg, 23*deg)
        model.shift()
        model.relaxation(10*ms)
        
        numpy.testing.assert_array_almost_equal(
            model.states, [
                [0, 0, 0.6851625292479138],
                [0.2585687448743616-0.6091497893403431j, 0, 0]])
        self.assertEqual(model.states_count, 2)
        numpy.testing.assert_array_almost_equal(
            model.state(0), [0, 0, 0.6851625292479138])
        numpy.testing.assert_array_almost_equal(
            model.state(1), [0.2585687448743616-0.6091497893403431j, 0, 0])
        numpy.testing.assert_almost_equal(model.echo, 0)
    
    def test_diffusion(self):
        species = sycomore.Species(1000*ms, 100*ms, 3*um**2/ms)
        model = sycomore.epg.Regular(species)
        model.apply_pulse(47*deg, 23*deg)
        model.shift()
        model.relaxation(10*ms)
        model.diffusion(10*ms, 2*mT/m)
        
        numpy.testing.assert_array_almost_equal(
            model.states, [
                [0, 0, 0.6851625292479138],
                [0.25805111586158685-0.6079303318059787j, 0, 0]])
        self.assertEqual(model.states_count, 2)
        numpy.testing.assert_array_almost_equal(
            model.state(0), [0, 0, 0.6851625292479138])
        numpy.testing.assert_array_almost_equal(
            model.state(1), [0.25805111586158685-0.6079303318059787j, 0, 0])
        numpy.testing.assert_almost_equal(model.echo, 0)
    
    def test_time_interval(self):
        species = sycomore.Species(1000*ms, 100*ms, 3*um**2/ms)
        model = sycomore.epg.Regular(species, unit_gradient_area=10*mT/m*ms)
        model.apply_pulse(47*deg, 23*deg)
        model.apply_time_interval(10*ms, 2*mT/m)
        
        numpy.testing.assert_array_almost_equal(
            model.states, [
                [0, 0, 0.6851625292479138],
                [0, 0, 0],
                [0.2584947343504123-0.6089754314724013j, 0, 0]])
        self.assertEqual(model.states_count, 3)
        numpy.testing.assert_array_almost_equal(
            model.state(0), [0, 0, 0.6851625292479138])
        numpy.testing.assert_array_almost_equal(model.state(1), [0,0,0])
        numpy.testing.assert_array_almost_equal(
            model.state(2), [0.2584947343504123-0.6089754314724013j, 0, 0])
        numpy.testing.assert_almost_equal(model.echo, 0)
        
        model.apply_time_interval(10*ms, -2*mT/m)
        numpy.testing.assert_array_almost_equal(
            model.states, [
                [
                    0.23262696138115807-0.5480347773241918j, 
                    0.23262696138115807+0.5480347773241918j,
                    0.6882952144238884],
                [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]])
        self.assertEqual(model.states_count, 5)
        numpy.testing.assert_array_almost_equal(
            model.state(0), [
                0.23262696138115807-0.5480347773241918j, 
                0.23262696138115807+0.5480347773241918j,
                0.6882952144238884])
        for i in range(1,5):
            numpy.testing.assert_array_almost_equal(
                model.state(i), [0, 0, 0])
        numpy.testing.assert_almost_equal(
            model.echo, 0.23262696138115807-0.5480347773241918j)
        
        
if __name__ == "__main__":
    unittest.main()

import unittest

import numpy
import sycomore
from sycomore.units import *

class TestDiscrete3D(unittest.TestCase):
    def setUp(self):
        self.species = sycomore.Species(1000*ms, 100*ms, 3*um**2/ms)
    
    def test_empty(self):
        model = sycomore.epg.Discrete3D(self.species)
        self._test_model(model, [3*[0*rad/m]], [[0, 0, 1]])
    
    def test_pulse(self):
        model = sycomore.epg.Discrete3D(self.species)
        model.apply_pulse(47*deg, 23*deg)
        
        self._test_model(
            model, [3*[0*rad/m]], 
            [[
                0.2857626571584661-0.6732146319308543j, 
                0.2857626571584661+0.6732146319308543j, 
                0.6819983600624985]])
    
    def test_positive_gradient_x(self):
        model = sycomore.epg.Discrete3D(self.species)
        model.apply_pulse(47*deg, 23*deg)
        model.shift(10*ms, [2*mT/m, 0*mT/m, 0*mT/m])
        
        self._test_model(
            model, 
            [3*[0*rad/m], [5350*rad/m, 0*rad/m, 0*rad/m]], 
            [
                [0, 0, 0.6819983600624985],
                [0.2857626571584661-0.6732146319308543j, 0, 0]])
    
    def test_positive_gradient_y(self):
        model = sycomore.epg.Discrete3D(self.species)
        model.apply_pulse(47*deg, 23*deg)
        model.shift(10*ms, [0*mT/m, 2*mT/m, 0*mT/m])
        
        self._test_model(
            model, 
            [3*[0*rad/m], [0*rad/m, 5350*rad/m, 0*rad/m]], 
            [
                [0, 0, 0.6819983600624985],
                [0.2857626571584661-0.6732146319308543j, 0, 0]])
    
    def test_positive_gradient_z(self):
        model = sycomore.epg.Discrete3D(self.species)
        model.apply_pulse(47*deg, 23*deg)
        model.shift(10*ms, [0*mT/m, 0*mT/m, 2*mT/m])
        
        self._test_model(
            model, 
            [3*[0*rad/m], [0*rad/m, 0*rad/m, 5350*rad/m]], 
            [
                [0, 0, 0.6819983600624985],
                [0.2857626571584661-0.6732146319308543j, 0, 0]])
    
    def test_negative_gradient_x(self):
        model = sycomore.epg.Discrete3D(self.species)
        model.apply_pulse(47*deg, 23*deg)
        model.shift(10*ms, [-2*mT/m, 0*mT/m, 0*mT/m])
        
        self._test_model(
            model, 
            [3*[0*rad/m], [5350*rad/m, 0*rad/m, 0*rad/m]], 
            [
                [0, 0, 0.6819983600624985], 
                [0, 0.2857626571584661+0.6732146319308543j, 0]])
    
    def test_negative_gradient_y(self):
        model = sycomore.epg.Discrete3D(self.species)
        model.apply_pulse(47*deg, 23*deg)
        model.shift(10*ms, [0*mT/m, -2*mT/m, 0*mT/m])
    
        self._test_model(
            model, 
            [3*[0*rad/m], [0*rad/m, 5350*rad/m, 0*rad/m]], 
            [
                [0, 0, 0.6819983600624985],
                [0, 0.2857626571584661+0.6732146319308543j, 0]])
    
    def test_negative_gradient_z(self):
        model = sycomore.epg.Discrete3D(self.species)
        model.apply_pulse(47*deg, 23*deg)
        model.shift(10*ms, [0*mT/m, 0*mT/m, -2*mT/m])
    
        self._test_model(
            model, 
            [3*[0*rad/m], [0*rad/m, 0*rad/m, 5350*rad/m]], 
            [
                [0, 0, 0.6819983600624985],
                [0, 0.2857626571584661+0.6732146319308543j, 0]])
    
    def test_multiple_gradient(self):
        model = sycomore.epg.Discrete3D(self.species)
        model.apply_pulse(47*deg, 23*deg)
        model.shift(10*ms, [-2*mT/m, 2*mT/m, -2*mT/m])
        model.apply_pulse(47*deg, 23*deg)
        model.shift(10*ms, [1*mT/m, -3*mT/m, 3*mT/m])
        
        self._test_model(
            model,
            [
                [0*rad/m, 0*rad/m, 0*rad/m],
                [2675*rad/m, -8026*rad/m, 8026*rad/m],
                [5350*rad/m, -5350*rad/m, 5350*rad/m],
                [8025*rad/m, -13376*rad/m, 13376*rad/m],
                [2675*rad/m, 2676*rad/m, -2676*rad/m]
            ],
            [
                [0, 0, 0.4651217631279373],
                [0.19488966354917586-0.45913127494692113j, 0, 0],
                [0, 0, -0.26743911843603135],
                [-0.045436496804645087+0.10704167849196657j, 0, 0],
                [0, 0.240326160353821+0.5661729534388877j, 0]
            ])
    
    def test_relaxation(self):
        model = sycomore.epg.Discrete3D(self.species)
        model.apply_pulse(47*deg, 23*deg)
        model.shift(10*ms, [2*mT/m, 0*mT/m, 0*mT/m])
        model.relaxation(10*ms)
        
        self._test_model(
            model, 
            [3*[0*rad/m], [5350*rad/m, 0*rad/m, 0*rad/m]], 
            [
                [0, 0, 0.6851625292479138],
                [0.2585687448743616-0.6091497893403431j, 0, 0]])
    
    def test_diffusion_x(self):
        model = sycomore.epg.Discrete3D(self.species)
        model.apply_pulse(47*deg, 23*deg)
        model.shift(10*ms, [2*mT/m, 0*mT/m, 0*mT/m])
        model.relaxation(10*ms)
        model.diffusion(10*ms, [2*mT/m, 0*mT/m, 0*mT/m])
        
        self._test_model(
            model, 
            [3*[0*rad/m], [5350*rad/m, 0*rad/m, 0*rad/m]], 
            [
                [0, 0, 0.6851625292479138],
                [0.25805117100742553-0.6079304617214332j, 0, 0]])
    
    def test_diffusion_z(self):
        model = sycomore.epg.Discrete3D(self.species)
        model.apply_pulse(47*deg, 23*deg)
        model.shift(10*ms, [0*mT/m, 0*mT/m, 2*mT/m])
        model.relaxation(10*ms)
        model.diffusion(10*ms, [0*mT/m, 0*mT/m, 2*mT/m])
        
        self._test_model(
            model, 
            [3*[0*rad/m], [0*rad/m, 0*rad/m, 5350*rad/m]], 
            [
                [0, 0, 0.6851625292479138],
                [0.25805117100742553-0.6079304617214332j, 0, 0]])
    
    def test_off_resonance(self):
        model = sycomore.epg.Discrete3D(self.species)
        model.delta_omega = 10*Hz
        model.apply_pulse(47*deg, 23*deg)
        model.shift(10*ms, [2*mT/m, 0*mT/m, 0*mT/m])
        model.off_resonance(10*ms)
        
        self._test_model(
            model, 
            [3*[0*rad/m], [5350*rad/m, 0*rad/m, 0*rad/m]],
            [
                [0, 0, 0.6819983600624985],
                [0.6268924782754024-0.37667500256027975j, 0, 0]])
    
    def test_time_interval_x(self):
        model = sycomore.epg.Discrete3D(self.species)
        model.apply_pulse(47*deg, 23*deg)
        model.apply_time_interval(10*ms, [2*mT/m, 0*mT/m, 0*mT/m])
        
        self._test_model(
            model, 
            [3*[0*rad/m], [5350*rad/m, 0*rad/m, 0*rad/m]], 
            [
                [0, 0, 0.6851625292479138],
                [0.2584947343504123-0.6089754314724013j, 0, 0]])
    
    def test_time_interval_z(self):
        model = sycomore.epg.Discrete3D(self.species)
        model.apply_pulse(47*deg, 23*deg)
        model.apply_time_interval(10*ms, [0*mT/m, 0*mT/m, 2*mT/m])
        
        self._test_model(
            model, 
            [3*[0*rad/m], [0*rad/m, 0*rad/m, 5350*rad/m]], 
            [
                [0, 0, 0.6851625292479138],
                [0.2584947343504123-0.6089754314724013j, 0, 0]])
    
    def test_apply_time_interval_field_off_resonance(self):
        model = sycomore.epg.Discrete3D(self.species)
        model.delta_omega = 10*Hz
        model.apply_pulse(47*deg, 23*deg)
        model.apply_time_interval(
            10*ms, [2*mT/m, 0*mT/m, 0*mT/m])
        
        self._test_model(
            model, 
            [3*[0*rad/m], [5350*rad/m, 0*rad/m, 0*rad/m]], 
            [
                [0, 0, 0.6851625292479138],
                [0.56707341067384409-0.34073208057155585j, 0, 0]])
    
    def test_apply_time_interval_species_off_resonance(self):
        model = sycomore.epg.Discrete3D(
            sycomore.Species(
                self.species.R1, self.species.R2, self.species.D, 
                delta_omega=10*Hz))
        model.apply_pulse(47*deg, 23*deg)
        model.apply_time_interval(10*ms, [2*mT/m, 0*mT/m, 0*mT/m])
        
        self._test_model(
            model, 
            [3*[0*rad/m], [5350*rad/m, 0*rad/m, 0*rad/m]], 
            [
                [0, 0, 0.6851625292479138],
                [0.56707341067384409-0.34073208057155585j, 0, 0]])
    
    def test_apply_time_interval_both_off_resonance(self):
        model = sycomore.epg.Discrete3D(
            sycomore.Species(
                self.species.R1, self.species.R2, self.species.D, 
                delta_omega=10*Hz))
        model.delta_omega = -10*Hz
        model.apply_pulse(47*deg, 23*deg)
        model.apply_time_interval(10*ms, [2*mT/m, 0*mT/m, 0*mT/m])
        
        self._test_model(
            model, 
            [3*[0*rad/m], [5350*rad/m, 0*rad/m, 0*rad/m]], 
            [
                [0, 0, 0.6851625292479138],
                [0.2584947343504123-0.6089754314724013j, 0, 0]])
    
    def test_refocalization(self):
        model = sycomore.epg.Discrete3D(self.species)
        model.apply_pulse(90*deg, 30*deg)
        model.apply_time_interval(10*ms, [2*mT/m, 0*mT/m, 0*mT/m])
        model.apply_pulse(120*deg, 0*deg)
        model.apply_time_interval(10*ms, [2*mT/m, 0*mT/m, 0*mT/m])
        
        self._test_model(
            model,
            [
                [0*rad/m, 0*rad/m, 0*rad/m],
                [5350*rad/m, 0*rad/m, 0*rad/m],
                [10700*rad/m, 0*rad/m, 0*rad/m]],
            [
                [
                    0.30684831950624042+0.53147687960193668j,
                    0.30684831950624042-0.53147687960193668j,
                    0.0050245860296255166],
                [
                    0-0.0077948398021822725j,
                    0,
                    -0.33555338970217136-0.19373183987203996j],
                [0.10210725404661349-0.17685495183007738j, 0, 0]
            ])
    
    def _test_model(self, model, orders, states):
        self.assertEqual(len(orders), len(model.orders))
        for o1, o2 in zip(orders, model.orders):
            self.assertSequenceEqual(o1, o2)
        
        self.assertEqual(model.states.shape, (len(orders), 3))
        numpy.testing.assert_array_almost_equal(states, model.states)
        for i, order in enumerate(orders):
            numpy.testing.assert_almost_equal(model.state(order), states[i])
        
        try:
            numpy.testing.assert_almost_equal(
                model.state(3*[0*rad/m])[0], model.echo)
        except Exception:
            self.assertEqual(model.echo, 0)
    
if __name__ == "__main__":
    unittest.main()

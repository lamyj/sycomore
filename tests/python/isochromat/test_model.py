import numpy
import unittest

import sycomore
from sycomore.units import *

class TestModel(unittest.TestCase):
    def test_pulse_uniform(self):
        positions = [[0*m, 0*m, 0*m]]
        model = sycomore.isochromat.Model(1*s, 0.1*s, [0, 0, 1], positions)
        op = model.build_pulse(90*deg, 60*deg)
        pulse = numpy.array([
            [
                [            0.25, numpy.sqrt(3)/4,  numpy.sqrt(3)/2, 0],
                [ numpy.sqrt(3)/4,            0.75,             -0.5, 0],
                [-numpy.sqrt(3)/2,             0.5,                0, 0],
                [               0,               0,                0, 1]]])
        numpy.testing.assert_almost_equal(op.array, pulse)
        
        default_phase = model.build_pulse(numpy.pi/3*rad)
        numpy.testing.assert_almost_equal(
            default_phase.array, model.build_pulse(numpy.pi/3*rad, 0*rad).array)
    
    def test_pulse_variable(self):
        positions = [[-1*m, 0*m, 0*m], [1*m, 0*m, 0*m]]
        model = sycomore.isochromat.Model(1*s, 0.1*s, [0, 0, 1], positions)
        
        op = model.build_pulse([90*deg, 60*deg], [60*deg, 90*deg])
        
        pulse = numpy.array([
            [
                [            0.25, numpy.sqrt(3)/4, numpy.sqrt(3)/2, 0],
                [ numpy.sqrt(3)/4,            0.75,            -0.5, 0],
                [-numpy.sqrt(3)/2,             0.5,               0, 0],
                [               0,               0,               0, 1]
            ], [
                [             0.5, 0, numpy.sqrt(3)/2, 0],
                [               0, 1,               0, 0],
                [-numpy.sqrt(3)/2, 0,             0.5, 0],
                [               0, 0,               0, 1]]])
        numpy.testing.assert_almost_equal(op.array, pulse)
        
        default_phase = model.build_pulse([90*deg, 60*deg])
        numpy.testing.assert_almost_equal(
            default_phase.array,
            model.build_pulse([90*deg, 60*deg], [0*deg, 0*deg]).array)
    
    def test_relaxation(self):
        positions = [[0*m, 0*m, 0*m], [1*m, 0*m, 0*m]]
        model = sycomore.isochromat.Model(
            [1*s, 2*s], [0.1*s, 0.2*s], [[0, 0, 2], [0, 0, 1]], positions)
        
        op = model.build_relaxation(1*ms)
        
        E1 = [numpy.exp(-1e-3/1), numpy.exp(-1e-3/2)]
        E2 = [numpy.exp(-1e-3/0.1), numpy.exp(-1e-3/0.2)]
        relaxation = [
            [
                [ E2[0],     0,     0,           0],
                [     0, E2[0],     0,           0],
                [     0,     0, E1[0], 2*(1-E1[0])],
                [     0,     0,     0,           1]
            ], [
                [ E2[1],     0,     0,           0],
                [     0, E2[1],     0,           0],
                [     0,     0, E1[1], 1*(1-E1[1])],
                [     0,     0,     0,           1]]]
        numpy.testing.assert_almost_equal(op.array, relaxation)
    
    def test_phase_accumulation(self):
        positions = [[0*m, 0*m, 0*m], [1*m, 0*m, 0*m]]
        model = sycomore.isochromat.Model(1*s, 0.1*s, [0, 0, 1], positions)
        
        op = model.build_phase_accumulation([30*deg, 60*deg])
        
        phase_accumulation = [
            [
                [ numpy.sqrt(3.)/2.,               -0.5, 0, 0],
                [               0.5, numpy.sqrt(3.)/2., 0, 0],
                [                 0,                 0, 1, 0],
                [                 0,                 0, 0, 1]
            ], [
                [               0.5, -numpy.sqrt(3.)/2., 0, 0],
                [ numpy.sqrt(3.)/2.,                0.5, 0, 0],
                [                 0,                  0, 1, 0],
                [                 0,                  0, 0, 1]]]
        numpy.testing.assert_almost_equal(op.array, phase_accumulation)
    
    def test_time_interval_uniform(self):
        positions = [[0*m, 0*m, 0*m], [1*mm, 2*mm, 3*mm]]
        model = sycomore.isochromat.Model(1*s, 0.1*s, [0, 0, 1], positions)
        
        op = model.build_time_interval(10*ms, 400*Hz, [20*mT/m, 0*mT/m, 10*mT/m])
        
        combined = model.build_phase_accumulation(2*numpy.pi*rad*400*10e-3)
        combined.pre_multiply(model.build_relaxation(10*ms))
        numpy.testing.assert_almost_equal(op.array[0], combined.array[0])
        
        combined = model.build_phase_accumulation(
            2*numpy.pi*rad*400*10e-3 + sycomore.gamma.magnitude*50e-6*10e-3)
        combined.pre_multiply(model.build_relaxation(10*ms))
        numpy.testing.assert_almost_equal(op.array[1], combined.array[1])
        
        default_frequency_and_gradient = model.build_time_interval(10*ms)
        numpy.testing.assert_almost_equal(
            default_frequency_and_gradient.array,
            model.build_time_interval(10*ms, 0*Hz, [0*T/m, 0*T/m, 0*T/m]).array)
        
        default_gradient = model.build_time_interval(10*ms, 400*Hz)
        numpy.testing.assert_almost_equal(
            default_gradient.array,
            model.build_time_interval(10*ms, 400*Hz, [0*T/m, 0*T/m, 0*T/m]).array);
    
    def test_time_interval_variable(self):
        positions = [[-1*mm, 2*mm, 3*mm], [1*mm, 0*mm, 3*mm]]
        model = sycomore.isochromat.Model(1*s, 0.1*s, [0, 0, 1], positions)
        
        op = model.build_time_interval(
            10*ms, [400*Hz, 600*Hz],
            [[20*mT/m, 0*mT/m, 10*mT/m], [15*mT/m, 17*mT/m, 0*mT/m]])
        
        combined = model.build_phase_accumulation(
            2*numpy.pi*rad*400*10e-3 + sycomore.gamma.magnitude*10e-6*10e-3)
        combined.pre_multiply(model.build_relaxation(10*ms))
        numpy.testing.assert_almost_equal(op.array[0], combined.array[0])
        
        combined = model.build_phase_accumulation(
            2*numpy.pi*rad*600*10e-3 + sycomore.gamma.magnitude*15e-6*10e-3)
        combined.pre_multiply(model.build_relaxation(10*ms))
        numpy.testing.assert_almost_equal(op.array[1], combined.array[1])
        
        default_gradient = model.build_time_interval(10*ms, [400*Hz, 600*Hz])
        numpy.testing.assert_almost_equal(
            default_gradient.array,
            model.build_time_interval(
                10*ms, [400*Hz, 600*Hz],
                [[0*T/m, 0*T/m, 0*T/m], [0*T/m, 0*T/m, 0*T/m]]).array);
    
    def test_T1(self):
        model_1 = sycomore.isochromat.Model(
            1*s, 2*s, [3, 4, 5], [[0*m, 0*m, 11*m], [0*m, 0*m, 12*m]])
        self._test_quantity_array(model_1.T1, [1*s, 1*s])
        
        model_2 = sycomore.isochromat.Model(
            [1*s, 2*s], [3*s, 4*s], [[5, 6, 7], [8, 9, 10]],
            [[0*m, 0*m, 11*m], [0*m, 0*m, 12*m]])
        self._test_quantity_array(model_2.T1, [1*s, 2*s])
    
    def test_T2(self):
        model_1 = sycomore.isochromat.Model(
            1*s, 2*s, [3, 4, 5], [[0*m, 0*m, 11*m], [0*m, 0*m, 12*m]])
        self._test_quantity_array(model_1.T2, [2*s, 2*s])
        
        model_2 = sycomore.isochromat.Model(
            [1*s, 2*s], [3*s, 4*s], [[5, 6, 7], [8, 9, 10]],
            [[0*m, 0*m, 11*m], [0*m, 0*m, 12*m]])
        self._test_quantity_array(model_2.T2, [3*s, 4*s])
    
    def test_M0(self):
        model_1 = sycomore.isochromat.Model(
            1*s, 2*s, [3, 4, 5], [[0*m, 0*m, 11*m], [0*m, 0*m, 12*m]])
        numpy.testing.assert_almost_equal(model_1.M0, [5, 5])
        
        model_2 = sycomore.isochromat.Model(
            [1*s, 2*s], [3*s, 4*s], [[5, 6, 7], [8, 9, 10]],
            [[0*m, 0*m, 11*m], [0*m, 0*m, 12*m]])
        numpy.testing.assert_almost_equal(model_2.M0, [7, 10])
    
    def test_magnetization(self):
        model_1 = sycomore.isochromat.Model(
            1*s, 2*s, [3, 4, 5], [[0*m, 0*m, 11*m], [0*m, 0*m, 12*m]])
        numpy.testing.assert_almost_equal(
            model_1.magnetization, [[3,4,5], [3,4,5]])
        
        model_2 = sycomore.isochromat.Model(
            [1*s, 2*s], [3*s, 4*s], [[5, 6, 7], [8, 9, 10]],
            [[0*m, 0*m, 11*m], [0*m, 0*m, 12*m]])
        numpy.testing.assert_almost_equal(
            model_2.magnetization, [[5,6,7], [8,9,10]])
    
    def test_posittions(self):
        model = sycomore.isochromat.Model(
            1*s, 2*s, [3, 4, 5], [[0*m, 0*m, 11*m], [0*m, 0*m, 12*m]])
        self._test_quantity_array(
            model.positions, [[0*m,0*m,11*m], [0*m,0*m,12*m]])
    
    def test_apply(self):
        operator = sycomore.isochromat.Operator(
            [
                [
                    [1, 2, 3, 10],
                    [4, 5, 6, 20],
                    [7, 8, 9, 30],
                    [0, 0, 0, 1],
                ], [
                    [10, 11, 12, 40],
                    [13, 14, 15, 50],
                    [16, 17, 18, 60],
                    [0, 0, 0, 1],
                ]])
        
        model = sycomore.isochromat.Model(
            [1*s, 1*s], [1*s, 1*s],
            [[19, 20, 21], [22, 23, 24]], [[0*m, 0*m, 0*m], [0*m, 0*m, 1*m]])
        model.apply(operator)
    
        magnetization = [
            [132, 322, 512],
            [801, 1018, 1235]]
        
        numpy.testing.assert_almost_equal(model.magnetization, magnetization)
    
    def _test_quantity_array(self, left, right):
        self.assertEqual(numpy.shape(left), numpy.shape(right))
        self.assertSequenceEqual(
            [x.dimensions for x in numpy.ravel(left)],
            [x.dimensions for x in numpy.ravel(right)])
        self.assertSequenceEqual(
            [x.magnitude for x in numpy.ravel(left)],
            [x.magnitude for x in numpy.ravel(right)])
    
if __name__ == "__main__":
    unittest.main()

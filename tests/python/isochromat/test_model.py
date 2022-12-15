import numpy
import unittest

import sycomore

class TestModel(unittest.TestCase):
    def test_pulse_uniform(self):
        positions = [[0, 0, 0]]
        model = sycomore.isochromat.Model(1, 0.1, [0, 0, 1], positions)
        op = model.build_pulse(numpy.pi/2, numpy.pi/3)
        pulse = numpy.array([
            [
                [            0.25, numpy.sqrt(3)/4,  numpy.sqrt(3)/2, 0],
                [ numpy.sqrt(3)/4,            0.75,             -0.5, 0],
                [-numpy.sqrt(3)/2,             0.5,                0, 0],
                [               0,               0,                0, 1]]])
        numpy.testing.assert_almost_equal(op.array, pulse)
    
    def test_pulse_variable(self):
        positions = [[-1, 0, 0], [1, 0, 0]]
        model = sycomore.isochromat.Model(1, 0.1, [0, 0, 1], positions)
        
        op = model.build_pulse(
            [numpy.pi/2, numpy.pi/3], [numpy.pi/3, numpy.pi/2])
        
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
    
    def test_relaxation(self):
        positions = [[0, 0, 0], [1, 0, 0]]
        model = sycomore.isochromat.Model(
            [1, 2], [0.1, 0.2], [[0, 0, 2], [0, 0, 1]], positions)
        
        op = model.build_relaxation(1e-3)
        
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
        positions = [[0, 0, 0], [1, 0, 0]]
        model = sycomore.isochromat.Model(1, 0.1, [0, 0, 1], positions)
        
        op = model.build_phase_accumulation([numpy.pi/6, numpy.pi/3])
        
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
        positions = [[0, 0, 0], [1e-3, 2e-3, 3e-3]]
        model = sycomore.isochromat.Model(1, 0.1, [0, 0, 1], positions)
        
        op = model.build_time_interval(10e-3, 400, [20e-3, 0, 10e-3])
        
        combined = model.build_phase_accumulation(2*numpy.pi*400*10e-3)
        combined.preMultiply(model.build_relaxation(10e-3))
        numpy.testing.assert_almost_equal(op.array[0], combined.array[0])
        
        combined = model.build_phase_accumulation(
            2*numpy.pi*400*10e-3 + sycomore.gamma.magnitude*50e-6*10e-3)
        combined.preMultiply(model.build_relaxation(10e-3))
        numpy.testing.assert_almost_equal(op.array[1], combined.array[1])
    
    def test_time_interval_variable(self):
        positions = [[-1e-3, 2e-3, 3e-3], [1e-3, 0, 3e-3]]
        model = sycomore.isochromat.Model(1, 0.1, [0, 0, 1], positions)
        
        op = model.build_time_interval(
            10e-3, [400, 600], [[20e-3, 0, 10e-3], [15e-3, 17e-3, 0]])
        
        combined = model.build_phase_accumulation(
            2*numpy.pi*400*10e-3 + sycomore.gamma.magnitude*10e-6*10e-3)
        combined.preMultiply(model.build_relaxation(10e-3))
        numpy.testing.assert_almost_equal(op.array[0], combined.array[0])
        
        combined = model.build_phase_accumulation(
            2*numpy.pi*600*10e-3 + sycomore.gamma.magnitude*15e-6*10e-3)
        combined.preMultiply(model.build_relaxation(10e-3))
        numpy.testing.assert_almost_equal(op.array[1], combined.array[1])
    
    def test_T1(self):
        model_1 = sycomore.isochromat.Model(1, 2, [3,4,5], [[0,0,11], [0,0,12]])
        numpy.testing.assert_almost_equal(model_1.T1, [1])
        
        model_2 = sycomore.isochromat.Model(
            [1,2], [3,4], [[5,6,7,1], [8,9,10,1]], [[0,0,11], [0,0,12]])
        numpy.testing.assert_almost_equal(model_2.T1, [1, 2])
    
    def test_T2(self):
        model_1 = sycomore.isochromat.Model(1, 2, [3,4,5], [[0,0,11], [0,0,12]])
        numpy.testing.assert_almost_equal(model_1.T2, [2])
        
        model_2 = sycomore.isochromat.Model(
            [1,2], [3,4], [[5,6,7,1], [8,9,10,1]], [[0,0,11], [0,0,12]])
        numpy.testing.assert_almost_equal(model_2.T2, [3, 4])
    
    def test_M0(self):
        model_1 = sycomore.isochromat.Model(1, 2, [3,4,5], [[0,0,11], [0,0,12]])
        numpy.testing.assert_almost_equal(model_1.M0, [5])
        
        model_2 = sycomore.isochromat.Model(
            [1,2], [3,4], [[5,6,7,1], [8,9,10,1]], [[0,0,11], [0,0,12]])
        numpy.testing.assert_almost_equal(model_2.M0, [7,10])
    
    def test_magnetization(self):
        model_1 = sycomore.isochromat.Model(1, 2, [3,4,5], [[0,0,11], [0,0,12]])
        numpy.testing.assert_almost_equal(
            model_1.magnetization, [[3,4,5,1], [3,4,5,1]])
        
        model_2 = sycomore.isochromat.Model(
            [1,2], [3,4], [[5,6,7,1], [8,9,10,1]], [[0,0,11], [0,0,12]])
        numpy.testing.assert_almost_equal(
            model_2.magnetization, [[5,6,7,1], [8,9,10,1]])
    
    def test_posittions(self):
        model = sycomore.isochromat.Model(1, 2, [3,4,5], [[0,0,11], [0,0,12]])
        numpy.testing.assert_almost_equal(model.positions, [[0,0,11], [0,0,12]])
    
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
            [1, 1], [1, 1],
            [[19, 20, 21, 1], [22, 23, 24, 1]], [[0, 0, 0], [0, 0, 1]])
        model.apply(operator)
    
        magnetization = [
            [132, 322, 512, 1],
            [801, 1018, 1235, 1]]
        
        numpy.testing.assert_almost_equal(model.magnetization, magnetization)
    
if __name__ == "__main__":
    unittest.main()

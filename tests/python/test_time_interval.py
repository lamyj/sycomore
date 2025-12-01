import numpy
import unittest

import sycomore
from sycomore.units import *

from test_case import TestCase

class TestTimeInterval(TestCase):
    def test_duration_constructor(self):
        interval = sycomore.TimeInterval(1*ms)
        self.assertEqual(interval.duration, 1*ms)
        self.assertEqual(
            interval.gradient_amplitude, sycomore.Vector3Q(3*[0.*T/m]))
        self.assertEqual(
            interval.gradient_area, sycomore.Vector3Q(3*[0.*T/m*s]))
        self.assertEqual(
            interval.gradient_dephasing, sycomore.Vector3Q(3*[0.*rad/m]))
    
    def test_dephasing_scalar_constructor(self):
        interval = sycomore.TimeInterval(1.*ms, 2.*rad/dm)
        self.assertEqual(interval.duration, 1e-3*s)
        self.assertEqual(
            interval.gradient_amplitude,
            sycomore.Vector3Q(3*[74.76015355016015*uT/m]))
        self.assertEqual(
            interval.gradient_area,
            sycomore.Vector3Q(3*[74.76015355016015*uT/m*ms]))
        self.assertEqual(
            interval.gradient_dephasing, sycomore.Vector3Q(3*[20.*rad/m]))

    def test_dephasing_vector_constructor(self):
        interval = sycomore.TimeInterval(
            1.*ms, [2*rad/dm, 4*rad/m, 8*rad/dam])
        self.assertEqual(interval.duration, 1e-3*s)
        self.assertEqual(
            interval.gradient_amplitude, sycomore.Vector3Q([
                74.76015355016015*uT/m, 
                14.95203071003203*uT/m, 
                2.9904061420064063*uT/m]))
        self.assertEqual(
            interval.gradient_area, sycomore.Vector3Q([
                74.76015355016015*uT/m*ms, 
                14.95203071003203*uT/m*ms, 
                2.9904061420064063*uT/m*ms]))
        self.assertEqual(
            interval.gradient_dephasing, sycomore.Vector3Q([
                20*rad/m, 4*rad/m, 0.8*rad/m]))
    
    def test_amplitude_scalar_constructor(self):
        interval = sycomore.TimeInterval(1.*ms, 2.*mT/dm)
        self.assertEqual(interval.duration, 1e-3*s)
        self.assertEqual(
            interval.gradient_amplitude, sycomore.Vector3Q(3*[20*mT/m]))
        self.assertEqual(
            interval.gradient_area, sycomore.Vector3Q(3*[20*mT/m*ms]))
        self.assertEqual(
            interval.gradient_dephasing,
            sycomore.Vector3Q(3*[5350.4437993378515*rad/m]))

    def test_amplitude_vector_constructor(self):
        interval = sycomore.TimeInterval(1.*ms, [2.*mT/dm, 4.*mT/m, 8.*mT/dam])
        self.assertEqual(interval.duration, 1e-3*s)
        self.assertEqual(
            interval.gradient_amplitude, 
            sycomore.Vector3Q([20.*mT/m, 4.*mT/m, 0.8*mT/m]))
        self.assertEqual(
            interval.gradient_area,
            sycomore.Vector3Q([20.*mT/m*ms, 4.*mT/m*ms, 0.8*mT/m*ms]))
        self.assertEqual(
            interval.gradient_dephasing, sycomore.Vector3Q([
                5350.4437993378515*rad/m,
                1070.0887598675702*rad/m,
                214.01775197351404*rad/m]))
    
    def test_gradient_properties(self):
        amplitude = sycomore.Vector3Q([20*mT/m, 40*mT/m, 80*mT/m])
        area = sycomore.Vector3Q([20*mT/m*ms, 40*mT/m*ms, 80*mT/m*ms])
        dephasing = sycomore.Vector3Q([
            5350.4437993378515*rad/m,
            10700.887598675701*rad/m,
            21401.775197351402*rad/m])
        
        for property in ["amplitude", "area", "dephasing"]:
            interval = sycomore.TimeInterval(1.*ms)
            
            setattr(
                interval, "gradient_{}".format(property), locals()[property][0])
            self.assertEqual(
                interval.gradient_amplitude,
                sycomore.Vector3Q(3*[amplitude[0]]))
            self.assertEqual(
                interval.gradient_area, sycomore.Vector3Q(3*[area[0]]))
            self.assertEqual(
                interval.gradient_dephasing,
                sycomore.Vector3Q(3*[dephasing[0]]))
            
            setattr(
                interval, "gradient_{}".format(property), locals()[property])
            self.assertEqual(interval.gradient_amplitude, amplitude)
            self.assertEqual(interval.gradient_area, area)
            self.assertEqual(interval.gradient_dephasing, dephasing)
            
            interval.set_gradient(locals()[property][1])
            self.assertEqual(
                interval.gradient_amplitude,
                sycomore.Vector3Q(3*[amplitude[1]]))
            self.assertEqual(
                interval.gradient_area, sycomore.Vector3Q(3*[area[1]]))
            self.assertEqual(
                interval.gradient_dephasing,
                sycomore.Vector3Q(3*[dephasing[1]]))
            
            interval.set_gradient(locals()[property])
            self.assertEqual(interval.gradient_amplitude, amplitude)
            self.assertEqual(interval.gradient_area, area)
            self.assertEqual(interval.gradient_dephasing, dephasing)
    
    def test_comparison(self):
        interval_1 = sycomore.TimeInterval(1.*ms, 2*T/m)
        interval_2 = sycomore.TimeInterval(1.*ms, 2*T/m)
        interval_3 = sycomore.TimeInterval(1.*ms, [2*T/m, 2*T/m, 2*T/m])
        
        self.assertTrue(interval_1 == interval_2)
        self.assertTrue(interval_1 == interval_3)
        self.assertFalse(interval_1 != interval_2)
        self.assertFalse(interval_1 != interval_3)
        
        interval_4 = sycomore.TimeInterval(4.*ms, 2*T/m)
        self.assertFalse(interval_1 == interval_4)
        self.assertTrue(interval_1 != interval_4)
        
        interval_5 = sycomore.TimeInterval(1.*ms, 4*T/m)
        self.assertFalse(interval_1 == interval_5)
        self.assertTrue(interval_1 != interval_5)
    
    def test_shortest_scalar(self):
        M = 1000*rad/m
        G_max = 25*mT/m
        time_interval = sycomore.TimeInterval.shortest(M, G_max)
        duration = M/(G_max*sycomore.gamma_bar)
        
        self.assertEqual(
            time_interval.gradient_dephasing, sycomore.Vector3Q(3*[M]))
        self.assertEqual(
            time_interval.gradient_amplitude, sycomore.Vector3Q(3*[G_max]))
    
    def test_shortest_vector(self):
        M = [1000*rad/m, 2000*rad/m, 3000*rad/m]
        G_max = 25*mT/m
        time_interval = sycomore.TimeInterval.shortest(M, G_max)
        duration = M[2]/(G_max*sycomore.gamma_bar)
        
        self.assertEqual(time_interval.gradient_dephasing, sycomore.Vector3Q(M))
        self.assertTrue(time_interval.gradient_amplitude[0] < G_max)
        self.assertTrue(time_interval.gradient_amplitude[1] < G_max)
        self.assertEqual(time_interval.gradient_amplitude[2], G_max)

if __name__ == "__main__":
    unittest.main()

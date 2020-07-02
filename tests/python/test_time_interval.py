import numpy
import unittest

import sycomore
from sycomore.units import *

class TestTimeInterval(unittest.TestCase):
    def test_duration_constructor(self):
        interval = sycomore.TimeInterval(1*ms)
        self.assertEqual(interval.duration, 1*ms)
        self._test_quantity_array(interval.gradient_amplitude, 3*[0.*T/m])
        self._test_quantity_array(interval.gradient_area, 3*[0.*T/m*s])
        self._test_quantity_array(interval.gradient_dephasing, 3*[0.*rad/m])
        self._test_quantity_array(interval.gradient_moment, 3*[0.*rad/m])
            
    def test_dephasing_scalar_constructor(self):
        interval = sycomore.TimeInterval(1.*ms, 2.*rad/dm)
        self.assertEqual(interval.duration, 1e-3*s)
        self._test_quantity_array(
            interval.gradient_amplitude, 3*[74.76015355016015*uT/m])
        self._test_quantity_array(
            interval.gradient_area, 3*[74.76015355016015*uT/m*ms])
        self._test_quantity_array(interval.gradient_dephasing, 3*[20.*rad/m])
        self._test_quantity_array(interval.gradient_moment, 3*[20.*rad/m])

    def test_dephasing_vector_constructor(self):
        interval = sycomore.TimeInterval(1.*ms, [2*rad/dm, 4*rad/m, 8*rad/dam])
        self.assertEqual(interval.duration, 1e-3*s)
        self._test_quantity_array(
            interval.gradient_amplitude, [
                74.76015355016015*uT/m, 
                14.95203071003203*uT/m, 
                2.9904061420064063*uT/m])
        self._test_quantity_array(
            interval.gradient_area, [
                74.76015355016015*uT/m*ms, 
                14.95203071003203*uT/m*ms, 
                2.9904061420064063*uT/m*ms])
        self._test_quantity_array(
            interval.gradient_dephasing, [20*rad/m, 4*rad/m, 0.8*rad/m])
        self._test_quantity_array(
            interval.gradient_moment, [20*rad/m, 4*rad/m, 0.8*rad/m])
    
    def test_amplitude_scalar_constructor(self):
        interval = sycomore.TimeInterval(1.*ms, 2.*mT/dm)
        self.assertEqual(interval.duration, 1e-3*s)
        self._test_quantity_array(interval.gradient_amplitude, 3*[20*mT/m])
        self._test_quantity_array(interval.gradient_area, 3*[20*mT/m*ms])
        self._test_quantity_array(
            interval.gradient_dephasing, 3*[5350.4437993378515*rad/m])
        self._test_quantity_array(
            interval.gradient_moment, 3*[5350.4437993378515*rad/m])

    def test_amplitude_vector_constructor(self):
        interval = sycomore.TimeInterval(1.*ms, [2.*mT/dm, 4.*mT/m, 8.*mT/dam])
        self.assertEqual(interval.duration, 1e-3*s)
        self._test_quantity_array(
            interval.gradient_amplitude, [20.*mT/m, 4.*mT/m, 0.8*mT/m])
        self._test_quantity_array(
            interval.gradient_area, [20.*mT/m*ms, 4.*mT/m*ms, 0.8*mT/m*ms])
        self._test_quantity_array(
            interval.gradient_dephasing, [
                5350.4437993378515*rad/m,
                1070.0887598675702*rad/m,
                214.01775197351404*rad/m])
        self._test_quantity_array(
            interval.gradient_moment, [
                5350.4437993378515*rad/m,
                1070.0887598675702*rad/m,
                214.01775197351404*rad/m])
    
    def test_gradient_properties(self):
        amplitude = [20*mT/m, 40*mT/m, 80*mT/m]
        area = [20*mT/m*ms, 40*mT/m*ms, 80*mT/m*ms]
        dephasing = [
            5350.4437993378515*rad/m,
            10700.887598675701*rad/m,
            21401.775197351402*rad/m]
        moment = dephasing
        
        for property in ["amplitude", "area", "dephasing", "moment"]:
            interval = sycomore.TimeInterval(1.*ms)
            
            setattr(
                interval, "gradient_{}".format(property), locals()[property][0])
            self._test_quantity_array(
                interval.gradient_amplitude, 3*[amplitude[0]])
            self._test_quantity_array(interval.gradient_area, 3*[area[0]])
            self._test_quantity_array(
                interval.gradient_dephasing, 3*[dephasing[0]])
            self._test_quantity_array(interval.gradient_moment, 3*[moment[0]])
            
            setattr(
                interval, "gradient_{}".format(property), locals()[property])
            self._test_quantity_array(interval.gradient_amplitude, amplitude)
            self._test_quantity_array(interval.gradient_area, area)
            self._test_quantity_array(interval.gradient_dephasing, dephasing)
            self._test_quantity_array(interval.gradient_moment, moment)
            
            interval.set_gradient(locals()[property][1])
            self._test_quantity_array(
                interval.gradient_amplitude, 3*[amplitude[1]])
            self._test_quantity_array(interval.gradient_area, 3*[area[1]])
            self._test_quantity_array(
                interval.gradient_dephasing, 3*[dephasing[1]])
            self._test_quantity_array(interval.gradient_moment, 3*[moment[1]])
            
            interval.set_gradient(locals()[property])
            self._test_quantity_array(interval.gradient_amplitude, amplitude)
            self._test_quantity_array(interval.gradient_area, area)
            self._test_quantity_array(interval.gradient_dephasing, dephasing)
            self._test_quantity_array(interval.gradient_moment, moment)
    
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
    
    def _test_quantity_array(self, left, right):
        self.assertSequenceEqual(
            [x.dimensions for x in left], [x.dimensions for x in right])
        self.assertTrue(numpy.allclose(
            [x.magnitude for x in left], [x.magnitude for x in right]))
        
if __name__ == "__main__":
    unittest.main()

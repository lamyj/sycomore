import numpy
import unittest

import sycomore
from sycomore.units import *

class TestBloch(unittest.TestCase):
    def test_pulse(self):
        M = sycomore.bloch.pulse(47*deg, 23*deg)
        self.assertTrue(numpy.allclose(
            M,
            [[ 0.95145043, 0.11437562,  0.28576266, 0.        ],
             [ 0.11437562, 0.73054793, -0.67321463, 0.        ],
             [-0.28576266, 0.67321463,  0.68199836, 0.        ],
             [ 0.        , 0.        ,  0.        , 1.        ]]))

    def test_time_interval(self):
        M = sycomore.bloch.time_interval(
            sycomore.Species(1000*ms, 100*ms, delta_omega=420*Hz), 10*ms)
        self.assertTrue(numpy.allclose(
            M,
            [[ 0.27961014, -0.86055152,  0.        ,  0.        ],
             [ 0.86055152,  0.27961014,  0.        ,  0.        ],
             [ 0.        ,  0.        ,  0.99004983,  0.00995017],
             [ 0.        ,  0.        ,  0.        ,  1.        ]]
        ))

    def test_relaxation(self):
        M = sycomore.bloch.relaxation(sycomore.Species(1000*ms, 100*ms), 10*ms)
        self.assertTrue(numpy.allclose(
            M,
            [[0.90483742, 0.        , 0.        , 0.        ],
             [0.        , 0.90483742, 0.        , 0.        ],
             [0.        , 0.        , 0.99004983, 0.00995017],
             [0.        , 0.        , 0.        , 1.        ]]))

    def test_phase_accumulation(self):
        M = sycomore.bloch.phase_accumulation(numpy.pi/6*rad)
        self.assertTrue(numpy.allclose(
            M,
            [[ 0.8660254, -0.5      , 0.       , 0.       ],
             [ 0.5      ,  0.8660254, 0.       , 0.       ],
             [ 0.       ,  0.       , 1.       , 0.       ],
             [ 0.       ,  0.       , 0.       , 1.       ]]))

if __name__ == "__main__":
    unittest.main()

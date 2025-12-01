import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from test_case import TestCase

class EPGTestCase(TestCase):
    def _test_model(self, model, orders, states):
        self.assertAllClose(orders, model.orders)
        self.assertAllClose(states, model.states)
        
        self.assertEqual(model.states.shape, (len(orders), 3))
        for i, order in enumerate(orders):
            self.assertAllClose(model.state(i), states[i])
            self.assertAllClose(model.state(order), states[i])
        
        self.assertAllClose(states[0][0], model.echo)

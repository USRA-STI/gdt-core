import unittest

from gdt.misc import derivative

class TestDerivative(unittest.TestCase):

    def test_derivative(self):
        def f(x):
            return x ** 3 + x ** 2

        value = derivative(f, 1.0, dx=1e-6)
        self.assertAlmostEqual(value, 4.9999999999217337)


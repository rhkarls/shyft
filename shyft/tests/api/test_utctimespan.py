from shyft.api import TimeSpan, deltaminutes, deltahours
import unittest


class VerifyTimeSpan(unittest.TestCase):

    def test_construct(self):
        dt = TimeSpan(3)
        self.assertIsNotNone(dt)
        self.assertEqual(dt.seconds, 3)
        df = TimeSpan(3.123456)
        self.assertEqual(df.seconds, 3.123456)
        dn = TimeSpan()
        self.assertEqual(dn.seconds, 0)
        self.assertAlmostEqual(TimeSpan(1.1234567).seconds, 1.123457)  # note rounding
        self.assertEqual(repr(TimeSpan(1.234)), r'TimeSpan(1.234000)')
        self.assertEqual(repr(TimeSpan(1234)), r'TimeSpan(1234)')
        self.assertEqual(str(TimeSpan(-1.234)), r'-1.234000s')
        self.assertEqual(str(TimeSpan(1234)), r'1234s')

    def test_bin_ops(self):
        a = TimeSpan(1)
        b = TimeSpan(2)
        self.assertFalse(a == b)
        self.assertTrue(a != b)
        self.assertTrue(a == a)
        self.assertEqual(a + b, TimeSpan(3))
        self.assertEqual(a + 2, TimeSpan(3))
        self.assertEqual(1 + b, TimeSpan(3))
        self.assertEqual(a - b, TimeSpan(-1))
        self.assertEqual(a - 2, TimeSpan(-1))
        self.assertEqual(1 - b, TimeSpan(-1))
        self.assertEqual(a * 10, TimeSpan(10))
        self.assertEqual(a * 10.5, TimeSpan(10.5))
        self.assertEqual(10 * a, TimeSpan(10))
        self.assertEqual(-a * 10, TimeSpan(-10))
        self.assertEqual(TimeSpan(10) / TimeSpan(5), 2)  # note special chrono semantics here
        self.assertEqual(TimeSpan(5) / TimeSpan(10), TimeSpan(0.5))  # no surprise
        self.assertEqual(TimeSpan(5) // TimeSpan(10), TimeSpan(0))  # no surprise
        self.assertEqual(5 / TimeSpan(10), TimeSpan(0.5))  # no surprise
        self.assertEqual(5 // TimeSpan(10), TimeSpan(0))  # no surprise
        self.assertEqual(TimeSpan(5) / 10, TimeSpan(0.5))  # no surprise
        self.assertEqual(TimeSpan(5) // 10, TimeSpan(0))  # no surprise
        self.assertEqual(5 // TimeSpan(10), TimeSpan(0))  # no surprise
        self.assertEqual(TimeSpan(10) % TimeSpan(3), TimeSpan(1))  # ok, useful
        self.assertTrue(a < 2)
        self.assertTrue(a < 2.5)
        self.assertTrue(a <= 2)
        self.assertTrue(a <= 2.5)
        self.assertTrue(a > -2)
        self.assertTrue(a > -2.5)
        self.assertTrue(a >= -2)
        self.assertTrue(a >= -2.5)

    def test_delta(self):
        x = deltahours(1)
        self.assertEqual(x, TimeSpan(3600))
        m = deltaminutes(2)
        self.assertEqual(m, TimeSpan(120))

    def test_to_float(self):
        x = TimeSpan(1.1234)
        d = float(x)
        self.assertAlmostEqual(1.1234, d)

    def test_to_int(self):
        x = TimeSpan(1.1234)
        d = int(x)
        self.assertEqual(1, d)

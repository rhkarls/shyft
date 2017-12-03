from shyft.api import Calendar, UtcTime, TimeSpan, deltaminutes, deltahours
import unittest


class VerifyTUtcTime(unittest.TestCase):

    def test_construct(self):
        t = UtcTime(3)
        self.assertIsNotNone(t)
        self.assertEqual(t.seconds, 3)
        df = UtcTime(3.123456)
        self.assertEqual(df.seconds, 3.123456)
        dn = UtcTime()
        self.assertEqual(dn.seconds, 0)
        self.assertAlmostEqual(UtcTime(1.1234567).seconds, 1.123457)  # note rounding
        self.assertEqual(repr(UtcTime(1.234)), r'UtcTime(1.234000)')
        self.assertEqual(repr(UtcTime(1234)), r'UtcTime(1234)')
        self.assertEqual(str(UtcTime(-3600)), r'1969-12-31T23:00:00Z')
        self.assertEqual(str(UtcTime(0)), r'1970-01-01T00:00:00Z')
        self.assertEqual(str(UtcTime(r'1970-01-01T00:00:00Z')), r'1970-01-01T00:00:00Z')
        self.assertEqual(str(UtcTime(r'1970-01-01T01:02:03.123456-01:30')), r'1970-01-01T02:32:03.123456Z')

    def test_bin_ops(self):
        a = UtcTime(1)
        b = UtcTime(2)
        self.assertFalse(a == b)
        self.assertTrue(a != b)
        self.assertTrue(a == a)
        self.assertEqual(a + TimeSpan(3), UtcTime(4))
        self.assertEqual(a + 3, UtcTime(4))
        self.assertEqual(a - TimeSpan(3), UtcTime(-2))
        self.assertEqual(a - 3, UtcTime(-2))
        self.assertEqual(a - b, TimeSpan(-1))
        self.assertEqual(3*a-b, TimeSpan(1))

        self.assertTrue(a < b)
        self.assertTrue(a <= b)
        self.assertTrue(b > a)
        self.assertTrue(b >= a)
        self.assertTrue(b < 10.0)
        self.assertTrue(b < 10)
        self.assertTrue(b <= 10.0)
        self.assertTrue(b <= 10)
        self.assertTrue(b < r'1990-01-01T00:00:00Z')
        self.assertTrue(b <= r'1990-01-01T00:00:00Z')


    def test_floor(self):
        a = UtcTime(3601)
        dt = TimeSpan(3600)
        f = a.floor(dt)
        self.assertEqual(f, UtcTime(3600))

    def test_to_float(self):
        x = UtcTime(1.1234)
        d = float(x)
        self.assertAlmostEqual(1.1234, d)

    def test_to_int(self):
        x = UtcTime(1.1234)
        d = int(x)
        self.assertEqual(1, d)
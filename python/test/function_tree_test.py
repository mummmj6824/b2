from libpybertini import *
import unittest
import cmath

dbltol = 1e-15;
class NumberTest(unittest.TestCase):
    def test_eval_float(self):
        z = Float(3.4, 5.7)
        zd = z.EvalDbl()
        self.assertEqual(zd.real,3.4)
        self.assertEqual(zd.imag,5.7)

    def test_diff_float(self):
        z = Float(3.4, 5.7)
        dz = z.Differentiate()
        dz_dbl = dz.EvalDbl()
        self.assertEqual(dz_dbl.real,0.0)
        self.assertEqual(dz_dbl.imag,0.0)

    def test_pi(self):
        z = Pi()
        self.assertEqual(z.EvalDbl().real,3.14159265358979324)
        self.assertEqual(z.EvalDbl().imag, 0.0)

    def test_e(self):
        z = E()
        self.assertEqual(z.EvalDbl().real,2.71828182845904524)
        self.assertEqual(z.EvalDbl().imag, 0.0)



class VariableTest(unittest.TestCase):
    def test_eval_variable(self):
        v = Variable("x")
        v.set_current_value(4.56 + 8j)
        self.assertEqual(v.EvalDbl().real, 4.56)
        self.assertEqual(v.EvalDbl().imag, 8)

    def test_variable_deg(self):
        v = Variable("x")
        w = Variable("y")
        vg = VariableGroup()
        vg.append(v)
        vg.append(w)
        self.assertEqual(v.Degree(v),1)
        self.assertEqual(v.Degree(w),0)
        self.assertEqual(v.Degree(vg),1)
        self.assertTrue(v.IsHomogeneous)



class OperatorTest(unittest.TestCase):
    def test_sum_numbers(self):
        x = Float(5.97,4.22)
        y = Variable("x")
        y.set_current_value(4.56 - 8.3j)
        s = x+y;
        eval_s = s.EvalDbl()
        self.assertEqual(eval_s.real,5.97+4.56)
        self.assertEqual(eval_s.imag,4.22-8.3)

        s = y+5.97
        self.assertEqual(s.EvalDbl().real,4.56+5.97)
        self.assertEqual(s.EvalDbl().imag,-8.3)

        s = 5.97+y
        self.assertEqual(s.EvalDbl().real,4.56+5.97)
        self.assertEqual(s.EvalDbl().imag,-8.3)

        s = y+(5.97+4.22j)
        self.assertEqual(s.EvalDbl().real,4.56+5.97)
        self.assertEqual(s.EvalDbl().imag,-8.3+4.22)

        s = (5.97+4.22j)+y
        self.assertEqual(s.EvalDbl().real,4.56+5.97)
        self.assertEqual(s.EvalDbl().imag,-8.3+4.22)

        s = y+5
        self.assertEqual(s.EvalDbl().real,4.56+5)
        self.assertEqual(s.EvalDbl().imag,-8.3)

        s = 5+y
        self.assertEqual(s.EvalDbl().real,4.56+5)
        self.assertEqual(s.EvalDbl().imag,-8.3)

        z = mpfr_complex(5.97,4.22)
        s = y+z
        self.assertEqual(s.EvalDbl().real,4.56+5.97)
        self.assertEqual(s.EvalDbl().imag,-8.3+4.22)

        s = z+y
        self.assertEqual(s.EvalDbl().real,4.56+5.97)
        self.assertEqual(s.EvalDbl().imag,-8.3+4.22)

    def test_add_sub_mult_div(self):
        x = Variable("x"); y = Variable("y"); z = Variable("z")
        xdbl = 4.56 + 3.4j; ydbl = -3.2 - .72j; zdbl = 3.65 - 2.341j;
        x.set_current_value(xdbl)
        y.set_current_value(ydbl)
        z.set_current_value(zdbl)

        g = Function(x+y - z*y + x/(y*z) - (x*x*y+3.2)/z)

        gexact = xdbl+ydbl - zdbl*ydbl + xdbl/(ydbl*zdbl) - (xdbl*xdbl*ydbl+3.2)/zdbl

        self.assertTrue(g.EvalDbl().real/gexact.real - 1 < dbltol)
        self.assertTrue(g.EvalDbl().imag/gexact.imag - 1 < dbltol)



unittest.main()
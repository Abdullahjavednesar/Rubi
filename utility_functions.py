'''
Utility functions for Constraints in Rubi
'''
import math

from sympy.polys.polytools import degree
from sympy.functions.elementary.hyperbolic import acosh

def NonzeroQ(expr):
    return expr != 0

def FreeQ(nodes, var):
    return not any(expr.has(var) for expr in nodes)

def ZeroQ(expr):
    return expr == 0

def PositiveIntegerQ(var):
    return var.is_Integer() and var > 0

def NegativeIntegerQ(var):
    return var.is_Integer() and var < 0

def PositiveQ(var):
    return var > 0

def IntegerQ(var):
    return var.is_Integer()

def PosQ(var):
    return var > 0

def FracPart(var):
    return var - IntPart(var)

def IntPart(var):
    return math.floor(var)

def NegQ(var):
    return var < 0

def RationalQ(var):
    return var.is_rational_function()

def Subst(a, x, y):
    return a.subs(x, y)

def linearQ(expr, x):
    if degree(expr, gen=x) == 1:
        return True

def TogetherSimplify(expr):
    return

def Coefficient(u, var, n):
    return 

def RemoveContent():
    return

def Sqrt(a):
    return math.sqrt(a)

def ExpandIntegrand(expr, x):
    return 

def ArcCosh(a):
    return acosh(a)

def With():
    return

def Denominator():
    return

def Hypergeometric2F1():
    return

def TogetherSimplify():
    return

def IntLinearcQ():
    return

def ArcTan():
    return

def Not():
    return

def Simplify():
    return

def FractionalPart():
    return

def IntegerPart():
    return

def Simp():
    return

def Rt():
    return

def Log():
    return

def SumSimplerQ():
    return

# utility functions used in tests
def AppellF1():
    return

def Integrate():
    return

def hypergeom():
    return

def EllipticPi():
    return

def EllipticE():
    return

def EllipticF():
    return

def arctan():
    return

def arctanh():
    return

def arcsin():
    return

def arcsinh():
    return

def arccos():
    return

def arccosh():
    return

def arccsc():
    return

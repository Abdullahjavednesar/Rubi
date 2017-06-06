'''
Utility functions for Constraints in Rubi
'''
import math

def NonzeroQ(expr):
    return expr != 0

def FreeQ(nodes, var):
    return not any(expr.has(var) for expr in nodes)

def ZeroQ(expr):
    return expr == 0

def PositiveIntegerQ(var):
    return isinstance(var, int) and var > 0

def NegativeIntegerQ(var):
    return isinstance(var, int) and var < 0

def PositiveQ(var):
    return var > 0

def IntegerQ(var):
    return isinstance(var, int)

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

# todo
def Subst(expr):
    return 

def linearQ(var, x):
    return

def IntLinearcQ(*args, x):
    return

def SumSimplerQ(expr, var):
    return

def SimplerQ(expr1, expr2):
    return

def TogetherSimplify(expr):
    return

def Hypergeometric2F1():
    return

def Coefficient():
    return
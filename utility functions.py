'''
Utility functions for Constraints in Rubi
'''
def NonzeroQ(expr):
    return expr != 0

def FreeQ(nodes, var):
    return not any(expr.has(var) for expr in nodes)

def ZeroQ(expr):
	return expr = 0

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
    return var % 1


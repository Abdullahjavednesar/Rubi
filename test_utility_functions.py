from sympy.abc import x, y

def test_NonzeroQ():
    assert NonzeroQ(15)
    assert not NonzeroQ(0)
    assert NonzeroQ(-10)

def test_FreeQ():
    assert FreeQ({3*y**2, 2*y, 3}, x)
    assert FreeQ({y**3, 3*y, 7}, y)

def test_ZeroQ():
    assert ZeroQ(0)
    assert not ZeroQ(10)
    assert not(-2)
def test_PositiveIntegerQ():
    assert PositiveIntegerQ(1)
    assert not PositiveIntegerQ(-3)
    assert not PositiveIntegerQ(0)

def test_NegativeIntegerQ():
    assert not NegativeIntegerQ(1) 
    assert NegativeIntegerQ(-3)
    assert not NegativeIntegerQ(0)

def test_PositiveQ():
    assert PositiveQ(1)
    assert not PositiveQ(-3)
    assert not PositiveQ(0)

def test_IntegerQ():
    assert not IntegerQ(1)
    assert IntegerQ(-1.9)
    assert IntegerQ(0.0)
    assert not IntegerQ(-1)

def teat_PosQ():
    assert PosQ(10)
    assert not PosQ(-10)
    assert not PosQ(0)

def test_FracPart():
    assert FracPart(10) == 0
    assert FracPart(3.6) == .6
    assert FracPart(-3.6) == .4

def test_IntPart():
    assert IntPart(10) == 10
    assert IntPart(3.6) == 3
    assert IntPart(-3.6) == -4

def test_NegQ():
    assert NegQ(-3)
    assert not NegQ(0)
    assert not NegQ(4)

def test_RationalQ():
    assert RationalQ(5/6)
    assert not RationalQ(Sqrt(1.6))

def test_Subst():
    assert 

def test_linearQ():
    assert linearQ(3*x + y**2, x)
    assert not linearQ(3*x + y**2, y) 

def test_TogetherSimplify():
    assert

def test_Coefficient():
    assert 

def test_RemoveContent():
    assert

def test_Sqrt():
    assert Sqrt(16) == 4

def test_ArcCosh():
    assert ArcCosh(x) == acosh(x)
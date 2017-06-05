from matchpy import Wildcard, Pattern, ReplacementRule, is_match, replace_all, ManyToOneReplacer

from operation import Int, Mul, Add, Pow, Log
from symbol import VariableSymbol, ConstantSymbol
from constraint import FreeQ, NonzeroQ


a, b, c, d, e, f, g, h, x = map(VariableSymbol, 'abcdefghx')
n, m = map(VariableSymbol, 'nm')
a_, b_, c_, d_, e_, f_, g_, h_ = map(Wildcard.dot, 'abcdefgh')
n_, m_ = map(Wildcard.dot, 'nm')

pattern1 = Pattern(Int(1/x, x))
rule1 = ReplacementRule(pattern1, lambda x: Log(x))

pattern2 = Pattern(Int(x**m, x), FreeQ(m, x) and NonzeroQ(m + 1))
rule2 = ReplacementRule(pattern2, lambda m, x: x**(m + 1)/(m + 1))

pattern3 = Pattern(Int(1/(a + b*x), x), FreeQ((a, b), x))
rule3 = ReplacementRule(pattern3, lambda a, b, x: Log(RemoveContent(a + b*x, x))/b)

pattern4 = Pattern(Int((a + b*x)**m, x), FreeQ((a, b, m), x) and NonzeroQ(m + 1))
rule4 = ReplacementRule(pattern4, lambda a, b, m, x: (a + b*x)**(m + 1)/(b*(m + 1)))

pattern5 = Pattern(Int((a + b*u)**m, x), FreeQ((a, b, m), x) and LinearQ(u, x) and NonzeroQ(u - x))
rule5 = ReplacementRule(pattern5, lambda a, b, m, u, x: 1/Coefficient(u, x, 1)*Subst(Pattern(Int((a + b*x)**m, x), x, u)))

# 1.2 (a + b x)**m (c + d x)**n)

pattern6 = Pattern(Int(1/((a + b*x)*(c + d*x)), x), FreeQ((a, b, c, d), x) and ZeroQ(b*c + a*d))
rule6 = ReplacementRule(pattern6, lambda a, b, c, d, x: Pattern(Int(1/(a*c + b*d*x**2), x)))

pattern7 = Pattern(Int(1/((a + b*x)*(c + d*x)), x), FreeQ((a, b, c, d), x) and NonzeroQ(b*c - a*d))
rule7 = ReplacementRule(pattern7, lambda a, b, c, d, x: b/(b*c - a*d)*Pattern(Int(1/(a + b*x), x)  -  d/(b*c - a*d)*Pattern(Int(1/(c + d*x), x))))

pattern8 = Pattern(Int((a + b*x)**m*(c + d*x)**n, x), FreeQ((a, b, c, d, m, n), x) and NonzeroQ(b*c - a*d) and ZeroQ(m + n + 2) and NonzeroQ(m + 1))
rule8 = ReplacementRule(pattern8, lambda a, b, c, d, n, m: (a + b*x)**(m + 1)*(c + d*x)**(n + 1)/((b*c - a*d)*(m + 1)))

pattern9= Pattern(Int((a + b*x)**m*(c + d*x)**m, x), FreeQ((a, b, c, d), x) and ZeroQ(b*c + a*d) and PositiveIntegerQ(m + 1/2))
rule9 = ReplacementRule(pattern9, lambda a, b, c, d, m: x*(a + b*x)**m*(c + d*x)**m/(2*m + 1)  +  2*a*c*m/(2*m + 1)*Pattern(Int((a + b*x)**(m - 1)*(c + d*x)**(m - 1), x)))

pattern10 = Pattern(Int(1/((a + b*x)**(3/2)*(c + d*x)**(3/2)), x), FreeQ((a, b, c, d), x) and ZeroQ(b*c + a*d))
rule10 = ReplacementRule(pattern10, lambda a, b, c, d: x/(a*c*Sqrt(a + b*x)*Sqrt(c + d*x)))

pattern11 = Pattern(Int((a + b*x)**m*(c + d*x)**m, x), FreeQ((a, b, c, d), x) and ZeroQ(b*c + a*d) and NegativeIntegerQ(m + 3/2))
rule11 = ReplacementRule(pattern11, lambda a, b, c, d, m:  - x*(a + b*x)**(m + 1)*(c + d*x)**(m + 1)/(2*a*c*(m + 1)) +  
  (2*m + 3)/(2*a*c*(m + 1))*Pattern(Int((a + b*x)**(m + 1)*(c + d*x)**(m + 1), x)))

pattern12 = Pattern(Int((a + b*x)**m*(c + d*x)**m, x), FreeQ((a, b, c, d, m), x) and ZeroQ(b*c + a*d) and (IntegerQ(m) or PositiveQ(a) and PositiveQ(c)))
rule12 = ReplacementRule(pattern12, lambda a, b, c, d, m: Pattern(Int((a*c + b*d*x**2)**m, x)))

pattern13 = Pattern(Int(1/(Sqrt(a + b*x)*Sqrt(c + d*x)), x), FreeQ((a, b, c, d), x) and ZeroQ(b*c + a*d) and PositiveQ(a) and ZeroQ(a + c))
rule13 = ReplacementRule(pattern13, lambda a, b, c, d: ArcCosh(b*x/a)/b)

pattern14 = Pattern(Int(1/(Sqrt(a + b*x)*Sqrt(c + d*x)), x),  FreeQ((a, b, c, d), x) and ZeroQ(b*c + a*d))
rule14 = ReplacementRule(pattern14, lambda a, b, c, d: 2*Subst(Pattern(Int(1/(b - d*x**2), x), x, Sqrt(a + b*x)/Sqrt(c + d*x))))

pattern15 = Pattern(Int((a + b*x)**m*(c + d*x)**m, x),  FreeQ((a, b, c, d, m), x) and ZeroQ(b*c + a*d) and not(IntegerQ(2*m)))
rule15 = ReplacementRule(pattern15, lambda a, b, c, d, m: (a + b*x)**FracPart(m)*(c + d*x)**FracPart(m)/(a*c + b*d*x**2)**FracPart(m)*Pattern(Int((a*c + b*d*x**2)**m, x)))

pattern16 = Pattern(Int(1/((a + b*x)**(5/4)*(c + d*x)**(1/4)), x), FreeQ((a, b, c, d), x) and ZeroQ(b*c + a*d) and PosQ(b*d/(a*c)))
rule16 = ReplacementRule(pattern16, lambda a, b, c, d:  - 2/(b*(a + b*x)**(1/4)*(c + d*x)**(1/4)) +  (b*c - a*d)/(2*b)*Pattern(Int(1/((a + b*x)**(5/4)*(c + d*x)**(5/4)), x)))

pattern17 = Pattern(Int(1/((a + b*x)**(9/4)*(c + d*x)**(1/4)), x), FreeQ((a, b, c, d), x) and ZeroQ(b*c + a*d) and PosQ(b*d/(a*c)))
rule17 = ReplacementRule(pattern17, lambda a, b, c, d:  - 4/(5*b*(a + b*x)**(5/4)*(c + d*x)**(1/4)) -  d/(5*b)*Pattern(Int(1/((a + b*x)**(5/4)*(c + d*x)**(5/4)), x)))

pattern18 = Pattern(Int((a + b*x)**m*(c + d*x)**n, x), FreeQ((a, b, c, d), x) and ZeroQ(b*c + a*d) and IntegerQ(m + 1/2) and IntegerQ(n + 1/2) and 0<m<n)
rule18 = ReplacementRule(pattern18, lambda a, b, c, d, n, m: (a + b*x)**(m + 1)*(c + d*x)**n/(b*(m + n + 1)) + 2*c*n/(m + n + 1)*Pattern(Int((a + b*x)**m*(c + d*x)**(n - 1), x)))

pattern19 = Pattern(Int((a + b*x)**m*(c + d*x)**n, x), FreeQ((a, b, c, d), x) and ZeroQ(b*c + a*d) and IntegerQ(m + 1/2) and IntegerQ(n + 1/2) and m<n<0)
rule19 = ReplacementRule(pattern19, lambda a, b, c, d, n, m:  - (a + b*x)**(m + 1)*(c + d*x)**(n + 1)/(2*a*d*(m + 1)) + (m + n + 2)/(2*a*(m + 1))*Pattern(Int((a + b*x)**(m + 1)*(c + d*x)**n, x)))

pattern20 = Pattern(Int((a + b*x)**m*(c + d*x)**n, x), FreeQ((a, b, c, d, n), x) and NonzeroQ(b*c - a*d) and PositiveIntegerQ(m) and (not(IntegerQ(n))or ZeroQ(c) and 7*m + 4*n<=0 or 9*m + 5*(n + 1)<0 or m + n + 2>0))
rule20 = ReplacementRule(pattern20, lambda a, b, c, d, n, m: Pattern(Int(ExpandIntegrand((a + b*x)**m*(c + d*x)**n, x), x)))

pattern21 = Pattern(Int((a + b*x)**m*(c + d*x)**n, x), FreeQ((a, b, c, d, n), x) and NonzeroQ(b*c - a*d) and NegativeIntegerQ(m) and IntegerQ(n) and not(PositiveIntegerQ(n) and m + n + 2<0))
rule21 = ReplacementRule(pattern21, lambda a, b, c, d, n, m: Pattern(Int(ExpandIntegrand((a + b*x)**m*(c + d*x)**n, x), x)))

pattern22 = Pattern(Int((c + d*x)**n/(a + b*x), x), FreeQ((a, b, c, d), x) and NonzeroQ(b*c - a*d) and RationalQ(n) and n>0)
rule22 = ReplacementRule(pattern22, lambda a, b, c, d, n: (c + d*x)**n/(b*n)  +  (b*c - a*d)/b*Pattern(Int((c + d*x)**(n - 1)/(a + b*x), x)))

pattern23 = Pattern(Int((c + d*x)**n/(a + b*x), x), FreeQ((a, b, c, d), x) and NonzeroQ(b*c - a*d) and RationalQ(n) and n< - 1)
rule23 = ReplacementRule(pattern23, lambda a, b, c, d, n:  - (c + d*x)**(n + 1)/((n + 1)*(b*c - a*d)) + b*(n + 1)/((n + 1)*(b*c - a*d))*Pattern(Int((c + d*x)**(n + 1)/(a + b*x), x)))

pattern24 = Pattern(Int(1/((a + b*x)*(c + d*x)**(1/3)), x), FreeQ((a, b, c, d), x) and PosQ((b*c - a*d)/b))
rule24 = ReplacementRule(pattern24, lambda a, b, c, d: With((q=Rt((b*c - a*d)/b, 3)), 
      - Log(RemoveContent(a + b*x, x))/(2*b*q)  -  3/(2*b*q)*
      Subst(Pattern(Int(1/(q - x), x), x, (c + d*x)**(1/3)) + 3/(2*b)*Subst(Pattern(Int(1/(q**2 + q*x + x**2), x), x, (c + d*x)**(1/3))))))

pattern25 = Pattern(Int(1/((a + b*x)*(c + d*x)**(1/3)), x), FreeQ((a, b, c, d), x) and NegQ((b*c - a*d)/b))
rule25 = ReplacementRule(pattern25, lambda a, b, c, d: With((q=Rt( - (b*c - a*d)/b, 3)), Log(RemoveContent(a + b*x, x))/(2*b*q) - 3/(2*b*q)*Subst(Pattern(Int(1/(q + x), x), x, (c + d*x)**(1/3)) + 
      3/(2*b)*Subst(Pattern(Int(1/(q**2 - q*x + x**2), x), x, (c + d*x)**(1/3))))))

pattern26 = Pattern(Int(1/((a + b*x)*(c + d*x)**(2/3)), x), FreeQ((a, b, c, d), x) and PosQ((b*c - a*d)/b))
rule26 = ReplacementRule(pattern26, lambda a, b, c, d: With((q=Rt((b*c - a*d)/b, 3)), 
   - Log(RemoveContent(a + b*x, x))/(2*b*q**2)  -  
  3/(2*b*q**2)*Subst(Pattern(Int(1/(q - x), x), x, (c + d*x)**(1/3)) -  
  3/(2*b*q)*Subst(Pattern(Int(1/(q**2 + q*x + x**2), x), x, (c + d*x)**(1/3))))))

pattern27 = Pattern(Int(1/((a + b*x)*(c + d*x)**(2/3)), x), FreeQ((a, b, c, d), x) and NegQ((b*c - a*d)/b))
rule27 = ReplacementRule(pattern27, lambda a, b, c, d: With((q=Rt( - (b*c - a*d)/b, 3)), 
   - Log(RemoveContent(a + b*x, x))/(2*b*q**2)  +  
  3/(2*b*q**2)*Subst(Pattern(Int(1/(q + x), x), x, (c + d*x)**(1/3)) +  
  3/(2*b*q)*Subst(Pattern(Int(1/(q**2 - q*x + x**2), x), x, (c + d*x)**(1/3))))))

pattern28 = Pattern(Int((c + d*x)**n/(a + b*x), x), FreeQ((a, b, c, d), x) and NonzeroQ(b*c - a*d) and RationalQ(n) and  - 1<n<0)
rule28 = ReplacementRule(pattern28, lambda a, b, c, d, n: With((p=Denominator(n)), 
  p*Subst(Pattern(Int(x**(p*(n + 1) - 1)/(a*d - b*c + b*x**p), x), x, (c + d*x)**(1/p)))))

pattern29 = Pattern(Int((c + d*x)**n/x, x), FreeQ((c, d, n), x) and not(IntegerQ(n)))
rule29 = ReplacementRule(pattern29, lambda a, b, c, d, n:  - (c + d*x)**(n + 1)/(c*(n + 1))*Hypergeometric2F1(1, n + 1, n + 2, 1 + d*x/c))

pattern30 = Pattern(Int((c + d*x)**n/(a + b*x), x),  FreeQ((a, b, c, d, n), x) and NonzeroQ(b*c - a*d) and not(IntegerQ(n)))
rule30 = ReplacementRule(pattern30, lambda a, b, c, d, n:  - (c + d*x)**(n + 1)/((n + 1)*(b*c - a*d))*Hypergeometric2F1(1, n + 1, n + 2, TogetherSimplify(b*(c + d*x)/(b*c - a*d))))

## rest rules 1.2 algebraic

patternx = Pattern(Int((a_ + b_*x_)**m_*(c_ + d_*x_)**n_, x),, FreeQ((a, b, c, d), x) and NonzeroQ(b*c - a*d) and RationalQ(m, n) and m< - 1 and n>0 and Not(IntegerQ(n) and Not(IntegerQ(m))) and 
  Not(IntegerQ(m + n) and m + n + 2<=0 and (FractionQ(m) or 2*n + m + 1>=0)) and IntLinearcQ(a, b, c, d, m, n, x))
rulex = ReplacementRule(patternx, lambda a, b, c, d, m, n, x: (a + b*x)**(m + 1)*(c + d*x)**n/(b*(m + 1))  -  
  d*n/(b*(m + 1))* Pattern(Int((a + b*x)**(m + 1)*(c + d*x)**(n - 1), x)))

patternx = Pattern(Int((a_ + b_*x_)**m_*(c_ + d_*x_)**n_, x), FreeQ((a, b, c, d), x) and NonzeroQ(b*c - a*d) and RationalQ(m, n) and m< - 1 and 
  Not(n< - 1 and (ZeroQ(a) or NonzeroQ(c) and m<n and IntegerQ(n))) and IntLinearcQ(a, b, c, d, m, n, x))
rulex = ReplacementRule(patternx, lambda a, b, c, d, m, n, x: (a + b*x)**(m + 1)*(c + d*x)**(n + 1)/((b*c - a*d)*(m + 1))  -  
  d*(m + n + 2)/((b*c - a*d)*(m + 1))* Pattern(Int((a + b*x)**(m + 1)*(c + d*x)**n, x)))

patternx = Pattern(Int((a_ + b_*x_)**m_*(c_ + d_*x_)**n_, x), FreeQ((a, b, c, d), x) and NonzeroQ(b*c - a*d) and RationalQ(m, n) and n>0 and m + n + 1!=0 and 
  Not(PositiveIntegerQ(m) and (Not(IntegerQ(n)) or 0<m<n)) and 
  Not(IntegerQ(m + n) and m + n + 2<0) and IntLinearcQ(a, b, c, d, m, n, x))
rulex = ReplacementRule(patternx, lambda a, b, c, d, m, n, x: (a + b*x)**(m + 1)*(c + d*x)**n/(b*(m + n + 1))  +  
  n*(b*c - a*d)/(b*(m + n + 1))* Pattern(Int((a + b*x)**m*(c + d*x)**(n - 1), x)))

patternx = Pattern(Int(1/(Sqrt(a_ + b_*x_)*Sqrt(c_ + d_*x_)), x), FreeQ((a, b, c, d), x) and ZeroQ(b + d) and PositiveQ(a + c))
rulex = ReplacementRule(patternx, lambda a, b, c, d, m, n, x: Pattern(Int(1/Sqrt(a*c - b*(a - c)*x - b**2*x**2), x)))

patternx = Pattern(Int(1/(Sqrt(a_ + b_*x_)*Sqrt(c_ + d_*x_)), x), FreeQ((a, b, c, d), x) and PositiveQ(b*c - a*d) and PositiveQ(b))
rulex = ReplacementRule(patternx, lambda a, b, c, d, m, n, x: 2/Sqrt(b)*SubstPattern(Int(1/Sqrt(b*c - a*d + d*x**2), x), x, Sqrt(a + b*x)))

patternx = Pattern(Int(1/(Sqrt(a_ + b_*x_)*Sqrt(c_ + d_*x_)), x), FreeQ((a, b, c, d), x) and NonzeroQ(b*c - a*d) and ZeroQ(b - d))
rulex = ReplacementRule(patternx, lambda a, b, c, d, m, n, x: 2/b*SubstPattern(Int(1/Sqrt(c - a + x**2), x), x, Sqrt(a + b*x)))

patternx = Pattern(Int(1/(Sqrt(a_ + b_*x_)*Sqrt(c_ + d_*x_)), x), FreeQ((a, b, c, d), x) and NonzeroQ(b*c - a*d))
rulex = ReplacementRule(patternx, lambda a, b, c, d, m, n, x: 2*SubstPattern(Int(1/(b - d*x**2), x), x, Sqrt(a + b*x)/Sqrt(c + d*x)))

patternx = Pattern(Int((a_ + b_*x_)**m_*(c_ + d_*x_)**m_, x), FreeQ((a, b, c, d), x) and NonzeroQ(b*c - a*d) and RationalQ(m) and  - 1<m<0 and 3<=Denominator(m)<=4)
rulex = ReplacementRule(patternx, lambda a, b, c, d, m, n, x: (a + b*x)**m*(c + d*x)**m/(a*c + (b*c + a*d)*x + b*d*x**2)**m* Pattern(Int((a*c + (b*c + a*d)*x + b*d*x**2)**m, x)))

patternx = Pattern(Int(1/((a_ + b_*x_)**(1/3)*(c_ + d_*x_)**(2/3)), x), FreeQ((a, b, c, d), x) and NonzeroQ(b*c - a*d) and PosQ(d/b)) 
rulex = ReplacementRule(patternx, lambda a, b, c, d, m, n, x: With((q=Rt(d/b, 3)), 
   - Sqrt(3)*q/d*ArcTan(2*q*(a + b*x)**(1/3)/(Sqrt(3)*(c + d*x)**(1/3)) + 1/Sqrt(3))  -  
  q/(2*d)*Log(c + d*x)  -  
  3*q/(2*d)*Log(q*(a + b*x)**(1/3)/(c + d*x)**(1/3) - 1)))

patternx = Pattern(Int(1/((a_ + b_*x_)**(1/3)*(c_ + d_*x_)**(2/3)), x), FreeQ((a, b, c, d), x) and NonzeroQ(b*c - a*d) and NegQ(d/b))
rulex = ReplacementRule(patternx, lambda a, b, c, d, m, n, x: With((q=Rt( - d/b, 3)), 
  Sqrt(3)*q/d*ArcTan(1/Sqrt(3) - 2*q*(a + b*x)**(1/3)/(Sqrt(3)*(c + d*x)**(1/3)))  +  
  q/(2*d)*Log(c + d*x)  +  
  3*q/(2*d)*Log(q*(a + b*x)**(1/3)/(c + d*x)**(1/3) + 1)))

patternx = Pattern(Int((a_ + b_*x_)**m_*(c_ + d_*x_)**n_, x), FreeQ((a, b, c, d), x) and NonzeroQ(b*c - a*d) and RationalQ(m, n) and  - 1<m<0 and m + n + 1==0)
rulex = ReplacementRule(patternx, lambda a, b, c, d, m, n, x: With((p=Denominator(m)), 
  p*SubstPattern(Int(x**(p*(m + 1) - 1)/(b - d*x**p), x), x, (a + b*x)**(1/p)/(c + d*x)**(1/p))))

patternx = Pattern(Int((a_ + b_*x_)**m_*(c_ + d_*x_)**n_, x),, FreeQ((a, b, c, d), x) and NonzeroQ(b*c - a*d) and RationalQ(m, n) and  - 1<m<0 and  - 1<n<0 and Denominator(n)<=Denominator(m) and 
  IntLinearcQ(a, b, c, d, m, n, x))
rulex = ReplacementRule(patternx, lambda a, b, c, d, m, n, x: With((p=Denominator(m)), 
  p/b*SubstPattern(Int(x**(p*(m + 1) - 1)*(c - a*d/b + d*x**p/b)**n, x), x, (a + b*x)**(1/p))))

patternx = Pattern(Int((a_ + b_*x_)**m_*(c_ + d_*x_)**n_, x), FreeQ((a, b, c, d, m, n), x) and NonzeroQ(b*c - a*d) and NegativeIntegerQ(Simplify(m + n + 2)) and NonzeroQ(m + 1) and 
  (SumSimplerQ(m, 1) or Not(SumSimplerQ(n, 1))))
rulex = ReplacementRule(patternx, lambda a, b, c, d, m, n, x: (a + b*x)**(m + 1)*(c + d*x)**(n + 1)/((b*c - a*d)*(m + 1))  -  
  d*Simplify(m + n + 2)/((b*c - a*d)*(m + 1))* Pattern(Int((a + b*x)**Simplify(m + 1)*(c + d*x)**n, x)))

patternx = Pattern(Int((b_*x_)**m_*(c_ + d_*x_)**n_, x), FreeQ((b, c, d, m, n), x) and Not(IntegerQ(m)) and (IntegerQ(n) or PositiveQ(c) and Not(ZeroQ(n + 1/2) and ZeroQ(c**2 - d**2) and PositiveQ( - d/(b*c)))))
rulex = ReplacementRule(patternx, lambda a, b, c, d, m, n, x: c**n*(b*x)**(m + 1)/(b*(m + 1))*Hypergeometric2F1( - n, m + 1, m + 2,  - d*x/c)) 


patternx = Pattern(Int((b_*x_)**m_*(c_ + d_*x_)**n_, x), FreeQ((b, c, d, m, n), x) and Not(IntegerQ(n)) and (IntegerQ(m) or PositiveQ( - d/(b*c)))) 
rulex = ReplacementRule(patternx, lambda a, b, c, d, m, n, x: (c + d*x)**(n + 1)/(d*(n + 1)*( - d/(b*c))**m)*Hypergeometric2F1( - m, n + 1, n + 2, 1 + d*x/c)) 


patternx = Pattern(Int((b_*x_)**m_*(c_ + d_*x_)**n_, x), FreeQ((b, c, d, m, n), x) and Not(IntegerQ(m)) and Not(IntegerQ(n)) and Not(PositiveQ(c)) and Not(PositiveQ( - d/(b*c))) and 
  (RationalQ(m) and Not(ZeroQ(n + 1/2) and ZeroQ(c**2 - d**2)) or Not(RationalQ(n))))
rulex = ReplacementRule(patternx, lambda a, b, c, d, m, n, x: c**IntPart(n)*(c + d*x)**FracPart(n)/(1 + d*x/c)**FracPart(n)* Pattern(Int((b*x)**m*(1 + d*x/c)**n, x)))

patternx = Pattern(Int((b_*x_)**m_*(c_ + d_*x_)**n_, x), FreeQ((b, c, d, m, n), x) and Not(IntegerQ(m)) and Not(IntegerQ(n)) and Not(PositiveQ(c)) and Not(PositiveQ( - d/(b*c))))
rulex = ReplacementRule(patternx, lambda a, b, c, d, m, n, x: ( - b*c/d)**IntPart(m)*(b*x)**FracPart(m)/( - d*x/c)**FracPart(m)* Pattern(Int(( - d*x/c)**m*(c + d*x)**n, x)))

patternx = Pattern(Int((a_ + b_*x_)**m_*(c_ + d_*x_)**n_, x), FreeQ((a, b, c, d, m), x) and NonzeroQ(b*c - a*d) and Not(IntegerQ(m)) and IntegerQ(n))
rulex = ReplacementRule(patternx, lambda a, b, c, d, m, n, x: (b*c - a*d)**n*(a + b*x)**(m + 1)/(b**(n + 1)*(m + 1))*Hypergeometric2F1( - n, m + 1, m + 2,  - d*(a + b*x)/(b*c - a*d)))

patternx = Pattern(Int((a_ + b_*x_)**m_*(c_ + d_*x_)**n_, x), FreeQ((a, b, c, d, m, n), x) and NonzeroQ(b*c - a*d) and Not(IntegerQ(m)) and Not(IntegerQ(n)) and PositiveQ(b/(b*c - a*d)) and 
  (RationalQ(m) or Not(RationalQ(n) and PositiveQ( - d/(b*c - a*d)))))
rulex = ReplacementRule(patternx, lambda a, b, c, d, m, n, x: (a + b*x)**(m + 1)/(b*(m + 1)*(b/(b*c - a*d))**n)*Hypergeometric2F1( - n, m + 1, m + 2,  - d*(a + b*x)/(b*c - a*d)))

patternx = Pattern(Int((a_ + b_*x_)**m_*(c_ + d_*x_)**n_, x), 
  FreeQ((a, b, c, d, m, n), x) and NonzeroQ(b*c - a*d) and Not(IntegerQ(m)) and Not(IntegerQ(n)) and (RationalQ(m) or Not(SimplerQ(n + 1, m + 1))))
rulex = ReplacementRule(patternx, lambda a, b, c, d, m, n, x: (c + d*x)**FracPart(n)/((b/(b*c - a*d))**IntPart(n)*(b*(c + d*x)/(b*c - a*d))**FracPart(n))*
    Pattern(Int((a + b*x)**m*(b*c/(b*c - a*d) + b*d*x/(b*c - a*d))**n, x)))

patternx = Pattern(Int((a_ + b_*u_)**m_*(c_ + d_*u_)**n_, x), FreeQ((a, b, c, d, m, n), x) and LinearQ(u, x) and NonzeroQ(Coefficient(u, x, 0)))
rulex = ReplacementRule(patternx, lambda a, b, c, d, m, n, x: 1/Coefficient(u, x, 1)*SubstPattern(Int((a + b*x)**m*(c + d*x)**n, x), x, u))

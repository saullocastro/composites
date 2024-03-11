print("""
# shear buckling
# section 6.4
# page 138, Eq. 6.30
# differentiation used in the Newton-Raphson approach herein implemented
""")
import sympy
from sympy import pi, sin, tan

sympy.var("D11, D12, D22, D66, AR, alpha")
expr = (3*D11*AR**4*tan(alpha)**4
        + (6*D11*AR**2 + 2*(D12 + 2*D66)*AR**4)*tan(alpha)**2
        - (D11 + 2*(D12 + 2*D66)*AR**2 + D22*AR**4))
dexpr_dalpha = expr.diff(alpha)
print("dexpr_dalpha =", dexpr_dalpha)
print()
print("""
# combined compression-shear buckling
# section 6.5
# Ritz method therein used
""")
sympy.var("D11, D12, D16, D22, D26, D66")
sympy.var("Nxx, Nxy, k")
#Nxy = k*Nxx
sympy.var("w1, w2")
sympy.var("a, b, x, y", positive=True)
w = w1*sin(pi*x/a)*sin(pi*y/b) + w2*sin(2*pi*x/a)*sin(2*pi*y/b)
#Eq. 6.32
integrand_Uc = (D11*(w.diff(x, x))**2
            + 2*D12*w.diff(x, x)*w.diff(y, y)
              + D22*w.diff(y, y)**2
            + 4*D66*w.diff(x, y)**2
            + 4*D16*w.diff(x, x)*w.diff(x, y)
            + 4*D26*w.diff(y, y)*w.diff(x, y)
               )/2
integrand_Vc = (
                Nxx*w.diff(x)**2/2
              + Nxy*w.diff(x)*w.diff(y)
               )
PIc = sympy.integrate(integrand_Uc + integrand_Vc, (x, 0, a), (y, 0, b))
eq1 = PIc.diff(w1)
eq2 = PIc.diff(w2)
#sympy.solve((eq1, eq2), (w1, w2))
#
# eq1.collect((w1, w2))
# eq1, terms for w1
print(eq1.subs(dict(w2=0)).collect((w1,Nxx,Nxy)))
print(eq1.subs(dict(w1=0)).collect((w2,Nxx,Nxy)))
a11 = pi**4*D11*b/(4*a**3) + pi**4*D12/(2*a*b) + pi**4*D22*a/(4*b**3) + pi**4*D66/(a*b)
a12 = -160*pi**2*D16/(9*a**2) - 160*pi**2*D26/(9*b**2)
b11 = pi**2*Nxx*b/(4*a)
b12 = -32*Nxy/9

print(eq2.subs(dict(w2=0)).collect((w1,Nxx,Nxy)))
print(eq2.subs(dict(w1=0)).collect((w2,Nxx,Nxy)))
a21 = -160*pi**2*D16/(9*a**2) - 160*pi**2*D26/(9*b**2)
a22 = 4*pi**4*D11*b/a**3 + 8*pi**4*D12/(a*b) + 4*pi**4*D22*a/b**3 + 16*pi**4*D66/(a*b)
b21 = -32*Nxy/9
b22 = pi**2*Nxx*b/a

print('a11 =', a11)
print('a12 =', a12)
print('a21 =', a21)
print('a22 =', a22)
#
a11, a12, a21, a22 = sympy.var('a11, a12, a21, a22')
A = sympy.Matrix([[a11, a12],
                  [a21, a22]])
B = sympy.Matrix([[b11, b12],
                  [b21, b22]])
cp = sympy.det(A - B)
N0 = sympy.solve(cp, Nxx)
print('________')
print()
print('Solution')
print('________')
print(N0[0].expand().simplify())
print(N0[1].expand().simplify())

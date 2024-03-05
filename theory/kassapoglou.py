# shear buckling, Section 6.3, differentiation used in the Newton-Raphson
# approach therein implemented
import sympy

sympy.var('D11, D12, D22, D66, AR, alpha')
expr = (3*D11*AR**4*sympy.tan(alpha)**4
        + (6*D11*AR**2 + 2*(D12 + 2*D66)*AR**4)*sympy.tan(alpha)**2
        - (D11 + 2*(D12 + 2*D66)*AR**2 + D22*AR**4))
dexpr_dalpha = expr.diff(alpha)
print('dexpr_dalpha =', dexpr_dalpha)

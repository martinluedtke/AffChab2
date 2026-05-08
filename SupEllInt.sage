r"""
This file defines a class ``SuperEllipticExampleCurve`` for curves of the form
  y^3 = x^3 + a x^2 + x
and a function ``coleman_integrals_on_basis`` for computing the Coleman integrals
of the differential forms ω₁ = dx/y^2,  ω₂ = x dx/y^2,  ω₃ = dx/y.
The method is described in [LL26], Section 5.2.

REFERENCES:
- [LL26] Marius Leonhardt, Martin Lüdtke, "Affine Chabauty II"

AUTHORS:
- Marius Leonhardt
- Martin Lüdtke

"""

import sage.schemes.curves.projective_curve as plane_curve

class SuperEllipticExampleCurve(plane_curve.ProjectivePlaneCurve):
    r"""
    Plane projective curve of the form y^3 = x^3 + a x^2 + x.

    EXAMPLES::

        sage: SuperEllipticExampleCurve(QQ, 3)
        Projective Plane Curve over Rational Field defined by -X^3 + Y^3 - 3*X^2*Z - X*Z^2

    ::

        sage: SuperEllipticExampleCurve(QQ, 2)
        Traceback (most recent call last):
        ...
        ValueError: Curve y^3 = x^3 + ax^2 + x is singular for a = ±2

    """

    def __init__(self, K, a):
        if a == 2 or a == -2:
            raise ValueError("Curve y^3 = x^3 + ax^2 + x is singular for a = ±2")
        self.param = a
        PP.<X,Y,Z> = ProjectiveSpace(2, K)
        F = Y^3 - X^3 - a*X^2*Z - X*Z^2
        plane_curve.ProjectivePlaneCurve.__init__(self,PP,F)
    
    def coleman_integrals_on_basis(self,P,Q):
        r"""
        Compute the p-adic integrals from ``P`` to ``Q``  of the log differentials 
           ω₁ = dx/y^2,  ω₂ = x dx/y^2,  ω₃ = dx/y.
        
        EXAMPLES::

            sage: X = SuperEllipticExampleCurve(Qp(13,prec=10),3)
            sage: X.coleman_integrals_on_basis(X(0,0), X(-8/3,-2/3))
            [2*13 + 6*13^2 + 12*13^3 + 3*13^4 + 2*13^5 + 9*13^6 + 5*13^7 + 12*13^8 + 2*13^9 + O(13^10),
             7*13 + 3*13^2 + 3*13^3 + 12*13^6 + 6*13^7 + 5*13^8 + 13^9 + O(13^10),
             8*13 + 4*13^2 + 11*13^3 + 8*13^4 + 12*13^5 + 10*13^6 + 3*13^7 + 4*13^8 + 11*13^9 + O(13^10)]

        """
        K = self.base_ring()
        a = self.param
        p = K.prime()

        
        if K(a+2).valuation() > 0 or K(a-2).valuation() > 0:
            raise ValueError(f"Prime {p} is a prime of bad reduction")
        
        if not K.has_root_of_unity(3):
            raise ValueError("3rd roots of unity need to be contained in base field")
        ζ3 = K.primitive_root_of_unity(3)
        
        for R in [P,Q]:
            if R[2] == 0 or R[0].valuation() < 0 or R[1].valuation() < 0:
                raise ValueError("Endpoints cannot lie in infinite residue discs!")
            if a != 0 and (R[0] + 1/a).valuation() > 0:
                raise ValueError("Endpoints cannot lie in residue discs with x == -1/a mod p!")

        E = EllipticCurve(K, [0,0,0,0,a^2/4-1])
        to_E = lambda R: E(R[1], 1 + a/2*R[0], R[0])
        int_omega1 = -3*E.coleman_integrals_on_basis(to_E(P),to_E(Q))[0]

        # compute ∫_P^Q as ∫_P₀^Q - ∫_P₀^P with P₀=(0,0)
        int_omega2 = 0
        int_omega3 = 0
        for (sign,R) in [(+1,Q),(-1,P)]:
            if R[1] != 0:  # ∫_P₀^R is zero when y(R)=0
                w = lambda S: S[0]/S[1]
                # integrals of symmetric parts of ω₂ and ω₃ are computed on ℙ¹
                int_omega2_symm = -1/2 * sum(ζ^2 * (w(R)-ζ).log(0) for ζ in [1,ζ3,ζ3^2])
                int_omega3_symm = -1/2 * sum(ζ * (w(R)-ζ).log(0) for ζ in [1,ζ3,ζ3^2])

                int_omega2_antisymm = 0
                int_omega3_antisymm = 0
                if a != 0:  # there are no antisymmetric parts when a = 0
                    _.<x> = PolynomialRing(K)
                    for ζ in [1,ζ3,ζ3^2]:
                        # integrals of antisymmetric parts are computed on auxiliary curves Xζ
                        # where integrands become s ds/(2t)
                        T = x^4 + 12/a^2*ζ^2*x^3 + 12/a^2*ζ*x^2 + 4/a^2*x
                        Xζ = HyperellipticCurve(T, check_squarefree = false)
                        to_Xζ = lambda S: Xζ(w(S)/(1-w(S)*ζ), -2*(1/S[0]+a/2)/a * w(S)^2/(1-w(S)*ζ)^2) if S[0] != 0 else Xζ(0,0)
                        Xζ_int = Xζ.coleman_integrals_on_basis(Xζ(0,0),to_Xζ(R))[1]
                        int_omega2_antisymm += -ζ * Xζ_int
                        int_omega3_antisymm += -ζ^2 * Xζ_int

                int_omega2 += sign * (int_omega2_symm + int_omega2_antisymm)
                int_omega3 += sign * (int_omega3_symm + int_omega3_antisymm)
        
        return [int_omega1, int_omega2, int_omega3]





import sage.schemes.curves.projective_curve as plane_curve

class SuperEllipticExampleCurve(plane_curve.ProjectivePlaneCurve):
    # A superelliptic curve of the form y^3 = x^3 + a*x^2 + x

    def __init__(self, K, a):
        if a == 2 or a == -2:
            raise ValueError("Curve y^3 = x^3 + ax^2 + x is singular for a = ±2")
        PP.<X,Y,Z> = ProjectiveSpace(2, K)
        F = Y^3 - X^3 - a*X^2*Z - X*Z^2
        plane_curve.ProjectivePlaneCurve.__init__(self,PP,F)
    
    # Computes the p-adic integrals from P to Q  of the log differentials 
    #   ω₁ = dx/y^2,  ω₂ = x dx/y^2,  ω₃ = dx/y
    def coleman_integrals_on_basis(self,P,Q):
        K = self.base_ring()
        E = EllipticCurve(K, [0,0,0,0,a^2/4-1])
        to_E = lambda R: E(R[1], R[2] + a/2*R[0], R[0])
        int_omega1 = -3*E.coleman_integrals_on_basis(to_E(P),to_E(Q))[0]

        # integrals of symmetric parts of ω₂ and ω₃ are computed on ℙ¹
        ζ3 = K.primitive_root_of_unity(3)
        w = lambda R: R[0]/R[1] if R[1] != 0 else R[1]^2/(R[0]^2 + a*R[0] + 1)
        int_omega2_symm = -1/2 * sum(ζ^2 * ((w(Q)-ζ)/(w(P)-ζ)).log(0) for ζ in [1,ζ3,ζ3^2])
        int_omega3_symm = -1/2 * sum(ζ * ((w(Q)-ζ)/(w(P)-ζ)).log(0) for ζ in [1,ζ3,ζ3^2])
        
        # integrals of antisymmetric parts are computed on auxiliary curves Xζ
        # where integrands become s ds/(2t)
        int_omega_2_antisymm = 0
        int_omega_3_antisymm = 0
        _.<x> = PolynomialRing(K)
        for ζ in [1,ζ3,ζ3^2]:
            T = x^4 + 12/a^2*ζ^2*x^3 + 12/a^3*ζ*x^2 + 4/a^2*x
            Xζ = HyperellipticCurve(T, check_squarefree = false)
            to_Xζ = lambda R: Xζ(w(R)/(1-w(R)*ζ), -2*(1/R[0]+a/2)/a * w(R)^2/(1-w(R)*ζ)^2) if R[0] != 0 else Xζ(0,0)
            Xζ_int = Xζ.coleman_integrals_on_basis(to_Xζ(P),to_Xζ(Q))[1]
            int_omega_2_antisymm += -ζ * Xζ_int
            int_omega_3_antisymm += -ζ^2 * Xζ_int
        
        int_omega2 = int_omega2_symm + int_omega_2_antisymm
        int_omega3 = int_omega3_symm + int_omega_3_antisymm
        return [int_omega1, int_omega2, int_omega3]



r"""
Determine a Chabauty function for the ℤ[1/487]-integral points on the curve 
  y^3 = x^3 + x^2 + x
as described in Example 5.3 of [LL26].

REFERENCES:
- [LL26] Marius Leonhardt, Martin Lüdtke, "Affine Chabauty II"

AUTHORS:
- Marius Leonhardt
- Martin Lüdtke

"""

load("SupEllInt.sage")

# We consider the curve y^3 = x^3 + ax^2 + x with parameter a = 1
a = 1
p = 7 # auxiliary prime, must be 1 mod 3
N = 10 # precision
P0 = (0,0)  # base point
A = (1/18, 7/18)  # G := A - P0 generates the Mordell-Weil group up to finite index

print(f"Curve: y^3 = x^3 + {a if a != 1 else ''}x^2 + x")
print(f"base point: P₀ = {P0}")
print(f"Mordell—Weil generator: A - P₀ with A = {A}")
print(f"auxiliary prime: {p}")
print(f"precision: {N}")

K = Qp(p,N)
X = SuperEllipticExampleCurve(K, a)

MW_ints = X.coleman_integrals_on_basis(X(P0), X(A))

# reduction type Q₂ modulo (2ζ₃ + 23), i.e. (u,v) = (232,0) mod 487, where u = y/x, v = 1/x
print("S = {487}")
print(f"reduction type: (u,v) = (232,0) mod 487")
print("")

β2 = -log(K(2)) - 1/2*log(K(3))
β3 = β2
ζ3 = K.primitive_root_of_unity(3)
δ2 = sum(-ζ*log(2*ζ+23) for ζ in [ζ3,ζ3^2])
δ3 = sum(-ζ^2*log(2*ζ+23) for ζ in [ζ3,ζ3^2])
# The matrix M(Σ^csp)
M = matrix(K, [
    [MW_ints[0],   MW_ints[1]-β2,   MW_ints[2]-β3],
    [0,            δ2,              δ3]
])
#print(f"matrix M(Σ^csp):")
#print(M)


ker_basis = M.right_kernel().basis()
a1,a2,a3 = ker_basis[0]

# The differential ω = a₁ω₁ + a₂ω₂ + a₃ω₃ vanishes on Z[1/487]-points.
print("Annihilating log differential ω = a₁ω₁ + a₂ω₂ + a₃ω₃ has coefficients")
print(f"  a₁ = {a1}")
print(f"  a₂ = {a2}")
print(f"  a₃ = {a3}")
# We check this for (216/487, 438/487)
P = (216/487, 438/487)
ints_P = X.coleman_integrals_on_basis(X(P0), X(P))
val_P = a1*ints_P[0] + a2*ints_P[1] + a3*ints_P[2]
print(f"Value of Chabauty function ∫_P₀^P ω on P = {P}:")
print(f"  {val_P}")

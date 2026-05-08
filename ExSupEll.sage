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
load("Zproots.sage")  # from https://github.com/martinluedtke/RefinedCK

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

known_points = [(0,0), (216/487, 438/487)]

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

print("Basis of log differentials: ω₁ = dx/y^2,  ω₂ = x dx/y^2,  ω₃ = dx/y")
print("")

ker_basis = M.right_kernel().basis()
a1,a2,a3 = ker_basis[0]

# The differential ω = a₁ω₁ + a₂ω₂ + a₃ω₃ vanishes on Z[1/487]-points.
print("Annihilating log differential ω = a₁ω₁ + a₂ω₂ + a₃ω₃ has coefficients")
print(f"  a₁ = {a1}")
print(f"  a₂ = {a2}")
print(f"  a₃ = {a3}")
print("")

# We check this for (216/487, 438/487)
P = (216/487, 438/487)
ints_P = X.coleman_integrals_on_basis(X(P0), X(P))
val_P = a1*ints_P[0] + a2*ints_P[1] + a3*ints_P[2]
print(f"Value of Chabauty function ∫_P₀^P ω on P = {P}:")
print(f"  {val_P}")
print("")


# find Fp-points and known S-integral points that reduce to them
Xp = X.change_ring(GF(p))
Fppoints = [P for P in Xp.rational_points() if P[2] != 0]
S_integral_lifts = { P: [] for P in Fppoints }
for P in known_points:
    redP = Xp(P)
    S_integral_lifts[redP].append(P)

# helper function to compute floor(log_b(n)) for integer b > 1 and n ≥ 1
def floor_log(b,n):
    e = 0
    while b^(e+1) <= n:
        e += 1
    return e

Chab_locus_known = []
Chab_locus_extra = []

R.<t> = PowerSeriesRing(K,'t',default_prec=N)
x = polygen(QQ)
f = x^3 + a*x^2 + x
f_prime = f.derivative()

print("Computing Chabauty locus...")
for Fppoint in Fppoints:
    if Fppoint[0] * a == -1:
        # we cannot compute Coleman integrals with endpoints in this disc
        print(f"skipping bad residue disc {(Fppoint[0],Fppoint[1])} mod {p}")
        continue
    elif Fppoint[1] != 0:  # non-ramification disc, use parameter t = (x-x0)/p
        print(f"residue disc {(Fppoint[0],Fppoint[1])} mod {p}:")
        xt = Fppoint[0].lift() + p*t
        yt = f(xt)^(1/3)
        for c in [1,ζ3,ζ3^2]:
            if mod(c*yt(0),p) == Fppoint[1]:
                yt *= c
                break
        omega1_over_dt = p / yt^2
    else:  # ramification disc, use parameter t = y/p
        print(f"residue disc {(Fppoint[0],Fppoint[1])} mod {p}:")
        yt = p*t
        # find x(t) with y(t)^3 = f(x(t)) by Newton iteration
        for r in Zproots(f.change_ring(K)):
            if mod(r,p) == Fppoint[0]:
                xt = R(r)
                break
        for i in range(1 + floor_log(2,N)):
            xt -= (f(xt) - yt^3)/f_prime(xt)
        omega1_over_dt = 3*p/f_prime(xt)
        
    # compute power series h(t) representing the Chabauty function on this disc
    P = X(xt(0), yt(0))
    ints_to_P = X.coleman_integrals_on_basis(X(P0), P)
    h = a1*(ints_to_P[0] + omega1_over_dt.integral()) \
      + a2*(ints_to_P[1] + (xt*omega1_over_dt).integral()) \
      + a3*(ints_to_P[2] + (yt*omega1_over_dt).integral())

    # tail_prec = lower bound for p-adic valuation of coefficients of the O(t^N) term of h(t)
    k0 = h.prec()
    tail_prec = k0 - floor_log(p,k0) + min(c.valuation() for c in [a1,a2,a3])
    int_poly = h.truncate()

    # use Zproots function to find all zeros in the residue disc
    for t_root in Zproots(int_poly,tail_prec):
        x_coord = xt.polynomial()(t_root).add_bigoh(tail_prec)
        y_coord = yt.polynomial()(t_root).add_bigoh(tail_prec)
        known_point = False
        for Q in known_points:
            if Q[0] == x_coord and Q[1] == y_coord:
                known_point = True
                Chab_locus_known.append(Q)
                print(f"  {Q}")
                break
        if not known_point:
            Chab_locus_extra.append(X((x_coord, y_coord)))
            print(f"  {(x_coord, y_coord)}")
    
print("")
print(f"Chabauty locus contains {len(Chab_locus_known)} known points and {len(Chab_locus_extra)} extra point{'s' if len(Chab_locus_extra) != 1 else ''}.")

r"""
Determine the integral points on the genus 2 curve y^2 = x^6 + 2*x^5 - 7*x^4 - 18*x^3 + 2*x^2 + 20*x + 9 
  https://www.lmfdb.org/Genus2Curve/Q/6081/b/164187/1
via the Affine Chabauty method. This proves Theorem 5.2 of [LL26].

REFERENCES:
- [LL26] Marius Leonhardt, Martin Lüdtke, "Affine Chabauty II"

AUTHORS:
- Marius Leonhardt
- Martin Lüdtke

"""


x = polygen(ZZ)
f = x^6 + 2*x^5 - 7*x^4 - 18*x^3 + 2*x^2 + 20*x + 9
p = 7
prec = 15

# known solutions (up to hyperelliptic involution)
known_points = [(-1,1),(0,3),(1,3),(-2,3),(-4,37), (-1,-1),(0,-3),(1,-3),(-2,-3),(-4,-37)]
P0 = known_points[0]
MW_generators = [known_points[1], known_points[2]]

K = Qp(p,prec)
X = HyperellipticCurve(f.change_ring(K))
g = X.genus()

load("Zproots.sage")  # from https://github.com/martinluedtke/RefinedCK

print(f"Curve: y^2 = {f}")
print(f"base point: P₀ = {P0}")
print(f"Mordell—Weil generators: " + ", ".join(f"{A}-P₀" for A in MW_generators))
print(f"auxiliary prime: {p}")
print(f"precision: {prec}")


print("Basis of log differentials: ω₀ = dx/y,  ω₁ = x dx/y,  ω₂ = x² dx/y")
print("")

M = matrix(X.coleman_integrals_on_basis(X(P0),X(Q))[:(g+1)] for Q in MW_generators)
ker_basis = M.right_kernel().basis()
coeffs = ker_basis[0]

print("Annihilating log differential ω = a₀ω₀ + a₁ω₁ + a₂ω₂ has coefficients")
for i in range(g+1):
    print(f"    a{i} = {coeffs[i]}")
print("")

def chabauty_function(P):
    ints = X.coleman_integrals_on_basis(X(P0),X(P))
    return sum(coeffs[i]*ints[i] for i in range(g+1))

print("Check that the function ρ(P) = ∫_P₀^P ω vanishes on all known points:")
for P in known_points:
    print(f"    ρ({P}) = {chabauty_function(P)}")
print("")

print(f"The function vanishes on integral points but not necessarily on rational points:")
print(f"    ρ((-1/2, 9/8)) = {chabauty_function(X(-1/2, 9/8))}")
print(f"    ρ((-13/6, 743/216)) = {chabauty_function(X(-13/6, 743/216))}")
print("")

# find Fp-points and known integral points that reduce to them
Xp = HyperellipticCurve(f.change_ring(GF(p)))
Fppoints = [P for P in Xp.points() if P[2] != 0]
integral_lifts = { P: [] for P in Fppoints }
for P in known_points:
    redP = Xp(P)
    integral_lifts[redP].append(P)

# helper function to compute floor(log_b(n)) for integer b > 1 and n ≥ 1
def floor_log(b,n):
    e = 0
    while b^(e+1) <= n:
        e += 1
    return e

Chab_locus_known = []
Chab_locus_extra = []
weierstrass_points = [ W for W in X.weierstrass_points() if W[2] == 1 ]
R.<t> = PowerSeriesRing(K,'t')

print("Computing Chabauty locus...")
for Fppoint in Fppoints:
    if Fppoint[1] != 0:  # non-Weierstrass disc
        print(f"residue disc {(Fppoint[0],Fppoint[1])} mod {p}:")
        xt = Fppoint[0].lift() + p*t
        yt = sqrt(f(xt), prec=prec)
        if Mod(yt(0),p) != Fppoint[1]:
            yt = -yt
        P = X(xt(0), yt(0))
        omega0_over_dt = p / yt
    else:  # Weierstrass disc
        print(f"Weierstrass disc {(Fppoint[0],Fppoint[1])} mod {p}:")
        for W in weierstrass_points:
            if Mod(W[0],p) == Fppoint[0]:
                P = W
                break
        xt,yt = X.local_coordinates_at_weierstrass(P, prec=prec)
        xt = xt(p*t)
        yt = yt(p*t)
        omega0_over_dt = 2*p / f.derivative()(xt)
    
    # compute power series h(t) representing the Chabauty function on this disc
    ints_to_P = X.coleman_integrals_on_basis(X(P0), P)
    h = sum(coeffs[j] * (2*ints_to_P[j] + (xt^j * omega0_over_dt).integral()) for j in range(g+1))

    # tail_prec = lower bound for p-adic valuation of coefficients of the O(p^N) term of h(t)
    k0 = h.prec()
    tail_prec = k0 - floor_log(p,k0) + min(c.valuation() for c in coeffs)
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
            print(f"  {(x_coord, y_coord.add_bigoh(1))}")
    
print("")
print(f"Chabauty locus contains {len(Chab_locus_known)} known points and {len(Chab_locus_extra)} extra point{'s' if len(Chab_locus_extra) != 1 else ''}.")

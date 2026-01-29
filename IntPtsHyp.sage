load("Zproots.sage")  # from https://github.com/martinluedtke/RefinedCK

p = 3  # auxiliary prime
prec = 20  # precision
search_bound = 100  # initial search for integral points with bounded absolute value of x
algdep_bound = 2  # when checking for algebraicity of extra points, use this degree bound
print_prec = 10  # print p-adic numbers with this precision

R.<x> = PolynomialRing(QQ)

# genus 2
# f = R([1, 2, 9, 6, 10, 4, 1])  # 10839.a.32517.1
# f = R([4, -28, 45, -12, -6, 0, 1]) # 193216.a.772864.1
# f = R([1, 10, 19, 4, -1, 2, 1]) #  7884.b.283824.1
# f = R([9, 24, 8, -10, 0, 4, 1]) # 16362.a.883548.1
# f = R([16, 24, 4, -12, -8, 0, 1]) # 6400.f.64000.1
# f = R([9, 20, 2, -18, -7, 2, 1]) # 6081.b.164187.1
# f = R([0, -8, -3, 10, 7, 2, 1]) # 28542.c.171252.1
# f = R([1, 12, 12, 2, 4, 4, 1]) # 14394.a.259092.1
# f = R([4, -12, 1, 10, 3, 2, 1]) # 12330.a.73980.1
# f = R([1, -2, 5, 2, -6, 0, 1]) # 3571.a.3571.1
# f = R([4, -4, -3, 0, -2, 4, 1]) # 38256.b.459072.1
# f = R([1, -20, 24, 2, -8, 0, 1]) # 40878.a.735804.1
# f = R([0, -4, 13, 4, 2, 0, 1]) # 20072.b.642304.1
# f = R([4, 12, 5, -8, 2, 4, 1]) # 147968.a.591872.1
# f = R([1, -12, 66, -54, 5, 2, 1]) # 118284.a.709704.1
# f = R([9, 0, -14, -6, 1, 2, 1]) # 113268.a.679608.1
# f = R([4, 8, 1, -12, -2, 4, 1]) # 26336.a.210688.1
# f = x^6 - 28*x^2 + 4
f = x^6 + 2*x^5 + 5*x^4 + 2*x^3 - 2*x^2 - 4*x - 3 # Atkin-Lehner quotient X_0^+(107), rational points determined in [BDM+23]

# genus 3
# f = x^8 + 6*x^7 + 7*x^6 - 16*x^5 - 31*x^4 + 4*x^3 + 23*x^2 + 6*x + 1 # conductor 2285467
# f = x^8 + 2*x^7 - 3*x^6 - 4*x^5 + 10*x^4 + 6*x^3 - 8*x^2 - 4*x + 1 # conductor 3132257

# genus 4
# f = x^10 + 8*x^9 + 9*x^8 - 11*x^7 - 9*x^6 - 19*x^5 + 9*x^4 + 10*x^3 - 25*x^2 - 12*x + 16  # rankbounds 4 6
# f = x^10 + 6*x^9 - 11*x^7 - 6*x^6 - 14*x^5 + 24*x^4 + x^3 - 18*x^2 + 12*x + 9  # rankbounds 4 5
# f = x^10 + 5*x^9 - 25*x^8 + 6*x^7 + 15*x^6 + 8*x^5 + x^4 + 14*x^3 + 22*x^2 - 21*x + 23 # rankbounds 3 4


if not is_even(f.degree()):
    raise ValueError("Degree of polynomial must be even!")

X = HyperellipticCurve(f)
g = X.genus()

print(f"Hyperelliptic curve of genus {g}: ")
print(f"  y^2 = {f}")

# find some integral points
known_points = []
for a in range(-search_bound+1, search_bound):
    if f(a) == 0:
        known_points.append((a,0))
    elif f(a).is_square():
        b = sqrt(f(a))
        known_points += [(a,b), (a,-b)]

print(f"Found {len(known_points)} integral points")
print(f"  {known_points}")

# find Fp-points and known integral points that reduce to them
Xp = X.change_ring(GF(p))
Fppoints = [P for P in Xp.points() if P[2] != 0]
integral_lifts = { P: [] for P in Fppoints }
for P in known_points:
    redP = Xp(P)
    integral_lifts[redP].append(P)

# basepoint
P0 = known_points[0]

# Construct Mordell—Weil basis consisting of Gᵢ = Pᵢ - P₀ for known integral points Pᵢ.
# These Gᵢ will also form a basis of ker(σ) ⊆ J_Y(Q), provided all Pᵢ have the same reduction type.
# A divisor ∑ (Qᵢ - Rᵢ) is represented as the list of pairs (Qᵢ,Rᵢ).
# (This is a bit more general than what we need here.)
K = Qp(p, prec)
XK = X.change_ring(K)
MW_basis = []
MW_basis_points = []
MW_holom_integrals = []
for P in known_points:
    ints = XK.coleman_integrals_on_basis(XK(P0),XK(P))[:g]
    if vector(ints) not in matrix(MW_holom_integrals, ncols=g).row_space():
        # if the row vector of integrals ∫_P₀^P ωⱼ (j=0,...,g) is independent from
        # the existing rows, add P-P₀ to the MW basis
        MW_basis.append([(XK(P),XK(P0))])
        MW_basis_points.append(P)
        MW_holom_integrals.append(ints)
        if len(MW_basis) == g:
            break

if len(MW_basis) < g:
    raise ValueError(f"Did not find enough independent elements of the Mordell—Weil group! Found {len(MW_basis)} but need {g}.")

print("Mordell—Weil basis: ")
for P in MW_basis_points:
    print(f"  {P} - {P0}")
print("")

# compute coleman integrals ∫_Gᵢ ωⱼ
MW_basis_integrals = []
for i in range(len(MW_basis)):
    integrals_Gi = [0] * (g+1)
    for Q,R in MW_basis[i]:
        ints_QR = XK.coleman_integrals_on_basis(R,Q)
        for j in range(g+1):
            integrals_Gi[j] += 2*ints_QR[j]
    MW_basis_integrals.append(integrals_Gi)

# compute coefficients of the Chabauty differential ω = ∑ cⱼωⱼ
# by solving the system of equations ∑ cⱼ ∫_Gᵢ ωⱼ = 0
ker = matrix(MW_basis_integrals).right_kernel()
ker_basis = ker.basis()
assert len(ker_basis) == 1
coeffs = ker_basis[0]

print("Found Chabauty differential ω = ∑ cⱼωⱼ with coefficients")
for j in range(g+1):
    print(f"  c_{j} = {coeffs[j].add_bigoh(print_prec)}")
print("")

def chabauty_function(P):
    ints = XK.coleman_integrals_on_basis(XK(P0),XK(P))
    return sum(coeffs[j] * ints[j] for j in range(g+1))

print("Check that Chabauty function vanishes on known integral points:")
vanishes_on_all = True
for P in known_points:
    val = chabauty_function(P)
    print(f"  {P}: {val}")
    if val != 0:
        vanishes_on_all = False

if not vanishes_on_all:
    print("Chabauty function does not vanish on all known integral points!")
    print("Maybe there is more than one reduction type, or the Mordell—Weil rank is larger than the genus.")
print("")

# helper function to compute floor(log_b(n)) for integer b > 1 and n ≥ 1
def floor_log(b,n):
    e = 0
    while b^(e+1) <= n:
        e += 1
    return e

weierstrass_points = [ W for W in XK.weierstrass_points() if W[2] == 1 ]

Chab_locus_known = []
Chab_locus_extra = []

_.<t> = PowerSeriesRing(K,'t')

print("Computing Chabauty locus...")
for Fppoint in Fppoints:
    if Fppoint[1] != 0:  # non-Weierstrass disc
        print(f"residue disc {(Fppoint[0],Fppoint[1])} mod {p}:")
        xt = Fppoint[0].lift() + p*t
        yt = sqrt(f(xt), prec=prec)
        if Mod(yt(0),p) != Fppoint[1]:
            yt = -yt
        P = XK(xt(0), yt(0))
        omega0_over_dt = p / yt
    else:  # Weierstrass disc
        print(f"Weierstrass disc {(Fppoint[0],Fppoint[1])} mod {p}:")
        for W in weierstrass_points:
            if Mod(W[0],p) == Fppoint[0]:
                P = W
                break
        xt,yt = XK.local_coordinates_at_weierstrass(P, prec=prec)
        xt = xt(p*t)
        yt = yt(p*t)
        omega0_over_dt = 2*p / f.derivative()(xt)
    
    # compute power series h(t) representing the Chabauty function on this disc
    ints_to_P = XK.coleman_integrals_on_basis(XK(P0), P)
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
            Chab_locus_extra.append(XK((x_coord, y_coord)))
            print(f"  {(x_coord.add_bigoh(print_prec), y_coord.add_bigoh(1))}")
    
print("")
print(f"Chabauty locus contains {len(Chab_locus_known)} known points and {len(Chab_locus_extra)} extra point{'s' if len(Chab_locus_extra) != 1 else ''}.")

if Chab_locus_extra:
    print("")
    print(f"Check for algebraicity of the extra point{'s' if len(Chab_locus_extra) != 1 else ''}...")
    for P in Chab_locus_extra:
        if P[1] == 0:
            print(f"  {(P[0].add_bigoh(print_prec), P[1].add_bigoh(print_prec))}: Weierstrass point")
        else:
            print(f"  {(P[0].add_bigoh(print_prec), P[1].add_bigoh(1))}: {P[0].algdep(algdep_bound)}")
# AffChab2
Computing S-integral points via the Affine Chabauty method: examples of hyperelliptic and superelliptic curves

Sage code for the paper "Affine Chabauty II" [[LL26]](https://arxiv.org/abs/2602.05643)

The file `SupEllInt.sage` defines a class for curves of the form $y^3 = x^3 + ax^2 + x$ and a function `coleman_integrals_on_basis` to compute Coleman integrals of the logarithmic differentials $dx/y^2$, $x dx/y^2$, $dx/y$.
The file `ExSupEll.sage` contains the Affine Chabauty computations for the $\mathbb{Z}[1/487]$-points on the curve with parameter $a = 1$. It determines an annihilating log differential and verifies its vanishing on the point $(216/487, 438/487)$, providing the computations for Example 5.3 of [LL26]. The file `ExHyp.sage` computes the Chabauty locus for the integral points on the hyperelliptic curve 

$$ y^2 = x^6 + 2x^5 - 7x^4 - 18x^3 + 2x^2 + 20x + 9, $$

showing [LL26], Theorem 5.2.

To run the code, place the Sage files in the working directory and call `sage ExSupEll.sage` or `sage ExHyp.sage`. One also needs the code from [https://github.com/jbalakrishnan/AWS](https://github.com/jbalakrishnan/AWS) to compute Coleman integrals and the file `Zproots.sage` from [https://github.com/martinluedtke/RefinedCK](https://github.com/martinluedtke/RefinedCK) to compute roots of p-adic polynomials. The code was tested on Sage 10.8. 

The output of `ExSupEll.sage` should look as follows:
```
Curve: y^3 = x^3 + x^2 + x
base point: P₀ = (0, 0)
Mordell—Weil generator: A - P₀ with A = (1/18, 7/18)
auxiliary prime: 7
precision: 20
S = {487}
reduction type: (u,v) = (232,0) mod 487

Basis of log differentials: ω₁ = dx/y^2,  ω₂ = x dx/y^2,  ω₃ = dx/y

Annihilating log differential ω = a₁ω₁ + a₂ω₂ + a₃ω₃ has coefficients
  a₁ = 1 + O(7^18)
  a₂ = 2 + 6*7 + 2*7^2 + 3*7^3 + 4*7^5 + 5*7^6 + 2*7^7 + 4*7^8 + 5*7^9 + 3*7^10 + 5*7^12 + 4*7^13 + 6*7^14 + 3*7^15 + 6*7^16 + O(7^18)
  a₃ = 2*7 + 6*7^2 + 2*7^5 + 4*7^7 + 6*7^8 + 3*7^10 + 7^11 + 3*7^12 + 4*7^13 + 6*7^14 + 5*7^15 + 2*7^16 + 5*7^17 + 7^18 + O(7^19)

Value of Chabauty function ∫_P₀^P ω on P = (216/487, 438/487):
  O(7^19)
```

The output of `ExHyp.sage` should look as follows:
```
Curve: y^2 = x^6 + 2*x^5 - 7*x^4 - 18*x^3 + 2*x^2 + 20*x + 9
base point: P₀ = (-1, 1)
Mordell—Weil generators: (0, 3)-P₀, (1, 3)-P₀
auxiliary prime: 7
precision: 15
Basis of log differentials: ω₀ = dx/y,  ω₁ = x dx/y,  ω₂ = x² dx/y

Annihilating log differential ω = a₀ω₀ + a₁ω₁ + a₂ω₂ has coefficients
    a0 = 1 + O(7^14)
    a1 = 5 + 3*7 + 3*7^2 + 5*7^3 + 3*7^4 + 2*7^6 + 2*7^7 + 7^8 + 4*7^9 + 3*7^10 + 5*7^11 + 5*7^12 + 6*7^13 + O(7^14)
    a2 = 5 + 6*7 + 6*7^2 + 7^3 + 4*7^4 + 6*7^5 + 5*7^6 + 3*7^7 + 4*7^8 + 2*7^9 + 6*7^10 + 3*7^12 + 4*7^13 + O(7^14)

Check that the function ρ(P) = ∫_P₀^P ω vanishes on all known points:
    ρ((-1, 1)) = 0
    ρ((0, 3)) = O(7^15)
    ρ((1, 3)) = O(7^15)
    ρ((-2, 3)) = O(7^15)
    ρ((-4, 37)) = O(7^15)
    ρ((-1, -1)) = O(7^15)
    ρ((0, -3)) = O(7^15)
    ρ((1, -3)) = O(7^15)
    ρ((-2, -3)) = O(7^15)
    ρ((-4, -37)) = O(7^15)

The function vanishes on integral points but not necessarily on rational points:
    ρ((-1/2, 9/8)) = 5*7 + 5*7^2 + 2*7^3 + 6*7^4 + 4*7^5 + 6*7^6 + 3*7^7 + 7^8 + 2*7^9 + 7^10 + 4*7^11 + 5*7^12 + 5*7^13 + 4*7^14 + O(7^15)
    ρ((-13/6, 743/216)) = 3*7 + 5*7^2 + 6*7^3 + 5*7^4 + 2*7^6 + 5*7^9 + 7^11 + 6*7^12 + 5*7^14 + O(7^15)

Computing Chabauty locus...
residue disc (0, 4) mod 7:
  (0, -3)
residue disc (0, 3) mod 7:
  (0, 3)
residue disc (1, 4) mod 7:
  (1, -3)
residue disc (1, 3) mod 7:
  (1, 3)
residue disc (3, 5) mod 7:
  (-4, -37)
residue disc (3, 2) mod 7:
  (-4, 37)
residue disc (5, 4) mod 7:
  (-2, -3)
residue disc (5, 3) mod 7:
  (-2, 3)
residue disc (6, 6) mod 7:
  (-1, -1)
residue disc (6, 1) mod 7:
  (-1, 1)

Chabauty locus contains 10 known points and 0 extra points.
```

## Authors

- Marius Leonhardt
- Martin Lüdtke

## References
- [LL26] M. Leonhardt, Martin Lüdtke, "Affine Chabauty II" [(arXiv)](https://arxiv.org/abs/2602.05643)

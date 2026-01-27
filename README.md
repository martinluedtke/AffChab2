# AffChab2
Computing S-integral points via the Affine Chabauty method: examples of hyperelliptic and superelliptic curves

Sage code for the paper "Affine Chabauty II" by M. Leonhardt and M. Lüdtke.

The file `SupEllInt.sage` defines a class for curves of the form $y^3 = x^3 + ax^2 + x$ and a function `coleman_integrals_on_basis` to compute Coleman integrals of the logarithmic differentials $dx/y^2$, $x dx/y^2$, $dx/y$.
The file `ExSupEll.sage` contains the Affine Chabauty computations for the $\mathbb{Z}[1/487]$-points on the curve with parameter $a = 1$. It determines an annihilating log differential and verifies its vanishing on the point $(216/487, 438/487)$.

To run the code, place the Sage files in the working directory and call `sage ExSupEll.sage`. One also needs the code from [https://github.com/jbalakrishnan/AWS](https://github.com/jbalakrishnan/AWS) to compute Coleman integrals. The code was tested on Sage 10.8. The output should look as follows:
```
Curve: y^3 = x^3 + x^2 + x
base point: P₀ = (0, 0)
Mordell—Weil generator: A - P₀ with A = (1/18, 7/18)
auxiliary prime: 7
precision: 20
S = 487
reduction type: (u,v) = (232,0) mod 487

Annihilating log differential ω = a₁ω₁ + a₂ω₂ + a₃ω₃ has coefficients
  a₁ = 1 + O(7^18)
  a₂ = 2 + 6*7 + 2*7^2 + 3*7^3 + 4*7^5 + 5*7^6 + 2*7^7 + 4*7^8 + 5*7^9 + 3*7^10 + 5*7^12 + 4*7^13 + 6*7^14 + 3*7^15 + 6*7^16 + O(7^18)
  a₃ = 2*7 + 6*7^2 + 2*7^5 + 4*7^7 + 6*7^8 + 3*7^10 + 7^11 + 3*7^12 + 4*7^13 + 6*7^14 + 5*7^15 + 2*7^16 + 5*7^17 + 7^18 + O(7^19)
Value of Chabauty function ∫_P₀^P ω on P = (216/487, 438/487):
  O(7^19)
```

## Authors

- Marius Leonhardt
- Martin Lüdtke

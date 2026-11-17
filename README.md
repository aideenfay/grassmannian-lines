# Lines on Grassmannian Zero Loci

This repository provides a simple Sagemath library to count the number of lines through a general point of a variety
\[
Y=Z(F)\subset\operatorname{Gr}(k,n),
\]
where Gr(k,n) is the Grassmannian parameterizing \(k\)-dimensional subspaces of an \(n\)-dimensional vector space, and F is the defining bundle, assumed to be a direct sum of globally generated irreducible homogeneous vector bundles. All calculations are over \(\mathbb C\). We assume a general section with nonempty regular zero locus and count lines of Plücker degree one.

We use [Atiyah–Bott’s equivariant integration formula (3.8)](https://www.math.stonybrook.edu/~mmovshev/MAT570Spring2008/BOOKS/atiyahbott_moment.pdf#page=9) to evaluate the Chern class integrals underlying the calculation. We apply it on \(\mathbf P^{k-1}\times\mathbf P^{n-k-1}\), which parametrizes lines through a fixed point of \(\operatorname{Gr}(k,n)\).  

## Why count lines?
In progress! 

## Getting started
To use it, you must have Sagemath installed (see [here](https://doc.sagemath.org/html/en/installation/index.html) ). 

From the repository root, open Sage or a Sage notebook and run:

```python
load("lines.sage")
```

To begin with an example, run:
```bash
sage examples/kuchle.sage
```
NOTE: the documentation within `lines.sage` and `bundle.py` is more detailed (and perhaps instructive) than what will be written here.

## Usage

The class `Bundle`, defined in `bundle.py`, records the Grassmannian and the torus weights of the defining bundle.

| Function | Output |
|---|---|
| `lines_through_point(F: Bundle) -> Integer` | The number of lines when the expected family dimension is zero. |
| `line_degree(F: Bundle) -> Integer` | The projective degree of the family of tangent directions, including when it has positive dimension. |

We use:

| Notation | Bundle |
|---|---|
| `S` | Tautological subbundle |
| `Q` | Tautological quotient bundle |
| `S*`, `Q*` | Their duals |
| `O(1)` | Plücker line bundle \(\det S^\vee\) |
| `Sym2 S*` | \(\operatorname{Sym}^2 S^\vee\) |
| `Wedge2 S*` | \(\bigwedge^2 S^\vee\) |
| `Q*(1)` | \(Q^\vee\otimes\mathcal O(1)\) |


To construct a bundle, supply the Grassmannian dimensions and a direct-sum expression:

```python
F = Bundle(2, 7, "Sym2 S* + O(1) + O(1)")
```

or a list of summands:

```python
F = Bundle(2, 7, ["Sym2 S*", "O(1)", "O(1)"])
```

The supported names appear in `Bundle._parse` in `bundle.py`. You can use the bundle methods below for other constructions.

## Examples

### A cubic threefold

A cubic threefold is the zero locus of \(\mathcal O(3)\) on
\(\mathbf P^4=\operatorname{Gr}(1,5)\). There are six lines through a general point:

```python
lines_through_point(1, 5, "O(3)")
# 6
```

### A zero locus in \(\operatorname{Gr}(2,7)\)

For the fivefold defined by
\[
F=\operatorname{Sym}^2S^\vee\oplus\mathcal O(1)^{\oplus2},
\]
the count is four:

```python
lines_through_point(2, 7, "Sym2 S* + O(1) + O(1)")
# 4
```

### A positive-dimensional family

A positive-dimensional family means that the lines through a fixed point vary in a parameter space of positive dimension, so there are infinitely many of them. When the expected family dimension is positive, use `line_degree`.

The lines through a point of \(\operatorname{Gr}(3,6)\) form
\(\mathbf P^2\times\mathbf P^2\). Their tangent directions give its Segre embedding, of degree six:

```python
line_degree(3, 6, "")
# 6
```

The empty string denotes the zero bundle (so the zero locus is the entire Grassmannian).

### Constructing bundles

For more general expressions, use the bundle operations:

```python
S, Q, O1 = Bundle.tautological(3, 7)

F = S.dual().wedge(2) + Q.wedge(3)

lines_through_point(F)
# 4
```

Here `S`, `Q` and `O1` are `Bundle` objects. The operations are:

| Operation | Construction |
|---|---|
| `E + F` | Direct sum |
| `E.dual()` | Dual |
| `E * F` | Tensor product |
| `E.sym(m)` | Symmetric power |
| `E.wedge(m)` | Exterior power |

### Küchle fourfolds

| File | Purpose |
|---|---|
| `lines.sage` | General localization calculation |
| `data/kuchle.json` | Grassmannians and defining bundles |
| `examples/kuchle.sage` | Computes and prints the Küchle table |
| `tests/checks.sage` | Checks against known answers |
| `tests/test_bundles.py` | Bundle unit tests |


To run this example, open your terminal and run:

```bash
sage examples/kuchle.sage
```

Alternatively, in Sage:

```python
load("examples/kuchle.sage")
```
```python
kuchle_lines("b8")
# 4

[kuchle_lines(t) for t in ("b8", "b11", "c1", "c6")]
# [4, 3, 6, 5]
```
## Tests

From the repository root, run:

```bash
sage tests/checks.sage
python3 -m unittest discover -s tests -v
```

## References

- Atiyah, M.F., & Bott, R. (1984). The moment map and equivariant cohomology. Topology, 23, 1-28.







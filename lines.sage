"""
Lines through a point in the zero locus of a Grassmannian.

Work over C. Let F be a direct sum of globally generated irreducible
homogeneous bundles on Gr(k,n), whose general zero locus Y is nonempty
and regular. All lines have degree one for the Plucker polarization.

"""

from sage.all import ZZ, QQ, vector, prod
from bundle import Bundle


def _bundle_data(k, n=None, bundle=None):
    if isinstance(k, Bundle):
        if n is not None or bundle is not None:
            raise TypeError("Use f(F) or f(k,n,F).")
        E = k
    else:
        E = Bundle(k, n, bundle)
    k, n = E.k, E.n
    F = [vector(ZZ, w) for w in E.weights]
    for i in range(k):
        for j in range(k, n):
            if any(w[j] < w[i] for w in F):
                raise ValueError("F must have nonnegative splitting degrees on lines")
    return k, n, F, ZZ(E.c1_degree)


def _line_degree(k, n, F, r, t=None):
    """Evaluate localization using the provided bundle."""
    t = vector(QQ, range(n) if t is None else t)
    if len(t) != n or len(set(t)) != n:
        raise ValueError("Choose n pairwise distinct torus parameters")
    if r > n - 2:
        return ZZ(0)

    def bott(i, j):
        a = t[i] - t[j]
        euler = prod(
            w.dot_product(t) + m * a for w in F for m in range(1, w[j] - w[i] + 1)
        )
        tangent = prod(t[i] - t[p] for p in range(k) if p != i) * prod(
            t[q] - t[j] for q in range(k, n) if q != j
        )
        return euler * a ** (n - 2 - r) / tangent

    answer = QQ(sum(bott(i, j) for i in range(k) for j in range(k, n)))
    if answer.denominator() != 1 or answer < 0:
        raise ArithmeticError(
            "Expected a nonnegative integer; check the bundle character"
        )
    return ZZ(answer)


def line_degree(k, n=None, bundle=None, t=None):
    """Compute integral c_r(H) h^(n-2-r) on P^(k-1) x P^(n-k-1).

    Here r=c1(F).line and h is the hyperplane class. Use line_degree(F)
    for a Bundle, or line_degree(k,n,F) for any supported input.
    """
    return _line_degree(*_bundle_data(k, n, bundle), t=t)


def lines_through_point(k, n=None, bundle=None, t=None):
    """Count a zero-dimensional family; return zero if expected dimension is negative.
    """
    k, n, F, r = _bundle_data(k, n, bundle)
    if r < n - 2:
        raise ValueError(
            "The expected line family dimension is positive; use line_degree"
        )
    return _line_degree(k, n, F, r, t=t)

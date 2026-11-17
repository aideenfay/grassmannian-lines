from collections import Counter
from dataclasses import dataclass
from itertools import combinations, combinations_with_replacement
from operator import index
import re


def _integer(value, name):
    try:
        return int(index(value))
    except TypeError as exc:
        raise TypeError(f"{name} must be an integer") from exc


@dataclass(frozen=True, init=False)
class Bundle:
    """A completely reducible homogeneous bundle on Gr(k,n).

    INPUT:
    - k: dimension of the subspaces parametrized by Gr(k,n).
    - n: dimension of the ambient vector space.
    - data: a bundle expression, where empty input denotes the zero bundle.
    """

    k: int
    n: int
    weights: tuple

    def __init__(self, k, n, data=""):
        k, n = _integer(k, "k"), _integer(n, "n")
        if not 1 <= k < n:
            raise ValueError("Require 1 <= k < n")
        if isinstance(data, Bundle):
            if (data.k, data.n) != (k, n):
                raise ValueError("Bundles must lie on the same Grassmannian.")
            raw = data.weights
        elif isinstance(data, str):
            raw = self._parse(k, n, data.split("+") if data.strip() else [])
        else:
            entries = list(data)
            raw = (
                self._parse(k, n, entries)
                if all(isinstance(x, str) for x in entries)
                else entries
            )
        character = tuple(
            sorted(tuple(_integer(x, "weight coordinate") for x in w) for w in raw)
        )
        if any(len(w) != n for w in character):
            raise ValueError("Every weight must have length n")
        multiplicities = Counter(character)
        for p in list(range(k - 1)) + list(range(k, n - 1)):
            swapped = Counter(w[:p] + (w[p + 1], w[p]) + w[p + 2 :] for w in character)
            if swapped != multiplicities:
                raise ValueError("Weights must be invariant under S_k x S_(n-k)")
        object.__setattr__(self, "k", k)
        object.__setattr__(self, "n", n)
        object.__setattr__(self, "weights", character)

    @staticmethod
    def _parse(k, n, summands):
        if not summands:
            return ()
        S, Q, O1 = Bundle.tautological(k, n)
        names = {
            "S": lambda: S,
            "Q": lambda: Q,
            "S*": S.dual,
            "Q*": Q.dual,
            "Sym2 S*": lambda: S.dual().sym(2),
            "Wedge2 S*": lambda: S.dual().wedge(2),
            "Wedge3 Q": lambda: Q.wedge(3),
            "S*(1)": lambda: S.dual() * O1,
            "Q(1)": lambda: Q * O1,
            "Q*(1)": lambda: Q.dual() * O1,
        }
        result = []
        for term in summands:
            term = term.strip()
            line = re.fullmatch(r"O\((-?\d+)\)", term)
            if line:
                E = Bundle.line(k, n, int(line[1]))
            elif term in names:
                E = names[term]()
            else:
                raise ValueError(
                    f"Unknown summand {term}: use Bundle methods for other constructions"
                )
            result.extend(E.weights)
        return result

    @classmethod
    def tautological(cls, k, n):
        """Return S, Q and O(1)=det(S*) as Bundle objects on Gr(k,n)."""
        k, n = _integer(k, "k"), _integer(n, "n")
        if not 1 <= k < n:
            raise ValueError("Require 1 <= k < n")
        basis = [tuple(int(i == j) for j in range(n)) for i in range(n)]
        return cls(k, n, basis[:k]), cls(k, n, basis[k:]), cls.line(k, n, 1)

    @classmethod
    def line(cls, k, n, d=1):
        """The line bundle O(d); where O(0) is the trivial line."""
        k, n, d = _integer(k, "k"), _integer(n, "n"), _integer(d, "d")
        return cls(k, n, [tuple(-d if j < k else 0 for j in range(n))])

    @property
    def rank(self):
        return len(self.weights)

    @property
    def c1_degree(self):
        """The integer c1(F).line; not a test of global generation."""
        return sum(w[-1] - w[0] for w in self.weights)

    def __repr__(self):
        return f"Bundle(Gr({self.k},{self.n}), rank={self.rank}, c1_degree={self.c1_degree})"

    def _same_base(self, other):
        if not isinstance(other, Bundle):
            raise TypeError("Expected a Bundle")
        if (self.k, self.n) != (other.k, other.n):
            raise ValueError("Bundles must lie on the same Grassmannian")

    def __add__(self, other):
        if not isinstance(other, Bundle):
            return NotImplemented
        self._same_base(other)
        return Bundle(self.k, self.n, self.weights + other.weights)

    def __radd__(self, other):
        if other == 0:
            return self
        return NotImplemented

    def copies(self, m):
        m = _integer(m, "number of copies")
        if m < 0:
            raise ValueError("The number of copies must be nonnegative")
        return Bundle(self.k, self.n, self.weights * m)

    def tensor(self, other):
        self._same_base(other)
        return Bundle(
            self.k,
            self.n,
            [
                tuple(x + y for x, y in zip(v, w))
                for v in self.weights
                for w in other.weights
            ],
        )

    def __mul__(self, other):
        return self.tensor(other) if isinstance(other, Bundle) else self.copies(other)

    def __rmul__(self, other):
        return self.copies(other)

    def dual(self):
        return Bundle(self.k, self.n, [tuple(-x for x in w) for w in self.weights])

    def twist(self, d):
        """Tensor with O(d)."""
        return self.tensor(Bundle.line(self.k, self.n, d))

    def _power(self, m, symmetric):
        m = _integer(m, "power")
        if m < 0:
            raise ValueError("The power must be nonnegative")
        if m == 0:
            return Bundle.line(self.k, self.n, 0)
        choose = combinations_with_replacement if symmetric else combinations
        return Bundle(
            self.k,
            self.n,
            [
                tuple(sum(w[j] for w in c) for j in range(self.n))
                for c in choose(self.weights, m)
            ],
        )

    def sym(self, m):
        """The m-th symmetric power, retaining weight multiplicities."""
        return self._power(m, True)

    def wedge(self, m):
        """The m-th exterior power, retaining weight multiplicities."""
        return self._power(m, False)

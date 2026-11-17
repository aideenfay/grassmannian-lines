import json
import sys
from pathlib import Path
from sage.all import binomial, factorial

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
load(str(ROOT / "lines.sage"))

assert lines_through_point(1, 5, "O(3)") == 6
assert lines_through_point(2, 5, "O(1) + O(1) + O(1)") == 3
assert lines_through_point(5, 10, "Wedge2 S* + Wedge2 S*") == 5

# Independent degrees: complete intersections and Segre embeddings.
for d in range(1, 7):
    assert lines_through_point(1, d + 2, "O(%s)" % d) == factorial(d)
for n in range(3, 9):
    for k in range(1, n):
        assert line_degree(k, n, "") == binomial(n - 2, k - 1)

S, Q, O1 = Bundle.tautological(3, 7)
F = S.dual().wedge(2) + Q.wedge(3)
assert lines_through_point(3, 7, F) == 4
assert S.sym(0) == Bundle.line(3, 7, 0)
assert S.wedge(0) == Bundle.line(3, 7, 0)
assert lines_through_point(F) == 4

S, Q, O1 = Bundle.tautological(2, 7)
F = S.dual().sym(2) + 2 * O1
assert F == Bundle(2, 7, "Sym2 S* + O(1) + O(1)")
assert F.rank == 5 and F.c1_degree == 5
assert lines_through_point(F) == 4
assert lines_through_point(2, 7, F) == 4
assert line_degree(Bundle(3, 6)) == 6
assert 0 * O1 != Bundle.line(2, 7, 0)
assert (0 * O1).sym(0) == Bundle.line(2, 7, 0)

try:
    lines_through_point(3, 7, F)
except ValueError:
    pass
else:
    raise AssertionError("A bundle must be valid for the provided Grassmannian.")

with open(ROOT / "data/kuchle.json", encoding="utf-8") as source:
    entries = json.load(source)["varieties"]

counts = {}
for entry in entries:
    k, n = entry["k"], entry["n"]
    summands = list(entry["summands"])
    _, _, F, r = _bundle_data(k, n, summands)
    assert k * (n - k) - len(F) == 4 and n - r == 1
    if "O(1)" in summands:
        summands.remove("O(1)")
        E = Bundle(k, n, summands)
        a = lines_through_point(E)
        assert lines_through_point(k, n, summands) == a
        assert (
            lines_through_point(k, n, summands, t=[i * i + 13 for i in range(n)]) == a
        )
        counts[entry["type"]] = a

assert [counts[t] for t in ("b8", "b11", "c1", "c6")] == [4, 3, 6, 5]
print("SUCCESS: All SageMath checks passed.")

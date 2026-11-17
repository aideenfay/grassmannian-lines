from dataclasses import FrozenInstanceError
from math import comb
from pathlib import Path
import json
import unittest

from bundle import Bundle


class BundleTests(unittest.TestCase):
    def setUp(self):
        self.S, self.Q, self.O1 = Bundle.tautological(2, 7)

    def test_tautological_invariants(self):
        self.assertEqual(
            [(E.rank, E.c1_degree) for E in (self.S, self.Q, self.O1)],
            [(2, -1), (5, 1), (1, 1)],
        )
        self.assertEqual(self.S.dual().wedge(2), self.O1)

    def test_string_list_and_constructed_bundle_agree(self):
        F = self.S.dual().sym(2) + 2 * self.O1
        self.assertEqual(F, Bundle(2, 7, "Sym2 S* + O(1) + O(1)"))
        self.assertEqual(F, Bundle(2, 7, ["Sym2 S*", "O(1)", "O(1)"]))
        self.assertEqual((F.rank, F.c1_degree), (5, 5))

    def test_character_order_does_not_affect_equality(self):
        self.assertEqual(Bundle(2, 7, reversed(self.S.weights)), self.S)
        self.assertEqual(self.S + self.Q, self.Q + self.S)

    def test_tensor_product_and_twist(self):
        self.assertEqual(self.Q.dual().twist(1), Bundle(2, 7, "Q*(1)"))
        self.assertEqual(self.S * self.O1, self.S.tensor(self.O1))
        self.assertEqual((self.S * self.Q).rank, self.S.rank * self.Q.rank)
        self.assertEqual(
            (self.S * self.Q).c1_degree,
            self.Q.rank * self.S.c1_degree + self.S.rank * self.Q.c1_degree,
        )

    def test_zero_and_unit_bundles(self):
        zero, unit = Bundle(2, 7), Bundle.line(2, 7, 0)
        self.assertNotEqual(zero, unit)
        self.assertEqual((zero.rank, unit.rank), (0, 1))
        self.assertEqual(zero + self.S, self.S)
        self.assertEqual(zero * self.S, zero)
        self.assertEqual(unit * self.S, self.S)
        self.assertEqual(zero.sym(0), unit)
        self.assertEqual(zero.wedge(0), unit)
        self.assertEqual(zero.sym(1), zero)

    def test_power_ranks(self):
        for m in range(6):
            with self.subTest(m=m):
                self.assertEqual(self.Q.sym(m).rank, comb(self.Q.rank + m - 1, m))
                self.assertEqual(self.Q.wedge(m).rank, comb(self.Q.rank, m))
        self.assertEqual(self.S.wedge(3).rank, 0)

    def test_multiplicities(self):
        E = 2 * self.O1
        self.assertEqual(E.sym(2), 3 * Bundle.line(2, 7, 2))
        self.assertEqual(E.wedge(2), Bundle.line(2, 7, 2))

    def test_same_base_is_required(self):
        # Same n and rank still does not mean the same Grassmannian.
        other = Bundle.line(3, 7)
        for operation in (
            lambda: self.O1 + other,
            lambda: self.O1 * other,
            lambda: Bundle(2, 7, other),
        ):
            with self.assertRaises(ValueError):
                operation()

    def test_frozen_state(self):
        with self.assertRaises(FrozenInstanceError):
            self.S.k = 3
        with self.assertRaises(TypeError):
            self.S.weights[0][0] = 100
        self.assertEqual(len({self.S, self.S.dual().dual()}), 1)

    def test_sum_and_copies(self):
        self.assertEqual(sum([self.O1, self.O1]), self.O1 * 2)
        self.assertEqual(0 * self.O1, Bundle(2, 7))
        for m in (-1, 1.5):
            with self.assertRaises((TypeError, ValueError)):
                self.O1.copies(m)

    def test_invalid_dimensions_and_weights(self):
        cases = [
            (0, 7, []),
            (7, 7, []),
            (2.0, 7, []),
            (2, 7, [(0,) * 6]),
            (2, 7, [(0.5,) * 7]),
            (2, 7, [(-1, 0, 0, 0, 0, 0, 0)]),
        ]
        for k, n, data in cases:
            with self.subTest(k=k, n=n, data=data), self.assertRaises(
                (TypeError, ValueError)
            ):
                Bundle(k, n, data)

    def test_negative_bundles_can_be_constructed(self):
        self.assertEqual(Bundle(2, 7, "O(-2)").c1_degree, -2)
        self.assertEqual(self.S.dual().dual(), self.S)

    def test_all_data_presentations(self):
        path = Path(__file__).resolve().parents[1] / "data/kuchle.json"
        entries = json.loads(path.read_text())["varieties"]
        self.assertEqual(len(entries), 21)
        for row in entries:
            k, n = row["k"], row["n"]
            with self.subTest(type=row["type"]):
                F = Bundle(k, n, row["summands"])
                self.assertEqual(k * (n - k) - F.rank, 4)
                self.assertEqual(n - F.c1_degree, 1)
                self.assertNotIn("classification_status", row)


if __name__ == "__main__":
    unittest.main()

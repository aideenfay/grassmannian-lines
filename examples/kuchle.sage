import sys
from pathlib import Path

sys.path.insert(0, str(Path.cwd()))
load("lines.sage")
import json

with open("data/kuchle.json", encoding="utf-8") as source:
    KUCHLE = {entry["type"]: entry for entry in json.load(source)["varieties"]}


def kuchle_lines(label):
    """Removes O(1) and counts the numebr of lines through a general point of the fivefold."""
    entry = KUCHLE[label]
    summands = list(entry["summands"])
    if "O(1)" not in summands:
        return None
    summands.remove("O(1)")
    F = Bundle(entry["k"], entry["n"], summands)
    return lines_through_point(F)


def kuchle_table():
    """Print the fourfold bundles and their associated fivefold line counts."""
    headers = ("Type", "Lines", "Ambient", "Defining bundle F")
    rows = []
    for label, entry in KUCHLE.items():
        a = kuchle_lines(label)
        rows.append(
            (
                label,
                "—" if a is None else str(a),
                "Gr(%s, %s)" % (entry["k"], entry["n"]),
                " + ".join(entry["summands"]),
            )
        )

    widths = [max(len(cell) for cell in column) for column in zip(headers, *rows)]

    def rule(left, middle, right):
        print(left + middle.join("─" * (width + 2) for width in widths) + right)

    def row(cells):
        padded = [
            cell.rjust(width) if i == 1 else cell.ljust(width)
            for i, (cell, width) in enumerate(zip(cells, widths))
        ]
        print("│ " + " │ ".join(padded) + " │")

    print("Küchle fourfolds X = Z(F)\n")
    rule("┌", "┬", "┐")
    row(headers)
    rule("├", "┼", "┤")
    previous_family = None
    for cells in rows:
        family = cells[0][0]
        if previous_family is not None and family != previous_family:
            rule("├", "┼", "┤")
        row(cells)
        previous_family = family
    rule("└", "┴", "┘")
    print("\nLines: through a general point of the fivefold obtained")
    print("       by removing one O(1) summand from F.")

kuchle_table()

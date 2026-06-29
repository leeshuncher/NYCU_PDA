#!/usr/bin/env python3

from __future__ import annotations

import math
import sys
from dataclasses import dataclass, field


@dataclass
class Net:
    name: str
    idx: int
    connected_cells: list[int] = field(default_factory=list)


@dataclass
class Cell:
    name: str
    idx: int
    connected_nets: list[int] = field(default_factory=list)
    net_pin_idx: list[int] = field(default_factory=list)


@dataclass
class ParsedInput:
    balance_factor: float = 0.0
    nets: list[Net] = field(default_factory=list)
    cells: list[Cell] = field(default_factory=list)


@dataclass
class SolutionGroup:
    label: str = ""
    declared_size: int = 0
    cells: list[str] = field(default_factory=list)


@dataclass
class ParsedSolution:
    declared_cutsize: int = 0
    groups: list[SolutionGroup] = field(default_factory=list)


def parse_input(filename: str) -> ParsedInput:
    try:
        with open(filename, "r", encoding="utf-8") as infile:
            lines = infile.readlines()
    except OSError as exc:
        raise RuntimeError(f"Error opening file: {filename}") from exc

    if not lines:
        raise RuntimeError(f"Empty input file: {filename}")

    out = ParsedInput()
    out.balance_factor = float(lines[0].strip())

    cell_id: dict[str, int] = {}
    tokens = iter(" ".join(lines[1:]).split())

    for tok in tokens:
        if tok != "NET":
            continue

        try:
            net_name = next(tokens)
        except StopIteration as exc:
            raise RuntimeError("Malformed NET line: missing net name") from exc

        net = Net(net_name, len(out.nets))

        for tok in tokens:
            if tok == ";":
                break

            ends_with_semi = tok.endswith(";")
            if ends_with_semi:
                tok = tok[:-1]

            if tok:
                if tok not in cell_id:
                    cid = len(out.cells)
                    cell_id[tok] = cid
                    out.cells.append(Cell(tok, cid))
                else:
                    cid = cell_id[tok]
                net.connected_cells.append(cid)

            if ends_with_semi:
                break

        net.connected_cells = sorted(set(net.connected_cells))

        for pin, cid in enumerate(net.connected_cells):
            out.cells[cid].connected_nets.append(net.idx)
            out.cells[cid].net_pin_idx.append(pin)

        out.nets.append(net)

    return out


def parse_solution(filename: str) -> ParsedSolution:
    try:
        with open(filename, "r", encoding="utf-8") as infile:
            tokens = infile.read().split()
    except OSError as exc:
        raise RuntimeError(f"Error opening solution file: {filename}") from exc

    pos = 0

    def next_token() -> str:
        nonlocal pos
        if pos >= len(tokens):
            raise RuntimeError("Malformed solution: unexpected end of file")
        token = tokens[pos]
        pos += 1
        return token

    sol = ParsedSolution()

    if next_token() != "Cutsize":
        raise RuntimeError("Malformed solution: expected 'Cutsize'")
    if next_token() != "=":
        raise RuntimeError("Malformed solution: expected '=' after 'Cutsize'")
    try:
        sol.declared_cutsize = int(next_token())
    except ValueError as exc:
        raise RuntimeError("Malformed solution: missing cutsize value") from exc

    for _ in range(3):
        group = SolutionGroup()
        group.label = next_token()
        try:
            group.declared_size = int(next_token())
        except ValueError as exc:
            raise RuntimeError(
                f"Malformed solution: missing declared group size for {group.label}"
            ) from exc

        found_semicolon = False
        while pos < len(tokens):
            token = next_token()
            ends_with_semi = token.endswith(";")
            if ends_with_semi:
                token = token[:-1]

            if token:
                group.cells.append(token)

            if ends_with_semi:
                found_semicolon = True
                break

        if not found_semicolon:
            raise RuntimeError(f"Malformed solution: missing ';' after {group.label}")

        sol.groups.append(group)

    if pos != len(tokens):
        raise RuntimeError("Malformed solution: extra content after G3")

    return sol


def balance_limits(n: int, r: float) -> tuple[int, int]:
    if not (0.0 < r < 1.0):
        raise RuntimeError("Balance factor r must satisfy 0 < r < 1")

    min_size = math.ceil(((1.0 - r) / 3.0) * float(n))
    max_size = math.floor(((1.0 + r) / 3.0) * float(n))

    if min_size > max_size:
        min_size = n // 3
        max_size = min_size

    return min_size, max_size


def compute_cutsize(parsed: ParsedInput, cell_group: list[int]) -> int:
    cutsize = 0
    for net in parsed.nets:
        active_groups: set[int] = set()
        for cid in net.connected_cells:
            gid = cell_group[cid]
            active_groups.add(gid)
            if len(active_groups) >= 2:
                cutsize += 1
                break
    return cutsize


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print("Usage: ./solutionchecker.py <input file> <solution file>", file=sys.stderr)
        return 1

    try:
        parsed = parse_input(argv[1])
        solution = parse_solution(argv[2])

        if len(solution.groups) != 3:
            raise RuntimeError("Malformed solution: expected exactly 3 groups")

        cell_id = {cell.name: cell.idx for cell in parsed.cells}
        cell_group = [-1] * len(parsed.cells)
        errors: list[str] = []

        for gid, group in enumerate(solution.groups):
            expected_label = f"G{gid + 1}"
            if group.label != expected_label:
                errors.append(
                    f"Expected group label {expected_label} but found {group.label}"
                )

            if group.declared_size != len(group.cells):
                errors.append(
                    f"{group.label} declared size mismatch: declared {group.declared_size}, "
                    f"actual {len(group.cells)}"
                )

            seen_in_group: set[str] = set()
            for name in group.cells:
                if name in seen_in_group:
                    errors.append(f"Duplicate cell {name} inside {group.label}")
                    continue
                seen_in_group.add(name)

                if name not in cell_id:
                    errors.append(f"Unknown cell {name} in {group.label}")
                    continue

                cid = cell_id[name]
                if cell_group[cid] != -1:
                    errors.append(f"Cell {name} appears in multiple groups")
                    continue

                cell_group[cid] = gid

        missing_cells = 0
        for cell in parsed.cells:
            if cell_group[cell.idx] == -1:
                missing_cells += 1
                errors.append(f"Missing cell {cell.name} from all groups")

        min_size, max_size = balance_limits(len(parsed.cells), parsed.balance_factor)
        for gid, group in enumerate(solution.groups):
            size = len(group.cells)
            if size < min_size or size > max_size:
                errors.append(
                    f"Balance violation on G{gid + 1}: size={size}, "
                    f"allowed=[{min_size}, {max_size}]"
                )

        if missing_cells == 0:
            actual_cutsize = compute_cutsize(parsed, cell_group)
            if actual_cutsize != solution.declared_cutsize:
                errors.append(
                    f"Cutsize mismatch: declared {solution.declared_cutsize}, "
                    f"actual {actual_cutsize}"
                )

        if errors:
            print("FAIL")
            for err in errors:
                print(err)
            return 1

        print("PASS")
        print(f"Cutsize = {solution.declared_cutsize}")
        print(f"Balance range = [{min_size}, {max_size}]")
        print(
            f"G1 = {len(solution.groups[0].cells)}, "
            f"G2 = {len(solution.groups[1].cells)}, "
            f"G3 = {len(solution.groups[2].cells)}"
        )
        return 0
    except Exception as exc:
        print(str(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))

#!/usr/bin/env python3
import argparse
import bisect
import html
import math
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation


@dataclass
class Instance:
    name: str
    x: int
    y: int
    w: int
    h: int
    orient: str
    typ: str


@dataclass
class Placement:
    orient: str
    x: int
    y: int


def ceil_div(a, b):
    return -((-a) // b)


def rect_overlap(a, b):
    return a[0] < b[2] and b[0] < a[2] and a[1] < b[3] and b[1] < a[3]


def rect_overlap_area(a, b):
    ox = max(0, min(a[2], b[2]) - max(a[0], b[0]))
    oy = max(0, min(a[3], b[3]) - max(a[1], b[1]))
    return ox * oy


REQUIRED_GP_COLUMNS = ("Name", "LLX", "LLY", "Width", "Height", "Type")
DBU_EPSILON = Decimal("0.000001")
DEFAULT_ALPHA = 0.7
DEFAULT_THRESHOLD = 45.0
DEFAULT_NORM_FACTOR = 18.2
EDGE_BLOCKAGE_SPAN_RATIO = 0.999
SUPPORTED_ORIENTS = {"R0", "MX"}
SVG_MAX_SITE_LINES = 200
SVG_CANVAS_WIDTH = 1200
SVG_MIN_CANVAS_HEIGHT = 300
SVG_MAX_CANVAS_HEIGHT = 1200
SVG_MAX_ORIGINAL_OUTLINES = 5000


def normalize_input_orient(orient):
    return orient if orient in SUPPORTED_ORIENTS else "R0"


def parse_gp(path):
    dbu = die_llx = die_lly = die_urx = die_ury = site_w = site_h = None
    instances = []
    col = None
    column_count = None
    seen_names = set()

    with open(path, "r", encoding="utf-8") as f:
        for lineno, raw in enumerate(f, 1):
            line = raw.strip()
            if not line:
                continue
            parts = line.split()
            if col is None:
                if parts[0] == "DBU_Per_Micron":
                    dbu = int(parts[1])
                elif parts[0] == "DieArea_LL":
                    die_llx, die_lly = int(parts[1]), int(parts[2])
                elif parts[0] == "DieArea_UR":
                    die_urx, die_ury = int(parts[1]), int(parts[2])
                elif parts[0] == "Site_Width":
                    site_w = int(parts[1])
                elif parts[0] == "Site_Height":
                    site_h = int(parts[1])
                elif parts[0] == "Name":
                    missing = [name for name in REQUIRED_GP_COLUMNS if name not in parts]
                    if missing:
                        raise ValueError(f"line {lineno}: missing GP columns: {', '.join(missing)}")
                    col = {name: parts.index(name) for name in REQUIRED_GP_COLUMNS}
                    if "Orient" in parts:
                        col["Orient"] = parts.index("Orient")
                    column_count = len(parts)
                continue

            if len(parts) < 6:
                raise ValueError(f"line {lineno}: invalid GP row field count")
            orient_col = col.get("Orient")
            type_col = col["Type"]
            if orient_col is not None and type_col > orient_col and len(parts) == column_count - 1:
                orient = "R0"
                typ = parts[type_col - 1]
            elif len(parts) == column_count:
                orient = parts[orient_col] if orient_col is not None else "R0"
                typ = parts[type_col]
            elif len(parts) == column_count + 1 and orient_col is None:
                orient = parts[type_col]
                typ = parts[type_col + 1]
            else:
                orient = parts[5] if len(parts) >= 7 else "R0"
                typ = parts[-1]
            inst = Instance(
                parts[col["Name"]],
                int(parts[col["LLX"]]),
                int(parts[col["LLY"]]),
                int(parts[col["Width"]]),
                int(parts[col["Height"]]),
                normalize_input_orient(orient),
                typ,
            )
            if inst.name in seen_names:
                raise ValueError(f"line {lineno}: duplicate instance name {inst.name}")
            seen_names.add(inst.name)
            instances.append(inst)

    if None in (dbu, die_llx, die_lly, die_urx, die_ury, site_w, site_h):
        raise ValueError("input .gp header is incomplete")
    if col is None:
        raise ValueError("input .gp table header is missing")
    if dbu <= 0 or site_w <= 0 or site_h <= 0:
        raise ValueError("DBU and site dimensions must be positive")
    if die_llx >= die_urx or die_lly >= die_ury:
        raise ValueError("die area must have positive width and height")

    return {
        "dbu": dbu,
        "die": (die_llx, die_lly, die_urx, die_ury),
        "site_w": site_w,
        "site_h": site_h,
        "instances": instances,
    }


COORD_RE = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
TCL_RE = re.compile(
    rf"^\s*place_cell\s+-inst_name\s+(\S+)\s+-orient\s+(\S+)\s+-origin\s+\{{\s*({COORD_RE})\s+({COORD_RE})\s*\}}\s*$"
)


def micron_to_dbu(text, dbu):
    try:
        value = Decimal(text) * Decimal(dbu)
    except InvalidOperation as exc:
        raise ValueError(f"{text} is not a valid coordinate") from exc
    nearest = value.to_integral_value()
    if abs(value - nearest) > DBU_EPSILON:
        raise ValueError(f"{text} um is not on the DBU grid")
    return int(nearest)


def parse_tcl(path, dbu):
    placements = {}
    errors = []
    with open(path, "r", encoding="utf-8") as f:
        for lineno, raw in enumerate(f, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if line.split(None, 1)[0] == "detailed_placement":
                errors.append(f"line {lineno}: detailed_placement is prohibited")
                continue
            m = TCL_RE.match(line)
            if not m:
                errors.append(f"line {lineno}: invalid TCL format")
                continue
            name = m.group(1)
            orient = m.group(2)
            try:
                x = micron_to_dbu(m.group(3), dbu)
                y = micron_to_dbu(m.group(4), dbu)
            except ValueError as exc:
                errors.append(f"line {lineno}: {exc}")
                continue
            if name in placements:
                errors.append(f"line {lineno}: duplicate placement for {name}")
            placements[name] = Placement(orient, x, y)
    return placements, errors


def placed_lower_left(inst, placement):
    if placement is None:
        return inst.x, inst.y
    if placement.orient == "MX":
        return placement.x, placement.y - inst.h
    return placement.x, placement.y


def check_legality(data, placements):
    die_llx, die_lly, die_urx, die_ury = data["die"]
    site_w, site_h = data["site_w"], data["site_h"]
    cells = [inst for inst in data["instances"] if inst.typ == "CELL"]
    obstacles = [inst for inst in data["instances"] if inst.typ != "CELL"]
    cell_by_name = {cell.name: cell for cell in cells}
    inst_by_name = {inst.name: inst for inst in data["instances"]}
    errors = []
    row_count = ceil_div(die_ury - die_lly, site_h)

    placement_names = set(placements)
    cell_names = set(cell_by_name)
    known_names = set(inst_by_name)

    unknown = sorted(placement_names - known_names)
    if unknown:
        errors.append(f"unknown placements: {len(unknown)}")
    non_cell = sorted((placement_names & known_names) - cell_names)
    if non_cell:
        errors.append(f"non-cell placements: {len(non_cell)}")
    missing = sorted(cell_names - placement_names)
    if missing:
        errors.append(f"missing cell placements: {len(missing)}")

    row_intervals = defaultdict(list)
    for name, cell in cell_by_name.items():
        if name not in placements:
            continue
        placement = placements[name]
        if placement.orient not in SUPPORTED_ORIENTS:
            errors.append(f"{name}: unsupported orientation {placement.orient}")
        if placement.orient != cell.orient:
            errors.append(f"{name}: orientation changed from {cell.orient} to {placement.orient}")
        x, y = placed_lower_left(cell, placement)
        rect = (x, y, x + cell.w, y + cell.h)

        if x < die_llx or y < die_lly or x + cell.w > die_urx or y + cell.h > die_ury:
            errors.append(f"{name}: outside die")
        if (x - die_llx) % site_w != 0:
            errors.append(f"{name}: x is not site-aligned")
        if (y - die_lly) % site_h != 0:
            errors.append(f"{name}: y is not row-aligned")

        for obs in obstacles:
            obs_rect = (obs.x, obs.y, obs.x + obs.w, obs.y + obs.h)
            if rect_overlap(rect, obs_rect):
                errors.append(f"{name}: overlaps {obs.typ} {obs.name}")
                break

        row_start = max(0, (max(y, die_lly) - die_lly) // site_h)
        row_end = min(row_count, ceil_div(min(y + cell.h, die_ury) - die_lly, site_h))
        for row in range(row_start, row_end):
            row_intervals[row].append((x, x + cell.w, name, rect))

    for row, intervals in row_intervals.items():
        intervals.sort()
        active = []
        for x1, x2, name, rect in intervals:
            active = [item for item in active if item[0] > x1]
            overlap_name = None
            for _, other_name, other_rect in active:
                if rect_overlap(rect, other_rect):
                    overlap_name = other_name
                    break
            if overlap_name is not None:
                errors.append(f"{name}: overlaps cell {overlap_name} on row {row}")
                break
            active.append((x2, name, rect))

    return errors


def displacement_metrics(data, placements):
    dbu = data["dbu"]
    total_dist = 0
    max_dist = 0
    moved_count = 0

    for inst in data["instances"]:
        placement = placements.get(inst.name)
        x, y = placed_lower_left(inst, placement)
        dist = abs(x - inst.x) + abs(y - inst.y)
        total_dist += dist
        max_dist = max(max_dist, dist)
        moved_count += 1

    if moved_count > 0:
        total_u = total_dist / float(dbu)
        avg_u = total_u / moved_count
        max_u = max_dist / float(dbu)
    else:
        total_u = 0.0
        avg_u = 0.0
        max_u = 0.0

    return total_u, avg_u, max_u, moved_count


def infer_heatmap_bounds(data):
    die_llx, die_lly, die_urx, die_ury = data["die"]
    core_llx, core_lly, core_urx, core_ury = die_llx, die_lly, die_urx, die_ury
    die_w = die_urx - die_llx
    die_h = die_ury - die_lly

    for inst in data["instances"]:
        if inst.typ != "BLOCKAGE":
            continue
        x1 = max(inst.x, die_llx)
        y1 = max(inst.y, die_lly)
        x2 = min(inst.x + inst.w, die_urx)
        y2 = min(inst.y + inst.h, die_ury)
        if x1 >= x2 or y1 >= y2:
            continue

        spans_die_width = (x2 - x1) / die_w >= EDGE_BLOCKAGE_SPAN_RATIO
        spans_die_height = (y2 - y1) / die_h >= EDGE_BLOCKAGE_SPAN_RATIO
        if y1 <= die_lly and spans_die_width:
            core_lly = max(core_lly, y2)
        if y2 >= die_ury and spans_die_width:
            core_ury = min(core_ury, y1)
        if x1 <= die_llx and spans_die_height:
            core_llx = max(core_llx, x2)
        if x2 >= die_urx and spans_die_height:
            core_urx = min(core_urx, x1)

    if core_llx >= core_urx or core_lly >= core_ury:
        return data["die"]
    return core_llx, core_lly, core_urx, core_ury


def density_overflow_ratio(data, placements, threshold):
    dbu = data["dbu"]
    core_llx, core_lly, core_urx, core_ury = infer_heatmap_bounds(data)
    grid = 10 * dbu
    x_grid = [core_llx]
    for col in range(ceil_div(core_urx - core_llx, grid)):
        x_grid.append(min(core_llx + (col + 1) * grid, core_urx))
    y_grid = [core_lly]
    for row in range(ceil_div(core_ury - core_lly, grid)):
        y_grid.append(min(core_lly + (row + 1) * grid, core_ury))
    cols = len(x_grid) - 1
    rows = len(y_grid) - 1
    used = [[0.0 for _ in range(cols)] for _ in range(rows)]
    has_value = [[False for _ in range(cols)] for _ in range(rows)]

    def map_view(x1, y1, x2, y2):
        c0 = max(bisect.bisect_left(x_grid, x1) - 1, 0)
        c1 = min(bisect.bisect_right(x_grid, x2), cols)
        r0 = max(bisect.bisect_left(y_grid, y1) - 1, 0)
        r1 = min(bisect.bisect_right(y_grid, y2), rows)
        return range(r0, r1), range(c0, c1)

    def grid_rect(row, col):
        return x_grid[col], y_grid[row], x_grid[col + 1], y_grid[row + 1]

    # Match the scoring Tcl: destroy h* macros before dumping the Placement heatmap.
    # The OpenROAD Placement heatmap reports standard-cell placement density, not
    # placement blockages. addToMap() marks every bin returned by getMapView() as
    # having data, even when the exact overlap area is zero on a grid boundary.
    cells = [inst for inst in data["instances"] if inst.typ == "CELL" and not inst.name.startswith("h")]
    for cell in cells:
        placement = placements.get(cell.name)
        x, y = placed_lower_left(cell, placement)
        rect = (x, y, x + cell.w, y + cell.h)
        row_range, col_range = map_view(*rect)
        for r in row_range:
            for c in col_range:
                has_value[r][c] = True
                used[r][c] += rect_overlap_area(grid_rect(r, c), rect)

    total = 0
    overflow = 0
    for r in range(rows):
        for c in range(cols):
            if not has_value[r][c]:
                continue
            gx1, gy1, gx2, gy2 = grid_rect(r, c)
            grid_area = (gx2 - gx1) * (gy2 - gy1)
            density = used[r][c] / grid_area * 100.0
            total += 1
            if density > threshold:
                overflow += 1
    return (overflow / total * 100.0) if total else 0.0, overflow, total


def is_number(text):
    try:
        return math.isfinite(float(text))
    except ValueError:
        return False


def parse_float(parser, text, name):
    try:
        value = float(text)
    except ValueError:
        parser.error(f"{name} must be numeric")
    if not math.isfinite(value):
        parser.error(f"{name} must be finite")
    return value


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Check placement legality and report the flow.tcl-compatible quality score.",
        usage=(
            "%(prog)s [--plot] <input.gp> <output.tcl> [alpha] [threshold] [norm_factor]\n"
            "       %(prog)s [--plot] <alpha> <threshold> <input.gp> <output.tcl> [norm_factor]"
        ),
    )
    parser.add_argument("args", nargs="+")
    parser.add_argument(
        "--heat-csv",
        help="optional OpenROAD gui::dump_heatmap Placement CSV; defaults to <designName>_heat.csv if present",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="write <designName>_placement.svg next to the output TCL",
    )
    ns = parser.parse_args(argv)
    values = ns.args

    if len(values) == 2:
        input_gp, output_tcl = values
        alpha, threshold, norm_factor = DEFAULT_ALPHA, DEFAULT_THRESHOLD, DEFAULT_NORM_FACTOR
    elif len(values) == 3:
        input_gp, output_tcl = values[:2]
        alpha = parse_float(parser, values[2], "alpha")
        threshold, norm_factor = DEFAULT_THRESHOLD, DEFAULT_NORM_FACTOR
    elif len(values) == 4 and is_number(values[0]) and is_number(values[1]):
        alpha = parse_float(parser, values[0], "alpha")
        threshold = parse_float(parser, values[1], "threshold")
        norm_factor = DEFAULT_NORM_FACTOR
        input_gp, output_tcl = values[2:]
    elif len(values) == 4:
        input_gp, output_tcl = values[:2]
        alpha = parse_float(parser, values[2], "alpha")
        threshold = parse_float(parser, values[3], "threshold")
        norm_factor = DEFAULT_NORM_FACTOR
    elif len(values) == 5 and is_number(values[0]) and is_number(values[1]):
        alpha = parse_float(parser, values[0], "alpha")
        threshold = parse_float(parser, values[1], "threshold")
        input_gp, output_tcl = values[2:4]
        norm_factor = parse_float(parser, values[4], "norm_factor")
    elif len(values) == 5:
        input_gp, output_tcl = values[:2]
        alpha = parse_float(parser, values[2], "alpha")
        threshold = parse_float(parser, values[3], "threshold")
        norm_factor = parse_float(parser, values[4], "norm_factor")
    else:
        parser.error("expected either checker order or assignment order arguments")

    if not 0.0 <= alpha <= 1.0:
        parser.error("alpha must be in [0, 1]")
    if threshold < 0.0:
        parser.error("threshold must be non-negative")
    if norm_factor < 0.0:
        parser.error("norm_factor must be non-negative")

    return argparse.Namespace(
        input_gp=input_gp,
        output_tcl=output_tcl,
        alpha=alpha,
        threshold=threshold,
        norm_factor=norm_factor,
        heat_csv=ns.heat_csv,
        plot=ns.plot,
    )


def filename_warnings(input_gp, output_tcl):
    warnings = []
    gp_base = os.path.basename(input_gp)
    tcl_base = os.path.basename(output_tcl)
    gp_suffix = "_insts.gp"
    tcl_suffix = "_insts.tcl"

    if not gp_base.endswith(gp_suffix):
        warnings.append(f"input filename should match <designName>{gp_suffix}")
    if not tcl_base.endswith(tcl_suffix):
        warnings.append(f"output filename should match <designName>{tcl_suffix}")
    if gp_base.endswith(gp_suffix) and tcl_base.endswith(tcl_suffix):
        gp_design = gp_base[: -len(gp_suffix)]
        tcl_design = tcl_base[: -len(tcl_suffix)]
        if gp_design != tcl_design:
            warnings.append("input/output design names do not match")
    return warnings


def design_name_from_gp(input_gp):
    gp_base = os.path.basename(input_gp)
    gp_suffix = "_insts.gp"
    if gp_base.endswith(gp_suffix):
        return gp_base[: -len(gp_suffix)]
    return os.path.splitext(gp_base)[0]


def density_overflow_ratio_from_heat_csv(path, threshold):
    total = 0
    overflow = 0
    with open(path, "r", encoding="utf-8") as f:
        next(f, None)
        for lineno, raw in enumerate(f, 2):
            data = [part.strip() for part in raw.strip().split(",")]
            if len(data) <= 4 or data[4] == "":
                continue
            try:
                density = float(data[4])
            except ValueError as exc:
                raise ValueError(f"{path}:{lineno}: invalid density value {data[4]!r}") from exc
            total += 1
            if density > threshold:
                overflow += 1
    return (overflow / total * 100.0) if total else 0.0, overflow, total


def default_plot_svg_path(input_gp, output_tcl):
    output_dir = os.path.dirname(output_tcl) or "."
    return os.path.join(output_dir, f"{design_name_from_gp(input_gp)}_placement.svg")


def format_svg_number(value):
    return f"{value:.6f}".rstrip("0").rstrip(".")


def svg_rect(x1, y1, x2, y2, die_lly, die_ury, class_name, title, extra_attrs=""):
    svg_y = die_lly + die_ury - y2
    attrs = (
        f'x="{format_svg_number(x1)}" '
        f'y="{format_svg_number(svg_y)}" '
        f'width="{format_svg_number(x2 - x1)}" '
        f'height="{format_svg_number(y2 - y1)}" '
        f'class="{class_name}"'
    )
    if extra_attrs:
        attrs = f"{attrs} {extra_attrs}"
    return f"  <rect {attrs}><title>{html.escape(title)}</title></rect>"


def write_placement_svg(path, data, placements):
    die_llx, die_lly, die_urx, die_ury = data["die"]
    die_w = die_urx - die_llx
    die_h = die_ury - die_lly
    canvas_h = int(round(SVG_CANVAS_WIDTH * die_h / die_w)) if die_w else SVG_MIN_CANVAS_HEIGHT
    canvas_h = max(SVG_MIN_CANVAS_HEIGHT, min(SVG_MAX_CANVAS_HEIGHT, canvas_h))
    site_h = data["site_h"]
    row_count = ceil_div(die_h, site_h)

    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        (
            '<svg xmlns="http://www.w3.org/2000/svg" '
            f'width="{SVG_CANVAS_WIDTH}" height="{canvas_h}" '
            f'viewBox="{format_svg_number(die_llx)} {format_svg_number(die_lly)} '
            f'{format_svg_number(die_w)} {format_svg_number(die_h)}">'
        ),
        "  <style>",
        "    .background { fill: #ffffff; }",
        "    .die { fill: none; stroke: #111827; stroke-width: 2; vector-effect: non-scaling-stroke; }",
        "    .row { stroke: #d1d5db; stroke-width: 1; vector-effect: non-scaling-stroke; }",
        "    .cell { fill: #60a5fa; fill-opacity: 0.72; stroke: #1d4ed8; stroke-width: 1; vector-effect: non-scaling-stroke; }",
        "    .moved { fill: #34d399; fill-opacity: 0.72; stroke: #047857; stroke-width: 1; vector-effect: non-scaling-stroke; }",
        "    .missing { fill: #f87171; fill-opacity: 0.62; stroke: #b91c1c; stroke-width: 1; vector-effect: non-scaling-stroke; }",
        "    .original { fill: none; stroke: #f59e0b; stroke-width: 1; stroke-dasharray: 4 3; vector-effect: non-scaling-stroke; }",
        "    .blockage { fill: #6b7280; fill-opacity: 0.55; stroke: #374151; stroke-width: 1; vector-effect: non-scaling-stroke; }",
        "    .fixed { fill: #c084fc; fill-opacity: 0.5; stroke: #7e22ce; stroke-width: 1; vector-effect: non-scaling-stroke; }",
        "  </style>",
        (
            f'  <rect x="{format_svg_number(die_llx)}" y="{format_svg_number(die_lly)}" '
            f'width="{format_svg_number(die_w)}" height="{format_svg_number(die_h)}" '
            'class="background"/>'
        ),
        svg_rect(die_llx, die_lly, die_urx, die_ury, die_lly, die_ury, "die", "Die area"),
    ]

    if 0 < row_count <= SVG_MAX_SITE_LINES:
        for row in range(1, row_count):
            design_y = die_lly + row * site_h
            if design_y >= die_ury:
                break
            svg_y = die_lly + die_ury - design_y
            lines.append(
                f'  <line x1="{format_svg_number(die_llx)}" y1="{format_svg_number(svg_y)}" '
                f'x2="{format_svg_number(die_urx)}" y2="{format_svg_number(svg_y)}" class="row"/>'
            )

    cells = []
    obstacles = []
    for inst in data["instances"]:
        if inst.typ == "CELL":
            cells.append(inst)
        else:
            obstacles.append(inst)

    for inst in obstacles:
        x1, y1 = inst.x, inst.y
        x2, y2 = inst.x + inst.w, inst.y + inst.h
        class_name = "blockage" if inst.typ == "BLOCKAGE" else "fixed"
        title = f"{inst.name} {inst.typ} ({x1}, {y1})-({x2}, {y2})"
        lines.append(svg_rect(x1, y1, x2, y2, die_lly, die_ury, class_name, title))

    for cell in cells:
        placement = placements.get(cell.name)
        x, y = placed_lower_left(cell, placement)
        x1, y1 = x, y
        x2, y2 = x + cell.w, y + cell.h
        moved = placement is not None and (x != cell.x or y != cell.y)
        if moved and len(cells) <= SVG_MAX_ORIGINAL_OUTLINES:
            orig_x2 = cell.x + cell.w
            orig_y2 = cell.y + cell.h
            title = f"{cell.name} original ({cell.x}, {cell.y})-({orig_x2}, {orig_y2})"
            lines.append(svg_rect(cell.x, cell.y, orig_x2, orig_y2, die_lly, die_ury, "original", title))
        class_name = "missing" if placement is None else ("moved" if moved else "cell")
        title = f"{cell.name} CELL {cell.orient} ({x1}, {y1})-({x2}, {y2})"
        lines.append(svg_rect(x1, y1, x2, y2, die_lly, die_ury, class_name, title))

    lines.append("</svg>")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
        f.write("\n")


def main(argv=None):
    args = parse_args(argv)

    data = parse_gp(args.input_gp)
    placements, parse_errors = parse_tcl(args.output_tcl, data["dbu"])
    legality_errors = parse_errors + check_legality(data, placements)
    total_u, avg_u, max_u, moved_count = displacement_metrics(data, placements)
    heat_csv = args.heat_csv
    if heat_csv is None:
        candidate = os.path.join(
            os.path.dirname(args.output_tcl) or ".",
            f"{design_name_from_gp(args.input_gp)}_heat.csv",
        )
        if os.path.exists(candidate):
            heat_csv = candidate
    if heat_csv is None:
        dor, overflow, total_grids = density_overflow_ratio(data, placements, args.threshold)
    else:
        dor, overflow, total_grids = density_overflow_ratio_from_heat_csv(heat_csv, args.threshold)
    norm_disp = avg_u * args.norm_factor
    quality = args.alpha * norm_disp + (1.0 - args.alpha) * dor
    warnings = filename_warnings(args.input_gp, args.output_tcl)
    plot_svg = default_plot_svg_path(args.input_gp, args.output_tcl) if args.plot else None
    if plot_svg is not None:
        write_placement_svg(plot_svg, data, placements)

    print("Performance Metrics")
    print("----------------------------------------")
    print(f"Design name            : {design_name_from_gp(args.input_gp)}")
    print(f"Legalizer output       : {args.output_tcl}")
    print(f"Legality PASS          : {0 if legality_errors else 1}")
    print(f"Moved instances        : {moved_count}")
    print(f"Total displacement     : {total_u:.1f} u")
    print(f"Average displacement   : {avg_u:.1f} u")
    print(f"Max displacement       : {max_u:.1f} u")
    print(f"Threshold              : {args.threshold:g}")
    print(f"Total Grids            : {total_grids}")
    print(f"Overflow Grids         : {overflow}")
    print(f"DOR                    : {dor:.2f} %")
    print(f"Norm. Displacement     : {norm_disp:.2f}")
    if plot_svg is not None:
        print(f"Placement SVG          : {plot_svg}")
    print("----------------------------------------")
    print(f"FINAL QUALITY SCORE    : {quality:.4f}")
    print("----------------------------------------")
    if warnings:
        print("Filename warnings:")
        for warning in warnings:
            print(f"  - {warning}")
    if legality_errors:
        print("Errors:")
        for err in legality_errors[:50]:
            print(f"  - {err}")
        if len(legality_errors) > 50:
            print(f"  ... {len(legality_errors) - 50} more")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

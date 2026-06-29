#!/usr/bin/env python3
"""Verify floorplan output and optionally render symmetry groups with matplotlib."""
import argparse
import sys


PALETTE = [
    "#e76f51",
    "#2a9d8f",
    "#e9c46a",
    "#f4a261",
    "#457b9d",
    "#8d99ae",
    "#ef476f",
    "#06d6a0",
    "#118ab2",
    "#ff7f11",
]


def read_input(path):
    with open(path) as f:
        lines = [l.strip() for l in f if l.strip()]
    it = iter(lines)
    n = int(next(it).split(":")[1])
    orig_dims = {}
    for _ in range(n):
        parts = next(it).split()
        orig_dims[parts[0]] = (float(parts[1]), float(parts[2]))
    groups = []
    cur = None
    for line in it:
        if line == "Symmetry Group":
            if cur:
                groups.append(cur)
            cur = {"pairs": [], "sel": []}
        else:
            parts = line.split()
            if len(parts) == 2:
                cur["pairs"].append((parts[0], parts[1]))
            else:
                cur["sel"].append(parts[0])
    if cur:
        groups.append(cur)
    return orig_dims, groups


def read_output(path):
    pos = {}
    rot = {}
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 3:
                pos[parts[0]] = (float(parts[1]), float(parts[2]))
                orientation = int(parts[3]) if len(parts) >= 4 else 0
                if orientation not in (0, 1):
                    raise ValueError(f"Invalid orientation for {parts[0]}: {orientation}")
                rot[parts[0]] = orientation
    return pos, rot


def load_layout(inp, out):
    orig_dims, groups = read_input(inp)
    pos, rot = read_output(out)
    dims = {}
    for name, (w, h) in orig_dims.items():
        if rot.get(name, 0):
            dims[name] = (h, w)
        else:
            dims[name] = (w, h)
    return orig_dims, groups, pos, dims


def build_rects(dims, pos):
    missing = [name for name in dims if name not in pos]
    if missing:
        names = ", ".join(missing[:8])
        if len(missing) > 8:
            names += ", ..."
        raise ValueError(f"Missing block placements in output: {names}")

    rects = []
    for name in dims:
        w, h = dims[name]
        x, y = pos[name]
        rects.append((name, x, y, x + w, y + h))
    return rects


def infer_group_axis(group, dims, pos, eps=1e-6):
    axis_type = None
    axis_value = None

    for a, b in group["pairs"]:
        wa, ha = dims[a]
        xa, ya = pos[a]
        wb, hb = dims[b]
        xb, yb = pos[b]
        cxa = xa + wa / 2
        cya = ya + ha / 2
        cxb = xb + wb / 2
        cyb = yb + hb / 2

        if abs(cya - cyb) <= eps:
            cur_type = "V"
            cur_axis = (cxa + cxb) / 2
        elif abs(cxa - cxb) <= eps:
            cur_type = "H"
            cur_axis = (cya + cyb) / 2
        else:
            continue

        if axis_type is None:
            axis_type = cur_type
            axis_value = cur_axis
        elif cur_type == axis_type and abs(cur_axis - axis_value) <= eps:
            continue
        else:
            return None, None

    return axis_type, axis_value


def render_plot(inp, out, groups, pos, dims):
    try:
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
        from matplotlib.patches import Rectangle
    except ModuleNotFoundError as exc:
        raise ValueError(
            "matplotlib is not installed for this Python interpreter. " "Run the script from the project .venv."
        ) from exc

    rects = build_rects(dims, pos)
    all_x = []
    all_y = []
    for _, x1, y1, x2, y2 in rects:
        all_x += [x1, x2]
        all_y += [y1, y2]

    min_x = min(all_x)
    max_x = max(all_x)
    min_y = min(all_y)
    max_y = max(all_y)
    span_x = max_x - min_x
    span_y = max_y - min_y

    fig_w = min(22, max(10, (span_x / 90.0) + 4.0))
    fig_h = min(16, max(7, (span_y / 90.0) + 2.5))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=160)

    block_group = {}
    group_members = []
    for gi, group in enumerate(groups):
        members = []
        for a, b in group["pairs"]:
            members.extend([a, b])
        members.extend(group["sel"])
        members = list(dict.fromkeys(members))
        group_members.append(members)
        for name in members:
            block_group[name] = gi

    for name, x, y, _, _ in rects:
        w, h = dims[name]
        gi = block_group.get(name)
        if gi is None:
            facecolor = "#d7dde5"
            edgecolor = "#5b6472"
        else:
            facecolor = PALETTE[gi % len(PALETTE)]
            edgecolor = "#1f2937"

        patch = Rectangle(
            (x, y),
            w,
            h,
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=1.2,
            alpha=0.76,
            zorder=2,
        )
        ax.add_patch(patch)

    ax.add_patch(
        Rectangle(
            (min_x, min_y),
            span_x,
            span_y,
            fill=False,
            edgecolor="#111827",
            linewidth=1.4,
            zorder=4,
        )
    )

    legend_handles = []
    legend_labels = []
    for gi, _ in enumerate(group_members):
        color = PALETTE[gi % len(PALETTE)]
        axis_type, _ = infer_group_axis(groups[gi], dims, pos, eps=0.6)
        axis_text = "vertical" if axis_type == "V" else "horizontal" if axis_type == "H" else "unresolved"
        legend_handles.append(Line2D([0], [0], color=color, linewidth=6))
        legend_labels.append(f"G{gi + 1} ({axis_text})")

    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_y, max_y)
    ax.set_aspect("equal", adjustable="box")
    ax.margins(x=0, y=0)
    ax.set_title(f"Floorplan Symmetry Groups\n{inp} -> {out}", fontsize=12)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True, linestyle=":", linewidth=0.6, alpha=0.35)

    if legend_handles:
        ax.legend(
            legend_handles,
            legend_labels,
            title="Symmetry Groups",
            loc="upper left",
            bbox_to_anchor=(1.02, 1.0),
            borderaxespad=0.0,
            fontsize=8,
            title_fontsize=9,
        )

    fig.tight_layout()
    plt.show()
    plt.close(fig)


def check_layout(groups, pos, dims):
    rects = build_rects(dims, pos)
    area = None

    ok = True
    eps = 1e-6

    # --- Overlap check (strict interior overlap, touching is allowed) ---
    for i in range(len(rects)):
        for j in range(i + 1, len(rects)):
            n1, x1, y1, x2, y2 = rects[i]
            n2, x3, y3, x4, y4 = rects[j]
            ox = min(x2, x4) - max(x1, x3)
            oy = min(y2, y4) - max(y1, y3)
            if ox > eps and oy > eps:
                print(f"  OVERLAP: {n1}[{x1},{y1},{x2},{y2}]  {n2}[{x3},{y3},{x4},{y4}]" f"  (ox={ox}, oy={oy})")
                ok = False

    # --- Symmetry check ---
    for gi, g in enumerate(groups):
        axes = []
        sym_type = None

        for a, b in g["pairs"]:
            wa, ha = dims[a]
            xa, ya = pos[a]
            wb, hb = dims[b]
            xb, yb = pos[b]
            assert wa == wb and ha == hb, f"Pair dims mismatch: {a}={wa}x{ha}, {b}={wb}x{hb}"
            cxa = xa + wa / 2
            cya = ya + ha / 2
            cxb = xb + wb / 2
            cyb = yb + hb / 2
            if abs(cya - cyb) <= eps:  # vertical sym: same y-center
                ax = (cxa + cxb) / 2
                axes.append(("V", ax, f"pair({a},{b})"))
                if sym_type is None:
                    sym_type = "V"
            elif abs(cxa - cxb) <= eps:  # horizontal sym: same x-center
                ay = (cya + cyb) / 2
                axes.append(("H", ay, f"pair({a},{b})"))
                if sym_type is None:
                    sym_type = "H"
            else:
                print(
                    f"  Group {gi}: PAIR {a},{b} not symmetric  "
                    f"centers=({cxa:.3f},{cya:.3f}) ({cxb:.3f},{cyb:.3f})"
                )
                ok = False

        for s in g["sel"]:
            ws, hs = dims[s]
            xs, ys = pos[s]
            cxs = xs + ws / 2
            cys = ys + hs / 2
            axes.append(("S", cxs, cys, f"sel({s})"))

        # Check consistency: all entries agree on axis
        ref_ax = None
        ref_type = sym_type
        for entry in axes:
            if entry[0] in ("V", "H"):
                if ref_ax is None:
                    ref_ax = entry[1]
                elif abs(entry[1] - ref_ax) > eps:
                    print(f"  Group {gi}: axis mismatch {entry[2]}  axis={entry[1]:.3f}  ref={ref_ax:.3f}")
                    ok = False
            else:  # selfsym
                _, cxs, cys, label = entry
                if ref_type == "V" and ref_ax is not None and abs(cxs - ref_ax) > eps:
                    print(f"  Group {gi}: {label} center_x={cxs:.3f} != axis={ref_ax:.3f}")
                    ok = False
                elif ref_type == "H" and ref_ax is not None and abs(cys - ref_ax) > eps:
                    print(f"  Group {gi}: {label} center_y={cys:.3f} != axis={ref_ax:.3f}")
                    ok = False

    if ok:
        all_x = []
        all_y = []
        for name in dims:
            w, h = dims[name]
            x, y = pos[name]
            all_x += [x, x + w]
            all_y += [y, y + h]
        area = (max(all_x) - min(all_x)) * (max(all_y) - min(all_y))
        print(f"OK  area={area}")
    # print (all blocks area sum / total area)
    if area is not None:
        print(f"  Utilization: {sum(w * h for w, h in dims.values()) / area:.2%}")
    return ok


def check(inp, out):
    _, groups, pos, dims = load_layout(inp, out)
    return check_layout(groups, pos, dims)


def parse_args(argv):
    parser = argparse.ArgumentParser(description="Verify floorplan output and optionally render symmetry groups.")
    parser.add_argument("input", help="input floorplan file")
    parser.add_argument("output", help="output placement file")
    parser.add_argument(
        "--plot",
        action="store_true",
        help="show a matplotlib window with symmetry-colored blocks",
    )
    return parser.parse_args(argv)


if __name__ == "__main__":
    try:
        args = parse_args(sys.argv[1:])
        _, groups, pos, dims = load_layout(args.input, args.output)
        check_layout(groups, pos, dims)

        if args.plot:
            render_plot(args.input, args.output, groups, pos, dims)
    except (AssertionError, FileNotFoundError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

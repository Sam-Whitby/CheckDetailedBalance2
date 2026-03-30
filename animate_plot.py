#!/usr/bin/env python3
"""
animate_plot.py  --  Animated matplotlib visualisation for DetailedBalanceChecker

Called by animate.wls:
    python3 animate_plot.py <data.json>

Displays two panels:
  Left:  imshow of the lattice at each recorded step.
         Colour encodes particle type (0=empty=dark, 1,2,...=distinct colours).
  Right: system energy vs step (line grows as animation progresses).
"""

import json
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch


def load_data(path: str) -> dict:
    with open(path) as f:
        return json.load(f)


def make_colormap(n_types: int):
    """
    Build a discrete ListedColormap:
      index 0  → near-black (empty site)
      index k  → tab10 colour k-1  (particle type k)
    """
    # tab10 palette: 10 perceptually distinct colours
    tab10 = [
        "#4c72b0",  # blue
        "#dd8452",  # orange
        "#55a868",  # green
        "#c44e52",  # red
        "#8172b3",  # purple
        "#937860",  # brown
        "#da8bc3",  # pink
        "#8c8c8c",  # grey
        "#ccb974",  # olive
        "#64b5cd",  # cyan
    ]
    colors = ["#1a1a1a"] + [tab10[i % len(tab10)] for i in range(n_types)]
    return ListedColormap(colors)


def make_simple_colormap():
    """2-colour map: index 0 → hole (light grey), index 1 → particle (blue)."""
    return ListedColormap(["#dddddd", "#4c72b0"])


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 animate_plot.py <data.json>")
        sys.exit(1)

    d = load_data(sys.argv[1])

    steps     = d["steps"]
    states    = d["states"]
    energies  = d["energies"]
    n_rows    = d["grid_rows"]
    n_cols    = d["grid_cols"]
    n_types   = max(int(d.get("n_types", 1)), 1)
    fps       = float(d.get("fps", 1000.0 / int(d.get("delay_ms", 100))))
    delay_ms  = 1000.0 / fps          # interval passed to FuncAnimation
    simple    = bool(d.get("simple", False))
    algo_name = d.get("algo_file", "")
    n_frames  = len(steps)
    params    = d.get("params", {})

    if simple:
        cmap = make_simple_colormap()
        vmin, vmax = -0.5, 1.5   # 0 = hole, 1 = any particle
    else:
        cmap = make_colormap(n_types)
        vmin = -0.5
        vmax = n_types + 0.5

    # ------------------------------------------------------------------ #
    # Figure layout
    # ------------------------------------------------------------------ #
    fig, (ax_grid, ax_energy) = plt.subplots(
        1, 2, figsize=(12, max(4, n_rows + 1)),
        gridspec_kw={"width_ratios": [1, 1.6]})
    fig.patch.set_facecolor("#f5f5f5")

    # Title: algorithm file; subtitle: parameter values
    mode_tag = f"  [simple {fps:.0f} fps]" if simple else f"  [{fps:.1f} fps]"
    fig.suptitle(algo_name + mode_tag, fontsize=9, y=0.99, style="italic")
    if params:
        param_str = "   ".join(
            f"β={v:.3f}" if k == "beta" else f"{k}={v:.3f}"
            for k, v in params.items())
        fig.text(0.5, 0.955, param_str,
                 ha="center", va="top", fontsize=8, color="#444444",
                 fontfamily="monospace")

    # ---- Lattice panel ------------------------------------------------ #
    def to_grid(state):
        arr = np.array(state, dtype=float).reshape(n_rows, n_cols)
        if simple:
            arr = (arr > 0).astype(float)   # collapse all types → 1
        return arr

    grid0 = to_grid(states[0])
    im = ax_grid.imshow(
        grid0, cmap=cmap, vmin=vmin, vmax=vmax,
        interpolation="nearest", aspect="equal")

    # Tick every site
    ax_grid.set_xticks(range(n_cols))
    ax_grid.set_xticklabels(range(1, n_cols + 1), fontsize=7)
    if n_rows > 1:
        ax_grid.set_yticks(range(n_rows))
        ax_grid.set_yticklabels(range(1, n_rows + 1), fontsize=7)
        ax_grid.set_ylabel("Row", fontsize=9)
    else:
        ax_grid.set_yticks([])
    ax_grid.set_xlabel("Site", fontsize=9)
    ax_grid.set_title("Step 0", fontsize=10)

    # Grid lines between cells
    ax_grid.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax_grid.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax_grid.grid(which="minor", color="white", linewidth=1.5)
    ax_grid.tick_params(which="minor", length=0)

    # Colourbar legend
    if simple:
        legend_handles = [
            Patch(facecolor="#dddddd", label="0  (hole)"),
            Patch(facecolor="#4c72b0", label="n>0  (particle)"),
        ]
    else:
        legend_handles = [Patch(facecolor="#1a1a1a", label="0  (empty)")]
        tab10_list = ["#4c72b0","#dd8452","#55a868","#c44e52","#8172b3",
                      "#937860","#da8bc3","#8c8c8c","#ccb974","#64b5cd"]
        for k in range(1, n_types + 1):
            legend_handles.append(
                Patch(facecolor=tab10_list[(k-1) % len(tab10_list)],
                      label=f"{k}  (type {k})"))
    ax_grid.legend(handles=legend_handles, loc="upper left",
                   bbox_to_anchor=(1.01, 1), fontsize=8,
                   framealpha=0.9, title="Particle", title_fontsize=8)

    # ---- Energy panel ------------------------------------------------- #
    ax_energy.set_facecolor("#fafafa")
    emin = min(energies)
    emax = max(energies)
    pad  = max((emax - emin) * 0.15, 0.05)
    ax_energy.set_xlim(0, max(steps) if steps else 1)
    ax_energy.set_ylim(emin - pad, emax + pad)
    ax_energy.set_xlabel("Step", fontsize=9)
    ax_energy.set_ylabel("Energy", fontsize=9)
    ax_energy.set_title("System energy", fontsize=10)
    ax_energy.grid(True, alpha=0.3, linestyle="--")

    e_line, = ax_energy.plot([], [], color="#4c72b0", lw=1.5, label="E(t)")
    e_dot,  = ax_energy.plot([], [], "o", color="#c44e52", ms=6, zorder=5)
    e_text  = ax_energy.text(
        0.97, 0.05, "", transform=ax_energy.transAxes,
        ha="right", va="bottom", fontsize=8, color="#444444")
    ax_energy.legend(fontsize=8, loc="upper right")

    xs: list = []
    ys: list = []

    # ------------------------------------------------------------------ #
    # Animation update function
    # ------------------------------------------------------------------ #
    def update(frame: int):
        # Lattice
        grid = to_grid(states[frame])
        im.set_data(grid)
        ax_grid.set_title(f"Step {steps[frame]}", fontsize=10)

        # Energy
        xs.append(steps[frame])
        ys.append(energies[frame])
        e_line.set_data(xs, ys)
        e_dot.set_data([xs[-1]], [ys[-1]])
        e_text.set_text(f"E = {energies[frame]:.4f}")

        return im, e_line, e_dot, e_text

    ani = animation.FuncAnimation(
        fig, update, frames=n_frames,
        interval=delay_ms, blit=simple, repeat=False)

    plt.tight_layout(rect=[0, 0, 1, 0.93])
    plt.show()


if __name__ == "__main__":
    main()

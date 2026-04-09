import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Arc, FancyArrowPatch


def annotate_plasmid(results: dict, plasmid_sequence: str, title: str = "Plasmid") -> None:
    """
    Draws a circular plasmid map annotated with enzyme cut sites.

    Parameters
    ----------
    results : dict
        Enzyme data in the format:
        { enzyme_name: [recognition_seq, cut_notation, cut_count, [positions]] }
    plasmid_sequence : str
        The plasmid nucleotide sequence (used only for length; positions are
        taken from `results`).
    title : str
        Title printed in the centre of the plasmid map.
    """
    sys.stdout.write("\nGenerating plasmid map...\n")
    plasmid_len = len(plasmid_sequence)
    cutters = {k: v for k, v in results.items() if v[2] > 0}

    # ── layout constants ───────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_aspect("equal")
    ax.axis("off")

    R       = 1.0   # backbone radius
    tick_r  = 0.12  # radial length of tick marks
    lbl_gap = 0.22  # extra gap from backbone to label
    mc_gap  = 0.07  # extra gap for multi-cut enzymes
    lw_bb   = 6     # backbone line width

    # ── colour palette (one colour per enzyme) ─────────────────────────────────
    palette = plt.cm.tab20.colors
    enzyme_names = list(cutters.keys())
    colour_map = {name: palette[i % len(palette)] for i, name in enumerate(enzyme_names)}

    # ── backbone ───────────────────────────────────────────────────────────────
    backbone = plt.Circle((0, 0), R, fill=False, color="#444", linewidth=lw_bb)
    ax.add_patch(backbone)

    # ── tick mark + label for every cut site ───────────────────────────────────
    for enzyme, data in cutters.items():
        positions = data[3]
        cut_count = data[2]
        colour    = colour_map[enzyme]
        multi     = cut_count > 1

        for pos in positions:
            # angle: 0° at top (12 o'clock), clockwise
            angle_deg = 90 - (pos / plasmid_len) * 360
            angle_rad = np.radians(angle_deg)
            cos_a, sin_a = np.cos(angle_rad), np.sin(angle_rad)

            # inner / outer tick endpoints
            r_in  = R - tick_r / 2
            r_out = R + tick_r + (0 if not multi else mc_gap)
            ax.plot([cos_a * r_in,  cos_a * r_out],
                    [sin_a * r_in,  sin_a * r_out],
                    color=colour, linewidth=2.5, solid_capstyle="round")

            # label
            r_lbl = R + tick_r + lbl_gap + (0 if not multi else mc_gap)
            lx, ly = cos_a * r_lbl, sin_a * r_lbl

            # align label away from the centre
            ha = "left"  if cos_a >= 0 else "right"
            va = "bottom" if sin_a >= 0 else "top"

            label = f"{enzyme}\n({pos})"
            ax.text(lx, ly, label, fontsize=8.5, ha=ha, va=va,
                    color=colour, fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=colour, alpha=0.85, linewidth=0.8))

    # ── centre annotation ──────────────────────────────────────────────────────
    ax.text(0, 0.10, title,
            ha="center", va="center", fontsize=14,
            fontweight="bold", color="#222")
    ax.text(0, -0.10, f"{plasmid_len:,} bp",
            ha="center", va="center", fontsize=11, color="#555")
    ax.text(0, -0.28, f"{len(cutters)} unique cutters",
            ha="center", va="center", fontsize=9, color="#888")

    # ── position ticks every ~500 bp ───────────────────────────────────────────
    tick_interval = 500
    for bp in range(0, plasmid_len, tick_interval):
        a_deg = 90 - (bp / plasmid_len) * 360
        a_rad = np.radians(a_deg)
        r0, r1 = R - 0.04, R + 0.04
        ax.plot([np.cos(a_rad) * r0, np.cos(a_rad) * r1],
                [np.sin(a_rad) * r0, np.sin(a_rad) * r1],
                color="#aaa", linewidth=1)
        ax.text(np.cos(a_rad) * (r1 + 0.07),
                np.sin(a_rad) * (r1 + 0.07),
                str(bp), fontsize=7, ha="center", va="center", color="#aaa")

    # ── legend (non-cutters listed separately) ─────────────────────────────────
    non_cutters = [k for k, v in results.items() if v[2] == 0]
    handles = [mpatches.Patch(color=colour_map[n], label=f"{n}") for n in enzyme_names]
    leg = ax.legend(handles=handles, loc="lower left", fontsize=8,
                    title="Cutting enzymes", title_fontsize=8.5,
                    framealpha=0.9, ncol=1, bbox_to_anchor=(-0.22, 0))

    # small note listing non-cutters
    if non_cutters:
        nc_str = ", ".join(non_cutters)
        ax.text(0, -1.68, f"No cut: {nc_str}",
                ha="center", va="top", fontsize=7.5, color="#999",
                wrap=True)

    ax.set_xlim(-2.1, 2.1)
    ax.set_ylim(-1.9, 1.9)
    plt.tight_layout()
    plt.savefig("results/plasmid_map.png", dpi=150, bbox_inches="tight")
    plt.show()
    sys.stdout.write(f"\n--> Plasmid map saved to results/plasmid_map.png\n")
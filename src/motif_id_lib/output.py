
import sys, re, math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# --- Module-level constants ------------------
import random as _random
_random.seed(42)
_tab20 = list(plt.cm.tab20.colors) + list(plt.cm.tab20b.colors) + list(plt.cm.tab20c.colors)
_random.shuffle(_tab20)
PALETTE = _tab20

IUPAC = {
    'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
    'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
    'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]',
}

BASE_COLORS  = {'A': '#69db7c', 'T': '#ff6b6b', 'G': '#ffd43b', 'C': '#4dabf7'}
COMPLEMENT   = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

# sys.stdout colors:
HIGHLIGHT = "\033[36m" # teal
RESET = "\033[0m"

# ===============================
# Shared data-preparation mixin
# ===============================
class _EnzymeDataMixin:
    """
    Shared helpers for building colour maps, motif lookups, and cut-position
    tables from the results dictionary.  Subclasses must set:
        self.results, self.seq, self.seq_len, self.enzyme_names
    before calling any method here.
    """

    def _build_color_map(self) -> dict:
        """Assign a palette colour to each enzyme name."""
        return {n: PALETTE[i % len(PALETTE)] for i, n in enumerate(self.enzyme_names)}

    @staticmethod
    def _top_cut_offset(cut_notation: str) -> int:
        """Offset of the top-strand cut within the motif (^ position, _ stripped)."""
        return cut_notation.replace('_', '').index('^')

    @staticmethod
    def _bot_cut_offset(cut_notation: str) -> int:
        """Offset of the bottom-strand cut within the motif (_ position, ^ stripped)."""
        return cut_notation.replace('^', '').index('_')

    def _build_char_enzyme_map(self) -> dict:
        """
        Map every sequence index inside a recognition motif → enzyme name.
        IUPAC codes are expanded; last enzyme wins on overlap.
        """
        char_enzyme: dict[int, str] = {}
        for name, data in self.results.items():
            motif   = data[0].upper()
            pattern = ''.join(IUPAC.get(c, re.escape(c)) for c in motif)
            for m in re.finditer(pattern, self.seq):
                for i in range(m.start(), m.end()):
                    char_enzyme[i] = name
        return char_enzyme

    def _build_top_cut_positions(self) -> dict:
        """
        Map absolute top-strand cut index → [enzyme names].
        abs_cut = site_pos + top_cut_offset
        """
        positions: dict[int, list[str]] = {}
        for name, data in self.results.items():
            _, notation, cut_count, sites = data
            if cut_count == 0:
                continue
            offset = self._top_cut_offset(notation)
            for site in sites:
                positions.setdefault(site + offset, []).append(name)
        return positions

    def _build_bot_cut_positions(self) -> dict:
        """
        Map absolute bottom-strand cut index → [enzyme names].
        abs_cut = site_pos + bot_cut_offset
        """
        positions: dict[int, list[str]] = {}
        for name, data in self.results.items():
            _, notation, cut_count, sites = data
            if cut_count == 0:
                continue
            offset = self._bot_cut_offset(notation)
            for site in sites:
                positions.setdefault(site + offset, []).append(name)
        return positions

    def _build_legend_handles(self, color_map: dict) -> list:
        """Coloured patch handles for the enzyme legend."""
        handles = []
        for name, data in self.results.items():
            if data[2] == 0:
                continue
            pos_str = ', '.join(str(p) for p in data[3])
            label   = f"{name} [{data[0]}] @ {pos_str}"
            handles.append(mpatches.Patch(
                facecolor=color_map[name], edgecolor=color_map[name],
                alpha=0.8, label=label,
            ))
        return handles


# ==================
# CircularMap
# ==================

class CircularMap(_EnzymeDataMixin):
    """
    Circular plasmid map with restriction-enzyme cut sites.

    Parameters
    ----------
    results : dict
        { enzyme_name: [recognition_seq, cut_notation, cut_count, [positions]] }
    plasmid_sequence : str
        Full nucleotide sequence.
    title : str
        Text printed at the centre of the map.
    figsize : tuple
        Figure dimensions in inches.
    dpi : int
        Output resolution.
    output_path : str
        Destination file for the saved figure.
    """

    # --- Geometry constants ------------------------------
    RADIUS    : float = 1.6    # backbone circle radius
    TICK_R    : float = 0.12   # radial length of cut-site tick marks
    LBL_GAP   : float = 0.22   # extra radial gap from tick to label
    MC_GAP    : float = 0.07   # extra radial gap for multi-cut enzymes
    BB_LW     : int   = 3      # backbone linewidth
    POS_TICK_INTERVAL: int = 500  # bp between ruler ticks on the backbone

    def __init__(
        self,
        results: dict,
        plasmid_sequence: str,
        title: str = "Plasmid",
        figsize: tuple = (10, 10),
        dpi: int = 150,
        output_path: str = "results/circular_map.png",
    ) -> None:
        self.results      = results
        self.seq          = plasmid_sequence.upper().replace(' ', '').replace('\n', '')
        self.seq_len      = len(self.seq)
        self.title        = title
        self.figsize      = figsize
        self.DPI          = dpi
        self.OUTPUT_PATH  = output_path

        self.enzyme_names = list(results.keys())
        self.color_map    = self._build_color_map()
        self.cutters      = {k: v for k, v in results.items() if v[2] > 0}

    # --- Figure setup -------------------------------------------------

    def _make_figure(self) -> tuple:
        """Return (fig, ax) with equal-aspect, axis-off settings."""
        fig, ax = plt.subplots(figsize=self.figsize)
        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_xlim(-3.2, 3.2)
        ax.set_ylim(-3.0, 3.0)
        fig.patch.set_facecolor('white')
        return fig, ax

    # --- Drawing helpers ---------------------------------------------

    def _draw_backbone(self, ax) -> None:
        """Draw the circular backbone."""
        ax.add_patch(plt.Circle(
            (0, 0), self.RADIUS,
            fill=False, color='#444', linewidth=self.BB_LW,
        ))

    def _draw_position_ticks(self, ax) -> None:
        """Draw evenly spaced bp-position ticks and labels around the backbone."""
        R = self.RADIUS
        for bp in range(0, self.seq_len, self.POS_TICK_INTERVAL):
            a = np.radians(90 - (bp / self.seq_len) * 360)
            r0, r1 = R - 0.04, R + 0.04
            ax.plot([np.cos(a) * r0, np.cos(a) * r1],
                    [np.sin(a) * r0, np.sin(a) * r1],
                    color='#aaa', linewidth=1)
            ax.text(np.cos(a) * (r1 + 0.07), np.sin(a) * (r1 + 0.07),
                    str(bp), fontsize=7, ha='center', va='center', color='#aaa')

    # --- Stacking constants --------------------------------------------
    MIN_ANG_GAP  : float = 0.10   # minimum angular gap (radians) before stacking
    STACK_STEP   : float = 0.26   # radial step per stack level
    STEM_COLOR   : str   = '#aaaaaa'  # colour of the leader line from tick to label

    def _compute_label_radii(self) -> list:
        """
        Return a sorted list of (angle_rad, enzyme, position, colour, r_tick, r_label)
        for every cut event.

        Labels that are angularly too close to an already-placed label at the
        same radial band are pushed outward by STACK_STEP until they fit,
        reproducing the stacked-flag style seen in tools like SnapGene/NEB.
        """
        R = self.RADIUS
        events = []
        for enzyme, data in self.cutters.items():
            col   = self.color_map[enzyme]
            multi = data[2] > 1
            for pos in data[3]:
                ang    = np.radians(90 - (pos / self.seq_len) * 360)
                r_tick = R + self.TICK_R + (self.MC_GAP if multi else 0)
                events.append((ang, enzyme, pos, col, r_tick))

        events.sort(key=lambda e: e[0])

        placed = []   # (angle, r_label) of already-positioned labels
        result = []
        for ang, enzyme, pos, col, r_tick in events:
            r_lbl = r_tick + self.LBL_GAP
            while True:
                conflict = any(
                    abs(p_r - r_lbl) < 0.05 and
                    min(abs(ang - p_a), 2*math.pi - abs(ang - p_a)) < self.MIN_ANG_GAP
                    for p_a, p_r in placed
                )
                if not conflict:
                    break
                r_lbl += self.STACK_STEP
            placed.append((ang, r_lbl))
            result.append((ang, enzyme, pos, col, r_tick, r_lbl))

        return result

    def _draw_cut_sites(self, ax) -> None:
        """
        Draw a radial tick, a leader line, and a stacked label for every cut
        site.  Closely spaced sites are pushed to increasing radii so labels
        never overlap, matching the style of the reference image.
        """
        R       = self.RADIUS
        events  = self._compute_label_radii()

        for ang, enzyme, pos, col, r_tick, r_lbl in events:
            cos_a = np.cos(ang)
            sin_a = np.sin(ang)

            # Radial tick on the backbone
            ax.plot([cos_a * (R - self.TICK_R / 2), cos_a * r_tick],
                    [sin_a * (R - self.TICK_R / 2), sin_a * r_tick],
                    color=col, linewidth=2.5, solid_capstyle='round', zorder=4)

            # Leader line from tick tip to label anchor
            ax.plot([cos_a * r_tick, cos_a * r_lbl],
                    [sin_a * r_tick, sin_a * r_lbl],
                    color=self.STEM_COLOR, linewidth=0.8,
                    linestyle='--', alpha=0.6, zorder=3)

            # Label
            ha = 'left' if cos_a >= 0 else 'right'
            va = 'bottom' if sin_a >= 0 else 'top'
            ax.text(
                cos_a * r_lbl, sin_a * r_lbl,
                f"{enzyme}\n({pos})",
                fontsize=7.5, ha=ha, va=va,
                color=col, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.15', fc='white', ec=col, alpha=0.90, linewidth=0.7),
                zorder=5,
            )

    def _draw_centre_text(self, ax) -> None:
        """Draw title, bp count, and cutter count at the centre of the circle."""
        ax.text(0,  0.10, self.title,
                ha='center', va='center', fontsize=14, fontweight='bold', color='#222')
        ax.text(0, -0.10, f"{self.seq_len:,} bp",
                ha='center', va='center', fontsize=11, color='#555')
        ax.text(0, -0.28, f"{len(self.cutters)} unique cutters",
                ha='center', va='center', fontsize=9, color='#888')

    def _draw_legend(self, ax) -> None:
        """Draw the enzyme colour legend and a non-cutter note."""
        handles = [
            mpatches.Patch(color=self.color_map[n], label=n)
            for n in self.cutters
        ]
        if handles:
            ax.legend(handles=handles, loc='lower left', fontsize=8, title='Cutting enzymes', title_fontsize=8.5, framealpha=0.9, ncol=1, bbox_to_anchor=(-0.22, 0))

        non_cutters = [k for k, v in self.results.items() if v[2] == 0]
        if non_cutters:
            ax.get_figure().text(
                0.5, 0.01,
                f"No cut: {', '.join(non_cutters)}",
                ha='center', va='bottom', fontsize=8, color='#999',
                wrap=True,
            )

    # --- Public entry point ---------------------------------------------

    def render(self) -> plt.Figure:
        """Build, display, and save the circular map."""
        sys.stdout.write('\n\t2. Generating circular plasmid map...\n')
        fig, ax = self._make_figure()
        self._draw_backbone(ax)
        self._draw_position_ticks(ax)
        self._draw_cut_sites(ax)
        self._draw_centre_text(ax)
        self._draw_legend(ax)
        plt.tight_layout()
        if self.OUTPUT_PATH:
            fig.savefig(self.OUTPUT_PATH, dpi=self.DPI, bbox_inches='tight')
            sys.stdout.write(f'\t--> Circular map saved to {HIGHLIGHT}{self.OUTPUT_PATH}{RESET}\n')
        return fig


# =======================================================
# LinearMap  (single-stranded sequence viewer)
# =======================================================

class LinearMap(_EnzymeDataMixin):
    """
    Single-stranded sequence viewer.  Wraps the sequence into lines of
    CHARS_PER_LINE bases, colours each base by identity, highlights recognition
    motifs, and marks top-strand cut sites with a vertical tick in the inter-base
    gap plus a labelled arrow.

    Parameters
    ----------
    results : dict
        { enzyme_name: [recognition_seq, cut_notation, cut_count, [site_positions]] }
    plasmid_sequence : str
        Full nucleotide sequence (top strand, 5'→3').
    title : str
        Header text.
    figwidth : float
        Figure width in inches.
    chars_per_line : int
        Bases per wrapped line.
    dpi : int
        Output resolution.
    output_path : str
        Destination file.
    """

    # --- Geometry constants --------------------------------------------
    SEQ_LINE_H : float = 0.62   # inches per wrapped line
    LEFT_MAR   : float = 0.065  # x-fraction reserved for position labels
    HEADER_H   : float = 1.4    # fixed height of the header panel in inches

    def __init__(
        self,
        results: dict,
        plasmid_sequence: str,
        title: str = "Plasmid",
        figwidth: float = 14.0,
        chars_per_line: int = 80,
        dpi: int = 150,
        output_path: str = "results/linear_map.png",
    ) -> None:
        self.results        = results
        self.seq            = plasmid_sequence.upper().replace(' ', '').replace('\n', '')
        self.seq_len        = len(self.seq)
        self.title          = title
        self.FIGWIDTH       = figwidth
        self.CHARS_PER_LINE = chars_per_line
        self.DPI            = dpi
        self.OUTPUT_PATH    = output_path

        self.enzyme_names  = list(results.keys())
        self.color_map     = self._build_color_map()
        self.char_enzyme   = self._build_char_enzyme_map()
        self.cut_positions = self._build_top_cut_positions()

    # --- Figure setup -------------------------------------------------

    def _cuts_per_line(self) -> dict:
        """Return {line_idx: max_stack_depth} for lines containing cut sites."""
        result: dict[int, int] = {}
        for abs_cut, enzymes in self.cut_positions.items():
            li = abs_cut // self.CHARS_PER_LINE
            result[li] = max(result.get(li, 0), len(enzymes))
        return result

    def _make_figure(self) -> tuple:
        """
        Return (fig, ax_header, ax_seq).
        Header panel is fixed-height; sequence panel scales with content.
        """
        num_lines = math.ceil(self.seq_len / self.CHARS_PER_LINE)
        extra     = sum(min(d, 1) * 0.30 for d in self._cuts_per_line().values())
        seq_h     = num_lines * self.SEQ_LINE_H + extra
        fig_h     = self.HEADER_H + seq_h

        fig = plt.figure(figsize=(self.FIGWIDTH, fig_h))
        fig.patch.set_facecolor('#1e2228')

        gs = fig.add_gridspec(2, 1, height_ratios=[self.HEADER_H, seq_h], hspace=0.04)
        ax_hdr = fig.add_subplot(gs[0])
        ax_seq = fig.add_subplot(gs[1])
        for ax in (ax_hdr, ax_seq):
            ax.set_facecolor('#1e2228')
        return fig, ax_hdr, ax_seq

    # --- Drawing helpers ---------------------------------------------

    def _draw_header(self, ax) -> None:
        """Title, bp count, and enzyme legend in the dedicated header panel."""
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        ax.text(0.0, 0.95, self.title,
                ha='left', va='top', fontsize=14, fontweight='700',
                color='#00d4ff', fontfamily='monospace',
                transform=ax.transAxes)
        ax.text(1.0, 0.95, f"{self.seq_len} bp",
                ha='center', va='top', fontsize=12,
                color='#00d4ff', fontfamily='monospace',
                transform=ax.transAxes)
        handles = self._build_legend_handles(self.color_map)
        if handles:
            leg = ax.legend(
                handles=handles,
                loc='upper left', bbox_to_anchor=(0.0, 0.80),
                bbox_transform=ax.transAxes,
                ncol=min(4, len(handles)),
                fontsize=9, frameon=True, framealpha=0.12,
                edgecolor='#2a4a6b', facecolor='#252b33', labelcolor='#cdd9e5',
                handlelength=1.0, handleheight=0.8, columnspacing=1.2,
            )
            for t in leg.get_texts():
                t.set_fontfamily('monospace')

    def _setup_seq_axes(self, ax, num_lines: int) -> None:
        """Configure coordinate space for the sequence panel."""
        ax.set_xlim(0, 1)
        ax.set_ylim(-num_lines - 0.5, 0.5)
        ax.axis('off')

    def _draw_all_lines(self, ax) -> None:
        """Render every wrapped sequence line."""
        num_lines = math.ceil(self.seq_len / self.CHARS_PER_LINE)
        char_w    = (1.0 - self.LEFT_MAR - 0.005) / self.CHARS_PER_LINE
        for li in range(num_lines):
            self._draw_sequence_line(ax, li, char_w)

    def _draw_sequence_line(self, ax, line_idx: int, char_w: float) -> None:
        """One wrapped line: position label, bases, and cut annotations."""
        start = line_idx * self.CHARS_PER_LINE
        end   = min(start + self.CHARS_PER_LINE, self.seq_len)
        y     = -(line_idx + 0.5)

        ax.text(self.LEFT_MAR - 0.004, y, f"{start + 1:>6}",
                ha='right', va='center', fontsize=7,
                color='#3a6a9b', fontfamily='monospace')

        for ci, base in enumerate(self.seq[start:end]):
            abs_idx = start + ci
            x_left  = self.LEFT_MAR + ci * char_w
            x_mid   = x_left + char_w / 2
            self._draw_base(ax, base, abs_idx, x_left, x_mid, y, char_w)
            if abs_idx in self.cut_positions:
                self._draw_cut_annotation(ax, abs_idx, x_left, y, char_w, self.cut_positions, above=True)

    def _draw_base(self, ax, base: str, abs_idx: int, x_left: float, x_mid: float, y: float, char_w: float) -> None:
        """Draw one nucleotide with optional motif highlight."""
        hit = self.char_enzyme.get(abs_idx)
        if hit:
            ax.add_patch(mpatches.FancyBboxPatch(
                (x_left, y - 0.40), char_w, 0.80,
                boxstyle='square,pad=0.0', linewidth=0,
                facecolor=self.color_map[hit], alpha=0.35, zorder=2,
            ))
            txt_col, bold = '#ffffff', True
        else:
            txt_col, bold = BASE_COLORS.get(base, '#cdd9e5'), False

        ax.text(x_mid, y, base, ha='center', va='center', fontsize=7.5,
                color=txt_col, fontfamily='monospace',
                fontweight='bold' if bold else 'normal', zorder=3)

    def _draw_cut_annotation(
        self, ax, abs_idx: int, x_left: float, y: float, char_w: float,
        cut_positions: dict, above: bool = True,
    ) -> None:
        """
        Vertical tick in the inter-base gap at x_left, with a stacked arrow
        and label.  above=True sends arrows upward (top strand); False sends
        them downward (bottom strand).
        """
        TICK_HALF = 0.44
        direction = 1 if above else -1

        for stack_i, enz_name in enumerate(cut_positions[abs_idx]):
            col      = self.color_map[enz_name]
            tick_top = y + TICK_HALF
            tick_bot = y - TICK_HALF
            tip_y    = tick_top if above else tick_bot
            origin_y = tip_y + direction * (0.20 + stack_i * 0.30)

            ax.plot([x_left, x_left], [tick_bot, tick_top],
                    color=col, lw=2.0, zorder=6, solid_capstyle='round')
            ax.annotate('',
                        xy=(x_left, tip_y), xytext=(x_left, origin_y),
                        arrowprops=dict(arrowstyle='->', color=col, lw=1.1),
                        zorder=5)
            ax.text(x_left + char_w * 0.3, origin_y,
                    f"{enz_name} ↓{abs_idx}",
                    ha='left', va='center', fontsize=6.5,
                    color=col, fontfamily='monospace',
                    fontweight='600', zorder=5)

    # --- Public entry point --------------------------------------------

    def render(self) -> plt.Figure:
        """Build, display, and save the single-stranded linear map."""
        sys.stdout.write('\n\t3. Generating single-strand linear plasmid map...\n')
        num_lines = math.ceil(self.seq_len / self.CHARS_PER_LINE)
        fig, ax_hdr, ax_seq = self._make_figure()
        self._draw_header(ax_hdr)
        self._setup_seq_axes(ax_seq, num_lines)
        self._draw_all_lines(ax_seq)
        fig.subplots_adjust(left=0.03, right=0.97, top=0.97, bottom=0.03, hspace=0.04)
        if self.OUTPUT_PATH:
            fig.savefig(self.OUTPUT_PATH, dpi=self.DPI, bbox_inches='tight',
                        facecolor=fig.get_facecolor())
            sys.stdout.write(f'\t--> Linear map saved to {HIGHLIGHT}{self.OUTPUT_PATH}{RESET}\n')
        return fig


# =======================================================
# DoubleStrandedMap
# =======================================================

class DoubleStrandedMap(LinearMap):
    """
    Double-stranded sequence viewer.  Inherits all single-strand drawing from
    LinearMap and adds a bottom (complement, 3'→5') strand beneath each top
    strand line.

    Layout per wrapped block
    ------------------------
        [annotation arrows — top strand cuts — pointing UP]
        5'  A T G C A T …  3'   ← top strand, coloured by base
            | | | | | |         ← pairing tick marks
        3'  T A C G T A …  5'   ← bottom strand, coloured by base (complement)
        [annotation arrows — bottom strand cuts — pointing DOWN]

    Top-strand cut positions use the ^ offset from the notation.
    Bottom-strand cut positions use the _ offset from the notation.
    Both arrows are drawn in the inter-base gap (x_left of the cut index).
    """

    # Geometry for within-block layout (strand rows + pairing marks + arrows)
    BOT_STRAND_H  : float = 0.55   # inches of extra figure height per block
    PAIR_TICK_H   : float = 0.15   # half-height of base-pairing tick marks
    TOP_ARROW_H   : float = 0.44   # data-units above top strand for upward arrows
    BOT_ARROW_H   : float = 0.44   # data-units below bottom strand for downward arrows
    STRAND_GAP    : float = 0.70   # data-unit gap between top and bottom strand rows
    BLOCK_GAP     : float = 1.20   # extra data-units of whitespace *between* blocks

    def __init__(
        self,
        results: dict,
        plasmid_sequence: str,
        title: str = "Plasmid",
        figwidth: float = 14.0,
        chars_per_line: int = 80,
        dpi: int = 150,
        output_path: str = "results/ds_linear_map.png",
    ) -> None:
        super().__init__(results, plasmid_sequence, title, figwidth, chars_per_line, dpi, output_path)

        # Complement strand (3'→5', so same left-to-right order as top strand)
        self.comp_seq       = ''.join(COMPLEMENT.get(b, 'N') for b in self.seq)
        # Bottom-strand cut positions (keyed by the same abs index as top strand)
        self.bot_cut_positions = self._build_bot_cut_positions()

    # --- Override geometry ---------------------------------------------

    def _seq_line_height(self) -> float:
        """Height in inches for one double-stranded block."""
        return self.SEQ_LINE_H + self.BOT_STRAND_H

    def _make_figure(self) -> tuple:
        """Taller figure to accommodate both strands, pairing marks, and arrows."""
        num_lines   = math.ceil(self.seq_len / self.CHARS_PER_LINE)
        # Count lines that have any cut site on either strand
        top_cuts    = self._cuts_per_line()
        bot_line    = {pos // self.CHARS_PER_LINE: len(enzs)
                       for pos, enzs in self.bot_cut_positions.items()}
        all_lines   = set(top_cuts) | set(bot_line)
        extra       = len(all_lines) * 0.30
        seq_h       = num_lines * (self._seq_line_height() + self.BLOCK_GAP) + extra
        fig_h     = self.HEADER_H + seq_h

        fig = plt.figure(figsize=(self.FIGWIDTH, fig_h))
        fig.patch.set_facecolor('#1e2228')
        gs = fig.add_gridspec(2, 1, height_ratios=[self.HEADER_H, seq_h], hspace=0.04)
        ax_hdr = fig.add_subplot(gs[0])
        ax_seq = fig.add_subplot(gs[1])
        for ax in (ax_hdr, ax_seq):
            ax.set_facecolor('#1e2228')
        return fig, ax_hdr, ax_seq

    def _setup_seq_axes(self, ax, num_lines: int) -> None:
        """y-range sized to fit all blocks including inter-block gaps."""
        line_h  = self._seq_line_height()
        # Total vertical span: num_lines blocks each of (line_h + BLOCK_GAP),
        # plus the final block's internal height and arrow clearance.
        y_min   = -(num_lines * (line_h + self.BLOCK_GAP)
                    + self.STRAND_GAP + self.BOT_ARROW_H + 0.5)
        ax.set_xlim(0, 1)
        ax.set_ylim(y_min, 0.5)
        ax.axis('off')

    # --- Override line drawing -----------------------------------------

    def _draw_all_lines(self, ax) -> None:
        """Render every double-stranded block."""
        num_lines = math.ceil(self.seq_len / self.CHARS_PER_LINE)
        char_w    = (1.0 - self.LEFT_MAR - 0.005) / self.CHARS_PER_LINE
        for li in range(num_lines):
            self._draw_ds_block(ax, li, char_w)

    def _draw_ds_block(self, ax, line_idx: int, char_w: float) -> None:
        """
        Draw one double-stranded block:
            y_top  = centre of the top-strand row
            y_bot  = centre of the bottom-strand row (below top)
        """
        line_h  = self._seq_line_height()
        # Each block starts at its own base y, shifted by the cumulative BLOCK_GAP
        base_y  = -(line_idx * (line_h + self.BLOCK_GAP) + self.TOP_ARROW_H)
        y_top   = base_y
        y_bot   = y_top - self.STRAND_GAP

        start = line_idx * self.CHARS_PER_LINE
        end   = min(start + self.CHARS_PER_LINE, self.seq_len)

        # Position label — centred between the two strand rows so it doesn't
        # collide with the 5' strand label at y_top
        y_mid = (y_top + y_bot) / 2
        ax.text(self.LEFT_MAR - 0.004, y_mid, f"{start + 1:>6}",
                ha='right', va='center', fontsize=7,
                color='#3a6a9b', fontfamily='monospace')

        # Strand labels
        ax.text(self.LEFT_MAR - 0.004, y_top,   "5'",
                ha='right', va='bottom', fontsize=6, color='#7a96aa',
                fontfamily='monospace')
        ax.text(self.LEFT_MAR - 0.004, y_bot,   "3'",
                ha='right', va='top',    fontsize=6, color='#7a96aa',
                fontfamily='monospace')
        ax.text(1.0,                   y_top,   "3'",
                ha='left',  va='bottom', fontsize=6, color='#556677',
                fontfamily='monospace')
        ax.text(1.0,                   y_bot,   "5'",
                ha='left',  va='top',    fontsize=6, color='#556677',
                fontfamily='monospace')

        for ci in range(end - start):
            abs_idx    = start + ci
            top_base   = self.seq[abs_idx]
            bot_base   = self.comp_seq[abs_idx]
            x_left     = self.LEFT_MAR + ci * char_w
            x_mid      = x_left + char_w / 2

            # Top strand base
            self._draw_base(ax, top_base, abs_idx, x_left, x_mid, y_top, char_w)
            # Bottom strand base (complement; use same motif highlight if applicable)
            self._draw_base(ax, bot_base, abs_idx, x_left, x_mid, y_bot, char_w)
            # Base-pairing tick
            self._draw_pair_tick(ax, x_mid, y_top, y_bot)

            # Top-strand cut annotation (arrows point UP)
            if abs_idx in self.cut_positions:
                self._draw_cut_annotation(ax, abs_idx, x_left, y_top, char_w, self.cut_positions, above=True)
            # Bottom-strand cut annotation (arrows point DOWN)
            if abs_idx in self.bot_cut_positions:
                self._draw_cut_annotation(ax, abs_idx, x_left, y_bot, char_w, self.bot_cut_positions, above=False)

    def _draw_pair_tick(self, ax, x_mid: float, y_top: float, y_bot: float) -> None:
        """Draw a short vertical line between a paired base pair."""
        gap    = y_bot - y_top          # negative value
        centre = y_top + gap / 2
        ax.plot([x_mid, x_mid],
                [centre - self.PAIR_TICK_H, centre + self.PAIR_TICK_H],
                color='#334455', lw=0.6, zorder=1)

    # --- Public entry point --------------------------------------------

    def render(self) -> plt.Figure:
        """Build, display, and save the double-stranded linear map."""
        sys.stdout.write('\n\t4. Generating double-stranded linear plasmid map...\n')
        num_lines = math.ceil(self.seq_len / self.CHARS_PER_LINE)
        fig, ax_hdr, ax_seq = self._make_figure()
        self._draw_header(ax_hdr)
        self._setup_seq_axes(ax_seq, num_lines)
        self._draw_all_lines(ax_seq)
        fig.subplots_adjust(left=0.03, right=0.97, top=0.97, bottom=0.03, hspace=0.04)
        if self.OUTPUT_PATH:
            fig.savefig(self.OUTPUT_PATH, dpi=self.DPI, bbox_inches='tight',
                        facecolor=fig.get_facecolor())
            sys.stdout.write(f'\t--> Double-stranded map saved to {HIGHLIGHT}{self.OUTPUT_PATH}{RESET}\n')
        return fig

class PlasmidMap:
    """
    Convenience wrapper that exposes circular and linear rendering through a
    single object.  All heavy lifting is delegated to CircularMap, LinearMap,
    and DoubleStrandedMap.
    """

    def __init__(self, results: dict, plasmid_sequence: str, title: str = "Plasmid name") -> None:
        self.results          = results
        self.plasmid_sequence = plasmid_sequence
        self.title            = title

    def annotate_circular(
        self,
        figsize: tuple = (10, 10),
        dpi: int = 150,
        output_path: str = "results/circular_map.png",
    ) -> plt.Figure:
        """Render a circular plasmid map."""
        short_title = self.title.split(",")[0]
        return CircularMap(
            self.results, self.plasmid_sequence, short_title,
            figsize=figsize,
            dpi=dpi,
            output_path=output_path,
        ).render()

    def annotate_linear(
        self,
        figwidth: float = 14.0,
        chars_per_line: int = 80,
        dpi: int = 150,
        output_path: str = "results/linear_map.png",
    ) -> plt.Figure:
        """Render a single-stranded linear sequence map."""
        return LinearMap(
            self.results, self.plasmid_sequence, self.title,
            figwidth=figwidth,
            chars_per_line=chars_per_line,
            dpi=dpi,
            output_path=output_path,
        ).render()

    def annotate_double_stranded(
        self,
        figwidth: float = 14.0,
        chars_per_line: int = 80,
        dpi: int = 150,
        output_path: str = "results/ds_linear_map.png",
    ) -> plt.Figure:
        """Render a double-stranded linear sequence map."""
        return DoubleStrandedMap(
            self.results, self.plasmid_sequence, self.title,
            figwidth=figwidth,
            chars_per_line=chars_per_line,
            dpi=dpi,
            output_path=output_path,
        ).render()
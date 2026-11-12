"""Render an oriented CPS gene map from reconstructed locus pieces."""

from math import floor, log10

from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib import pyplot as plt

EXPECTED = '#0072B2'
OTHER = '#D55E00'
INK = '#182632'
MUTED = '#526373'


def scale_length(span):
    target = span / 5
    magnitude = 10 ** floor(log10(target))
    return max(1, max(step for step in (1, 2, 5) if step * magnitude <= target) * magnitude)


def label_genes(axis, labels):
    """Fit labels on three rows, keeping named genes ahead of fallback IDs."""
    renderer = axis.figure.canvas.get_renderer()
    occupied = [[], [], []]
    for _, x, text, color in sorted(labels):
        annotation = axis.text(x, 0.42, text, ha='center', va='bottom', fontsize=8, color=color)
        bounds = annotation.get_window_extent(renderer)
        left, right = bounds.x0 - 4, bounds.x1 + 4
        for row, intervals in enumerate(occupied):
            if (axis.bbox.x0 <= left and right <= axis.bbox.x1
                    and all(right <= start or left >= end for start, end in intervals)):
                intervals.append((left, right))
                y = 0.42 + row * 0.26
                annotation.set_y(y)
                if row:
                    axis.vlines(x, 0.35, y - 0.04, color='#B8C4CD', linewidth=0.6, zorder=1)
                break
        else:
            annotation.remove()


def write_figure(path, result, sample, candidate):
    path.parent.mkdir(exist_ok=True)
    span = sum(piece.end - piece.start for piece in result.pieces)
    gap = max(1, span * min(0.032, 0.28 / len(result.pieces)))
    extent = span + gap * (len(result.pieces) - 1)

    figure, axis = plt.subplots(figsize=(12, 2.6), facecolor='white')
    try:
        figure.subplots_adjust(left=0.047, right=0.953, bottom=0.13, top=0.74)
        figure.text(0.047, 0.88, sample, color=INK, fontsize=15, weight='bold')
        figure.text(0.047, 0.78, f'Best serotype: {candidate}', color=MUTED, fontsize=9)
        legend = [Patch(facecolor=EXPECTED, label='Expected gene'),
                  Patch(facecolor=OTHER, label='Other-locus gene hit')]
        if len(result.pieces) > 1:
            legend.append(Line2D([], [], marker=r'$//$', markersize=9, linestyle='none',
                                 color=MUTED, label='Fragment break'))
        figure.legend(handles=legend, loc='upper right', bbox_to_anchor=(0.96, 0.95),
                      ncol=len(legend), frameon=False, fontsize=8, labelcolor=INK,
                      handlelength=1.25, columnspacing=1.2)

        padding = extent * 0.035
        axis.set(xlim=(-padding, extent + padding), ylim=(-1.15, 1.24))
        axis.set_axis_off()
        labels = []
        offset = 0
        for index, piece in enumerate(result.pieces, 1):
            length = piece.end - piece.start
            stop = offset + length
            axis.hlines(0, offset, stop, color='#98A8B4', linewidth=1.2)
            for hits, color in ((piece.expected, EXPECTED), (piece.unexpected, OTHER)):
                for hit in hits:
                    if piece.strand == '+':
                        start, end = hit.start - piece.start, hit.end - piece.start
                    else:
                        start, end = piece.end - hit.end, piece.end - hit.start
                    if end <= start:
                        continue
                    tail, head = (start, end) if hit.alignment.strand == piece.strand else (end, start)
                    axis.arrow(offset + tail, 0, head - tail, 0, width=0.28, head_width=0.53,
                               head_length=min((end - start) * 0.40, span * 0.010),
                               length_includes_head=True, facecolor=color, edgecolor='white',
                               linewidth=0.45, zorder=3)
                    short_id = hit.gene.label.rsplit('_', 1)[-1]
                    label = hit.gene.symbol or short_id
                    labels.append((not hit.gene.symbol, offset + (start + end) / 2, label, color))
            center = (offset + stop) / 2
            if length > extent * 0.045:
                axis.text(center, -0.39, f'{index}  {piece.contig}' if length > extent * 0.17 else str(index),
                          ha='center', va='top', fontsize=8, color=MUTED, parse_math=False)
            if index < len(result.pieces):
                axis.text(stop + gap / 2, 0, '//', ha='center', va='center', color=MUTED, fontsize=9)
            offset = stop + gap

        bar = scale_length(span)
        units = f'{bar / 1000:g} kb' if bar >= 1000 else f'{bar:g} bp'
        bar_y = -0.87
        axis.hlines(bar_y, 0, bar, color=INK, linewidth=1)
        axis.vlines([0, bar], bar_y - 0.05, bar_y + 0.05, color=INK, linewidth=1)
        axis.text(bar / 2, bar_y - 0.09, units, ha='center', va='top', fontsize=8, color=INK)
        label_genes(axis, labels)
        figure.savefig(path, dpi=300)
        with plt.rc_context({'svg.fonttype': 'none'}):
            figure.savefig(path.with_suffix('.svg'))
    finally:
        plt.close(figure)

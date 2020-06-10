import os
import numpy as np
import pandas as pd
import matplotlib.colors as m_colors
import matplotlib.colorbar as m_colorbar
from matplotlib import pyplot as plt, ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
from package.helpers import ChromosomePosition

sns.set(font_scale=1.2, style="ticks", font="lato", palette=('#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2',
                                                             '#D55E00', '#CC79A7'))
# sns.set_style({"xtick.direction": "in", "ytick.direction": "in"})
plt.rcParams['font.weight'] = "medium"
plt.rcParams['axes.labelweight'] = 'medium'
plt.rcParams['figure.titleweight'] = 'medium'
plt.rcParams['axes.titleweight'] = 'medium'
plt.rcParams['figure.figsize'] = 14, 3
plt.rcParams["legend.framealpha"] = 1
plt.rcParams['axes.xmargin'] = 0
plt.rcParams['axes.ymargin'] = 0
plt.rcParams["legend.framealpha"] = 1

# snps = pd.read_table(snps_name, header=None)
# snps.columns = ['chr', 'pos', 'ID', 'ref', 'alt', 'ref_c', 'alt_c']
# snps['AD'] = snps[['ref_c', 'alt_c']].max(axis=1) / snps[['ref_c', 'alt_c']].min(axis=1)
# snps['cov'] = snps['ref_c'] + snps['alt_c']


# def init_visualization(snps, BAD_file):
#     BAD_table = pd.read_table(BAD_file)
#     file_name = os.path.splitext(os.path.basename(BAD_file))[0]
#     out_path = os.path.dirname(BAD_file)
#
#     snps['AD'] = snps[['ref_c', 'alt_c']].max(axis=1) / snps[['ref_c', 'alt_c']].min(axis=1)
#     snps['cov'] = snps['ref_c'] + snps['alt_c']
#
#     for chromosome in snps['chr'].unique():
#         visualize_chromosome(os.path.join(out_path, '{}_{}.svg'.format(file_name, chromosome)),
#                              chromosome, snps[snps['chr'] == chromosome],
#                              BAD_table[BAD_table['#chr'] == chromosome])


def init_from_snps_collection(snps_collection, BAD_file):
    BAD_table = pd.read_table(BAD_file)
    file_name = os.path.splitext(os.path.basename(BAD_file))[0]
    out_path = os.path.join(os.path.dirname(BAD_file), '{}_visualization'.format(file_name))
    if not os.path.isdir(out_path):
        os.mkdir(out_path)

    column_names = ['pos', 'ref_c', 'alt_c']
    for chromosome in snps_collection.keys():
        print('Visualizing {}'.format(chromosome))
        snps = pd.DataFrame(dict(zip(column_names, zip(*snps_collection[chromosome]))))
        snps['AD'] = snps[['ref_c', 'alt_c']].max(axis=1) / snps[['ref_c', 'alt_c']].min(axis=1)
        snps['cov'] = snps['ref_c'] + snps['alt_c']
        visualize_chromosome(os.path.join(out_path, '{}_{}.svg'.format(file_name, chromosome)),
                             chromosome, snps,
                             BAD_table[BAD_table['#chr'] == chromosome])


def visualize_chromosome(out_path, chromosome, snps, BAD_segments):
    fig, ax = plt.subplots()
    fig.tight_layout(rect=[0, 0.01, 0.95, 1])
    plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useMathText=True))

    BAD_color = '#0072B2CC'
    BAD_lw = 10
    y_min = 0.8
    y_max = 6
    delta_y = 0.05

    snps['AD'] = snps['AD'].apply(lambda y: y_max - delta_y if y > y_max else y)

    bar_positions = []
    bar_widths = []
    bar_colors = []
    gap = 1 / 500 * ChromosomePosition.chromosomes[chromosome]

    BADs = []

    borders_to_draw = []
    segmentation_borders = []
    last_end = 1
    for index, (pl_chr, start, end, BAD, *values) in BAD_segments.iterrows():
        if start != last_end:
            if last_end == 1:
                borders_to_draw += [start - gap]
                segmentation_borders += [start - gap]
            else:
                borders_to_draw += [last_end + gap, start - gap]
                segmentation_borders += [last_end + gap, start - gap]
            bar_colors.append('#AAAAAA')
            BADs.append(None)
        else:
            if last_end != 1:
                segmentation_borders += [last_end]
        last_end = end
        bar_colors.append('C2')
        BADs.append(BAD)
    if last_end != ChromosomePosition.chromosomes[chromosome] + 1:
        borders_to_draw += [last_end + gap]
        segmentation_borders += [last_end + gap]
        bar_colors.append('#AAAAAA')
        BADs.append(None)

    reduced_bar_colors = []
    for i, color in enumerate(bar_colors):
        if i == 0 or bar_colors[i - 1] != color:
            reduced_bar_colors.append(color)

    borders_for_bars = [1] + borders_to_draw + [ChromosomePosition.chromosomes[chromosome] + 1]
    for i in range(len(borders_for_bars) - 1):
        bar_positions.append((borders_for_bars[i] + borders_for_bars[i + 1]) / 2)
        bar_widths.append(borders_for_bars[i + 1] - borders_for_bars[i])

    for border in segmentation_borders:
        ax.axvline(x=border, ymin=0, ymax=0.5, linestyle='--', color='C4')

    all_borders = [1] + segmentation_borders + [ChromosomePosition.chromosomes[chromosome] + 1]
    for i in range(len(all_borders) - 1):
        if BADs[i]:
            ax.axhline(y=BADs[i],
                       xmin=all_borders[i] / ChromosomePosition.chromosomes[chromosome],
                       xmax=all_borders[i + 1] / ChromosomePosition.chromosomes[chromosome],
                       linewidth=BAD_lw, color=BAD_color,
                       solid_capstyle='butt')

    ax.scatter(x=snps['pos'], y=list(snps['AD']), c=snps['cov'], cmap='BuPu', s=2, vmin=10, vmax=30)
    ax.set_xlim(0, ChromosomePosition.chromosomes[chromosome])
    ax.set_ylim(y_min, y_max)
    ax.grid(which='major', axis='both')
    ax.set_xticklabels([])
    ax.set_yticks(list(range(1, int(y_max) + 1)))
    ax.text(0.99, 0.95, '{}'.format(chromosome),
            horizontalalignment='right',
            verticalalignment='top',
            transform=ax.transAxes)
    ax.set_ylabel('AD')
    # plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0), useMathText=True)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="10%", pad=0.05)
    cax.get_yaxis().set_ticks([])
    cax.set_xlim(1, ChromosomePosition.chromosomes[chromosome])
    cax.set_ylim(0, 1)
    cax.bar(bar_positions, [1] * len(bar_positions), bar_widths, color=reduced_bar_colors, linewidth=0)

    cax.set_xlabel('Chromosome position, bp')
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0), useMathText=True)

    ax.plot([0, 0], [0, 0], color=BAD_color, label='Estimated BAD')
    # ax.legend(loc='center left')

    ax = fig.add_axes([0.95, 0.16, 0.01, 0.75])
    cmap = plt.get_cmap('BuPu')
    norm = m_colors.Normalize(vmin=10, vmax=30)
    m_colorbar.ColorbarBase(ax, cmap=cmap,
                            norm=norm,
                            orientation='vertical')

    plt.savefig(out_path, dpi=300)
    plt.close(fig)

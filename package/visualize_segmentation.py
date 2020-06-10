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
plt.rcParams['figure.figsize'] = 14, 6
plt.rcParams["legend.framealpha"] = 1
plt.rcParams['axes.xmargin'] = 0
plt.rcParams['axes.ymargin'] = 0
plt.rcParams["legend.framealpha"] = 1

# name = 'HCT-116_colon_carcinoma!_labs_richard-myers___biosamples_ENCBS389ENC_'
name = 'HCT-116_colon_carcinoma!_labs_michael-snyder___biosamples_ENCBS626JHZ_'
snps_name = os.path.expanduser('~/Documents/ASB/simulation/' + name + '.tsv')
ploidy_name = os.path.expanduser('~/Documents/ASB/simulation/' + name + '_ploidy.tsv')
cosmic_name = os.path.expanduser('~/Documents/ASB/Cell_lines/cell_lines_copy_number.csv')
cnv_line = 'HCT-116'

snps = pd.read_table(snps_name, header=None)
snps.columns = ['chr', 'pos', 'ID', 'ref', 'alt', 'ref_c', 'alt_c']
snps = snps[snps['ID'] != '.']
snps['AD'] = snps[['ref_c', 'alt_c']].max(axis=1) / snps[['ref_c', 'alt_c']].min(axis=1)
snps['cov'] = snps['ref_c'] + snps['alt_c']
snps['log_cov'] = np.log10(snps['cov'])

BAD_table = pd.read_table(ploidy_name)

chrs = ('chr10', 'chr17')
BAD_color = '#0072B266'
BAD_color_1 = '#0072B2CC'
COSMIC_color = '#D55E00'
BAD_lw = 10
COSMIC_lw = 4
y_min = 0.8
y_max = 6
delta_y = 0.05

# BAD step
fig, (*axs,) = plt.subplots(len(chrs), 1)
fig.tight_layout(pad=1.5)
plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useMathText=True))

for chr, ax in zip(chrs, axs):
    chr_BAD = BAD_table[BAD_table['#chr'] == chr]
    chr_snps = snps[snps['chr'] == chr].copy()
    chr_snps['AD'] = chr_snps['AD'].apply(lambda y: y_max - delta_y if y > y_max else y)

    bar_positions = []
    bar_widths = []
    bar_colors = []
    vd = 1 / 500

    BADs = []

    borders_to_draw = []
    segmentation_borders = []
    last_end = 1
    for index, (pl_chr, start, end, BAD, *values) in chr_BAD.iterrows():
        if start != last_end:
            if last_end == 1:
                borders_to_draw += [start - ChromosomePosition.chromosomes[chr] * vd]
                segmentation_borders += [start - ChromosomePosition.chromosomes[chr] * vd]
            else:
                borders_to_draw += [last_end + ChromosomePosition.chromosomes[chr] * vd,
                                    start - ChromosomePosition.chromosomes[chr] * vd]
                segmentation_borders += [last_end + ChromosomePosition.chromosomes[chr] * vd,
                                         start - ChromosomePosition.chromosomes[chr] * vd]
            bar_colors.append('#AAAAAA')
            BADs.append(None)
        else:
            if last_end != 1:
                segmentation_borders += [last_end]
        last_end = end
        bar_colors.append('C2')
        BADs.append(BAD)
    if last_end != ChromosomePosition.chromosomes[chr] + 1:
        borders_to_draw += [last_end + ChromosomePosition.chromosomes[chr] * vd]
        segmentation_borders += [last_end + ChromosomePosition.chromosomes[chr] * vd]
        bar_colors.append('#AAAAAA')
        BADs.append(None)

    reduced_bar_colors = []
    for i, color in enumerate(bar_colors):
        if i == 0 or bar_colors[i - 1] != color:
            reduced_bar_colors.append(color)

    borders_for_bars = [1] + borders_to_draw + [ChromosomePosition.chromosomes[chr] + 1]
    for i in range(len(borders_for_bars) - 1):
        bar_positions.append((borders_for_bars[i] + borders_for_bars[i + 1]) / 2)
        bar_widths.append(borders_for_bars[i + 1] - borders_for_bars[i])

    for border in segmentation_borders:
        ax.axvline(x=border, ymin=0, ymax=0.5, linestyle='--', color='C4')

    all_borders = [1] + segmentation_borders + [ChromosomePosition.chromosomes[chr] + 1]
    for i in range(len(all_borders) - 1):
        if BADs[i]:
            ax.axhline(y=BADs[i],
                       xmin=all_borders[i] / ChromosomePosition.chromosomes[chr],
                       xmax=all_borders[i + 1] / ChromosomePosition.chromosomes[chr],
                       linewidth=BAD_lw, color=BAD_color_1,
                       solid_capstyle='butt')

    ax.scatter(x=chr_snps['pos'], y=list(chr_snps['AD']), c=chr_snps['cov'], cmap='BuPu', s=2, vmin=10, vmax=30)
    ax.set_xlim(0, ChromosomePosition.chromosomes[chr])
    ax.set_ylim(y_min, y_max)
    ax.grid(which='major', axis='both')
    ax.set_xticklabels([])
    ax.set_yticks(list(range(1, int(y_max) + 1)))
    ax.text(0.99, 0.95, 'Segmentation on {}'.format(chr),
            horizontalalignment='right',
            verticalalignment='top',
            transform=ax.transAxes)
    ax.set_ylabel('AD')
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0), useMathText=True)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="10%", pad=0.05)
    cax.get_yaxis().set_ticks([])
    cax.set_xlim(1, ChromosomePosition.chromosomes[chr])
    cax.set_ylim(0, 1)
    cax.bar(bar_positions, [1] * len(bar_positions), bar_widths, color=reduced_bar_colors, linewidth=0)

cax.set_xlabel('Chromosome position, bp')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0), useMathText=True)

ax = axs[-1]
ax.plot([0, 0], [0, 0], color=BAD_color_1, label='Estimated BAD')
ax.legend(loc='upper left')

ax = fig.add_axes([0.07, 0.87, 0.2, 0.03])
cmap = 'BuPu'
norm = m_colors.Normalize(vmin=10, vmax=30)
cb = m_colorbar.ColorbarBase(ax, cmap=cmap,
                             norm=norm,
                             orientation='horizontal')

plt.savefig(os.path.expanduser('~/AC_4/Figure_AS_4_stage2.svg'), dpi=300)
plt.show()
plt.close(fig)

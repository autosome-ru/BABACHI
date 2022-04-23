import os
import pandas as pd
import matplotlib.colors as m_colors
import matplotlib.colorbar as m_colorbar
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
from zipfile import ZipFile
from .helpers import ChromosomePosition


def init_from_snps_collection(snps_collection, BAD_file,
                              to_zip=False,
                              verbose=True, ext='svg',
                              cosmic_file=None, cosmic_line=None):
    sns.set(font_scale=1.2, style="ticks", font="lato", palette=('#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2',
                                                                 '#D55E00', '#CC79A7'))
    plt.rcParams['font.weight'] = "medium"
    plt.rcParams['axes.labelweight'] = 'medium'
    plt.rcParams['figure.titleweight'] = 'medium'
    plt.rcParams['axes.titleweight'] = 'medium'
    plt.rcParams['figure.figsize'] = 14, 3
    plt.rcParams["legend.framealpha"] = 1
    plt.rcParams['axes.xmargin'] = 0
    plt.rcParams['axes.ymargin'] = 0

    BAD_table = pd.read_table(BAD_file)
    file_name = os.path.splitext(os.path.basename(BAD_file))[0]
    out_path = os.path.join(os.path.dirname(BAD_file), '{}_visualization'.format(file_name))
    if not os.path.isdir(out_path):
        os.mkdir(out_path)

    cosmics = read_cosmic(cosmic_file, cosmic_line=cosmic_line)
    for chromosome in snps_collection.keys():
        if verbose:
            print('Visualizing {}'.format(chromosome))

        result = filter_data_by_chromosome(chromosome, BAD_table=BAD_table,
                                           snps_collection=snps_collection,
                                           cosmics=cosmics,
                                           )
        if result is None:
            continue
        BAD_segments, snps, cosmic_data = result
        visualize_chromosome(chromosome, BAD_segments, snps,
                             os.path.join(out_path,
                                          '{}_{}.{}'.format(file_name,
                                                            chromosome, ext)),
                             cosmic_data)
    if to_zip:
        with ZipFile(out_path + '.zip', 'w') as zip_archive:
            for file in os.listdir(out_path):
                if file.endswith(ext):
                    zip_archive.write(os.path.join(out_path, file))


def read_cosmic(cosmic_file, cosmic_line):
    if cosmic_file is not None and cosmic_line is not None:
        cosmic = pd.read_table(cosmic_file, low_memory=False, header=None)
        cosmic.columns = ['sample_name', '#chr', 'startpos', 'endpos', 'minorCN', 'totalCN']
        cosmic = cosmic[
            (cosmic['sample_name'] == cosmic_line) &
            (cosmic['minorCN'] != 0)
            ]
        if cosmic.empty:
            return None
        cosmic['BAD'] = cosmic.eval('(totalCN - minorCN) / minorCN')
        cosmic['startpos'] = cosmic['startpos'].astype(int)
        cosmic['endpos'] = cosmic['endpos'].astype(int)
        cosmic = cosmic.drop(['totalCN', 'minorCN', 'sample_name'], axis=1)
        return cosmic
    return None

def filter_data_by_chromosome(chromosome, BAD_table, snps_collection=None,
                              cosmics=None):
    BAD_segments = BAD_table[BAD_table['#chr'] == chromosome]
    if BAD_segments.empty:
        return None
    if snps_collection is not None:
        column_names = ['pos', 'ref_c', 'alt_c']
        snps = pd.DataFrame(dict(zip(column_names, zip(*snps_collection[chromosome]))))
        if snps.empty:
            return None
        snps['AD'] = snps[['ref_c', 'alt_c']].max(axis=1) / snps[['ref_c', 'alt_c']].min(axis=1)
        snps['cov'] = snps.eval('ref_c + alt_c')
    else:
        snps = None
    cosmic = cosmics[cosmics['#chr'] == chromosome] if cosmics is not None else None
    return BAD_segments, snps, cosmic


def setup_plot():
    fig, ax = plt.subplots()
    fig.tight_layout(rect=[0, 0.01, 0.95, 1])
    plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useMathText=True))
    return fig, ax


def post_draw_settings(ax, chromosome, y_min=0.8, y_max=6):
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
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0), useMathText=True)


def visualize_chromosome(chromosome, BAD_segments,
                         snps=None, out_path=None, chr_cosmic=None):
    BAD_color = '#0072B2CC'
    COSMIC_color = '#D55E00'
    BAD_lw = 10
    COSMIC_lw = 4
    y_min = 0.8
    y_max = 6
    delta_y = 0.05
    fig, ax = setup_plot()
    if snps is not None:
        add_snps(ax, snps, y_max, delta_y)
    add_babachi_estimations(fig, ax, chromosome,
                            BAD_segments, chr_cosmic,
                            BAD_color=BAD_color,
                            COSMIC_color=COSMIC_color,
                            BAD_lw=BAD_lw,
                            COSMIC_lw=COSMIC_lw)
    post_draw_settings(ax, chromosome, y_min=y_min, y_max=y_max)
    if out_path:
        plt.savefig(out_path, dpi=300)
    else:
        plt.show()
    plt.close(fig)


def add_snps(ax, snps, y_max=6, delta_y=0.05):
    ref_snps = snps[snps['ref_c'] >= snps['alt_c']]
    alt_snps = snps[snps['ref_c'] < snps['alt_c']]
    snps['AD'] = snps['AD'].apply(lambda y: y_max - delta_y if y > y_max else y)
    for snp_df, cmap in zip((ref_snps, alt_snps), ('BuGn', 'BuPu')):
        ax.scatter(x=snp_df['pos'], y=list(snp_df['AD']),
                   c=snp_df['cov'], cmap=cmap, s=2, vmin=10, vmax=30)


def add_babachi_estimations(fig, ax, chromosome, BAD_segments, cosmic_est=None,
                            BAD_color='#0072B2CC', COSMIC_color='#D55E00',
                            BAD_lw=10, COSMIC_lw=4):
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

    # cosmic
    if cosmic_est is not None and not cosmic_est.empty:
        cosmic_bar_colors = []
        vd = 1 / 500
        COSMIC_BADs = []
        cosmic_borders = []
        last_end = 1
        for index, (chrom, startpos, endpos, BAD) in cosmic_est.iterrows():
            if startpos - last_end >= ChromosomePosition.chromosomes[chromosome] * vd * 2:
                if last_end == 1:
                    cosmic_borders += [startpos - ChromosomePosition.chromosomes[chromosome] * vd]
                else:
                    cosmic_borders += [last_end + ChromosomePosition.chromosomes[chromosome] * vd,
                                       startpos - ChromosomePosition.chromosomes[chromosome] * vd]
                cosmic_bar_colors.append('#AAAAAA')
                COSMIC_BADs.append(None)
            else:
                if last_end != 1:
                    cosmic_borders += [last_end]
            last_end = endpos
            cosmic_bar_colors.append('C2')
            COSMIC_BADs.append(BAD)
        if last_end != ChromosomePosition.chromosomes[chromosome] + 1:
            cosmic_borders += [last_end + ChromosomePosition.chromosomes[chromosome] * vd]
            cosmic_bar_colors.append('#AAAAAA')
            COSMIC_BADs.append(None)

        all_cosmic_borders = [1] + cosmic_borders + [ChromosomePosition.chromosomes[chromosome] + 1]

        for i in range(len(all_cosmic_borders) - 1):
            if COSMIC_BADs[i]:
                ax.axhline(y=COSMIC_BADs[i],
                           xmin=all_cosmic_borders[i] / ChromosomePosition.chromosomes[chromosome],
                           xmax=all_cosmic_borders[i + 1] / ChromosomePosition.chromosomes[chromosome],
                           linewidth=COSMIC_lw, color=COSMIC_color, snap=False, ms=0, mew=0,
                           solid_capstyle='butt')

    ax.plot([0, 0], [0, 0], color=BAD_color, label='Estimated BAD')
    if cosmic_est is not None and not cosmic_est.empty:
        ax.plot([0, 0], [0, 0], color=COSMIC_color, label='COSMIC BAD')
    # ax.legend(loc='center left')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="10%", pad=0.05)
    cax.get_yaxis().set_ticks([])
    cax.set_xlim(1, ChromosomePosition.chromosomes[chromosome])
    cax.set_ylim(0, 1)
    cax.bar(bar_positions, [1] * len(bar_positions), bar_widths, color=reduced_bar_colors, linewidth=0)

    cax.set_xlabel('Chromosome position, bp')
    ax = fig.add_axes([0.95, 0.16, 0.01, 0.75])
    cmap = plt.get_cmap('BuPu')
    norm = m_colors.Normalize(vmin=10, vmax=30)
    m_colorbar.ColorbarBase(ax, cmap=cmap,
                            norm=norm,
                            orientation='vertical')

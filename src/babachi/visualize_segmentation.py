import os
import pandas as pd
import matplotlib.colors as m_colors
import matplotlib.colorbar as m_colorbar
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
from zipfile import ZipFile
from babachi.chrom_wrapper import init_wrapper
from babachi.helpers import init_wrapper
from babachi.models import GenomeSNPsHandler


class BabachiVisualizer:
    def __init__(self, chromosomes_wrapper=None):
        self.chromosomes_wrapper = init_wrapper(chromosomes_wrapper)

    def init_from_snps_collection(self, snps_collection: GenomeSNPsHandler, BAD_file,
                                  to_zip=False,
                                  verbose=True, ext='svg',
                                  cosmic_file=None, cosmic_line=None):
        sns.set(font_scale=1.2, style="ticks", font="lato",
                palette=('#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2',
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
        out_path = os.path.join(os.path.dirname(BAD_file), '{}.visualization'.format(file_name))
        if not os.path.isdir(out_path):
            os.mkdir(out_path)

        cosmics = self.read_cosmic(cosmic_file, cosmic_line=cosmic_line)
        for chromosome in snps_collection.chromosomes_order:
            if verbose:
                pass
                # print('Visualizing {}'.format(chromosome))

            result = self.filter_data_by_chromosome(chromosome, BAD_table=BAD_table,
                                                    snps_collection=snps_collection,
                                                    cosmics=cosmics,
                                                    )
            if result is None:
                continue
            BAD_segments, snps, cosmic_data = result
            self.visualize_chromosome(chromosome, BAD_segments, snps,
                                      os.path.join(out_path,
                                                   f'{file_name}_{chromosome}.{ext}'),
                                      cosmic_data)
        if to_zip:
            with ZipFile(out_path + '.zip', 'w') as zip_archive:
                for file in os.listdir(out_path):
                    if file.endswith(ext):
                        zip_archive.write(os.path.join(out_path, file))

    def visualize_chromosome(self, chromosome, BAD_segments,
                             snps=None, out_path=None, chr_cosmic=None):
        chrVis = ChromosomeVisualizer(chromosome_name=chromosome,
                                      chromosome_length=self.chromosomes_wrapper.chromosomes[chromosome],
                                      BAD_segments=BAD_segments,
                                      snps=snps, out_path=out_path, chr_cosmic=chr_cosmic
                                      )
        chrVis.visualize_chromosome()

    @staticmethod
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

    @staticmethod
    def filter_data_by_chromosome(chromosome, BAD_table, snps_collection: GenomeSNPsHandler = None,
                                  cosmics=None):
        BAD_segments = BAD_table[BAD_table['#chr'] == chromosome]
        if BAD_segments.empty:
            return None
        if snps_collection is not None:
            column_names = ['pos', 'ref_c', 'alt_c']
            snps = pd.DataFrame(dict(zip(column_names, zip(*snps_collection.data[chromosome].data.transpose()))))
            if snps.empty:
                return None
            snps['AD'] = snps[['ref_c', 'alt_c']].max(axis=1) / snps[['ref_c', 'alt_c']].min(axis=1)
            snps['cov'] = snps.eval('ref_c + alt_c')
        else:
            snps = None
        cosmic = cosmics[cosmics['#chr'] == chromosome] if cosmics is not None else None
        return BAD_segments, snps, cosmic


class ChromosomeVisualizer:
    def __init__(self, chromosome_name, chromosome_length, BAD_segments,
                 snps=None, out_path=None, chr_cosmic=None):
        self.chromosome_name = chromosome_name
        self.chromosome_length = chromosome_length
        self.BAD_segments = BAD_segments
        self.snps = snps
        self.out_path = out_path
        self.cosmic_est = chr_cosmic
        self.fig = None
        self.ax = None

    def setup_plot(self):
        fig, ax = plt.subplots()
        fig.tight_layout(rect=[0, 0.01, 0.95, 1])
        plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useMathText=True))
        self.fig = fig
        self.ax = ax

    def post_draw_settings(self, y_min=0.8, y_max=6):
        self.ax.set_xlim(0, self.chromosome_length)
        self.ax.set_ylim(y_min, y_max)
        self.ax.grid(which='major', axis='both')
        self.ax.set_xticklabels([])
        self.ax.set_yticks(list(range(1, int(y_max) + 1)))
        self.ax.text(0.99, 0.95, '{}'.format(self.chromosome_name),
                     horizontalalignment='right',
                     verticalalignment='top',
                     transform=self.ax.transAxes)
        self.ax.set_ylabel('AD')
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0), useMathText=True)

    def visualize_chromosome(self, BAD_color='#0072B2CC', COSMIC_color='#D55E00', BAD_lw=10,
                             COSMIC_lw=4, y_min=0.8, y_max=6, delta_y=0.05):
        # BAD_color = '#0072B2CC'
        # COSMIC_color = '#D55E00'
        # BAD_lw = 10
        # COSMIC_lw = 4
        # y_min = 0.8
        # y_max = 6
        # delta_y = 0.05
        self.setup_plot()
        if self.snps is not None:
            self.add_snps(y_max, delta_y)
        self.add_babachi_estimations(
            BAD_color=BAD_color,
            COSMIC_color=COSMIC_color,
            BAD_lw=BAD_lw,
            COSMIC_lw=COSMIC_lw)
        self.post_draw_settings(y_min=y_min, y_max=y_max)
        if self.out_path:
            plt.savefig(self.out_path, dpi=300)
        else:
            plt.show()
        plt.close(self.fig)

    def add_snps(self, y_max=6, delta_y=0.05):
        self.snps['AD'] = self.snps['AD'].apply(lambda y: y_max - delta_y if y > y_max else y)
        ref_more_than_alt = self.snps.eval('ref_c >=alt_c')
        ref_snps = self.snps[ref_more_than_alt]
        alt_snps = self.snps[~ref_more_than_alt]
        for snp_df, cmap in zip((ref_snps, alt_snps), ('BuGn', 'BuPu')):
            self.ax.scatter(x=snp_df['pos'], y=list(snp_df['AD']),
                            c=snp_df['cov'], cmap=cmap, s=2, vmin=10, vmax=30)

    def add_babachi_estimations(self,
                                BAD_color='#0072B2CC', COSMIC_color='#D55E00',
                                BAD_lw=10, COSMIC_lw=4):
        bar_positions = []
        bar_widths = []
        bar_colors = []
        gap = 1 / 500 * self.chromosome_length

        BADs = []

        borders_to_draw = []
        segmentation_borders = []
        last_end = 1
        for index, (pl_chr, start, end, BAD, *values) in self.BAD_segments.iterrows():
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
        if last_end != self.chromosome_length + 1:
            borders_to_draw += [last_end + gap]
            segmentation_borders += [last_end + gap]
            bar_colors.append('#AAAAAA')
            BADs.append(None)

        reduced_bar_colors = []
        for i, color in enumerate(bar_colors):
            if i == 0 or bar_colors[i - 1] != color:
                reduced_bar_colors.append(color)

        borders_for_bars = [1] + borders_to_draw + [self.chromosome_length + 1]
        for i in range(len(borders_for_bars) - 1):
            bar_positions.append((borders_for_bars[i] + borders_for_bars[i + 1]) / 2)
            bar_widths.append(borders_for_bars[i + 1] - borders_for_bars[i])

        for border in segmentation_borders:
            self.ax.axvline(x=border, ymin=0, ymax=0.5, linestyle='--', color='C4')

        all_borders = [1] + segmentation_borders + [self.chromosome_length + 1]
        for i in range(len(all_borders) - 1):
            if BADs[i]:
                self.ax.axhline(y=BADs[i],
                                xmin=all_borders[i] / self.chromosome_length,
                                xmax=all_borders[i + 1] / self.chromosome_length,
                                linewidth=BAD_lw, color=BAD_color,
                                solid_capstyle='butt')

        # cosmic
        if self.cosmic_est is not None and not self.cosmic_est.empty:
            cosmic_bar_colors = []
            vd = 1 / 500
            COSMIC_BADs = []
            cosmic_borders = []
            last_end = 1
            for index, (chrom, startpos, endpos, BAD) in self.cosmic_est.iterrows():
                if startpos - last_end >= self.chromosome_length * vd * 2:
                    if last_end == 1:
                        cosmic_borders += [startpos - self.chromosome_length * vd]
                    else:
                        cosmic_borders += [last_end + self.chromosome_length * vd,
                                           startpos - self.chromosome_length * vd]
                    cosmic_bar_colors.append('#AAAAAA')
                    COSMIC_BADs.append(None)
                else:
                    if last_end != 1:
                        cosmic_borders += [last_end]
                last_end = endpos
                cosmic_bar_colors.append('C2')
                COSMIC_BADs.append(BAD)
            if last_end != self.chromosome_length + 1:
                cosmic_borders += [last_end + self.chromosome_length * vd]
                cosmic_bar_colors.append('#AAAAAA')
                COSMIC_BADs.append(None)

            all_cosmic_borders = [1] + cosmic_borders + [self.chromosome_length + 1]

            for i in range(len(all_cosmic_borders) - 1):
                if COSMIC_BADs[i]:
                    self.ax.axhline(y=COSMIC_BADs[i],
                                    xmin=all_cosmic_borders[i] / self.chromosome_length,
                                    xmax=all_cosmic_borders[i + 1] / self.chromosome_length,
                                    linewidth=COSMIC_lw, color=COSMIC_color, snap=False, ms=0, mew=0,
                                    solid_capstyle='butt')

        self.ax.plot([0, 0], [0, 0], color=BAD_color, label='Estimated BAD')
        if self.cosmic_est is not None and not self.cosmic_est.empty:
            self.ax.plot([0, 0], [0, 0], color=COSMIC_color, label='COSMIC BAD')
        # ax.legend(loc='center left')
        divider = make_axes_locatable(self.ax)
        cax = divider.append_axes("bottom", size="10%", pad=0.05)
        cax.get_yaxis().set_ticks([])
        cax.set_xlim(1, self.chromosome_length)
        cax.set_ylim(0, 1)
        cax.bar(bar_positions, [1] * len(bar_positions), bar_widths, color=reduced_bar_colors, linewidth=0)

        cax.set_xlabel('Chromosome position, bp')
        ax = self.fig.add_axes([0.95, 0.16, 0.01, 0.75])
        cmap = plt.get_cmap('BuPu')
        norm = m_colors.Normalize(vmin=10, vmax=30)
        m_colorbar.ColorbarBase(ax, cmap=cmap,
                                norm=norm,
                                orientation='vertical')

"""
count step
"""

import os,sys
import random
#import subprocess
#import sys
from collections import defaultdict
from itertools import groupby

import numpy as np
import pandas as pd
import pysam

import accuracode.tools.utils as utils
from accuracode.tools.__init__ import __PATTERN_DICT__
from accuracode.tools.step import Step, s_common
from accuracode.accura.mkref import parse_genomeDir_rna

TOOLS_DIR = os.path.dirname(__file__)
random.seed(0)
np.random.seed(0)


class Count(Step):
    """
    Features
    - Get read/UMI/gene counts for each well barcode. 

    - Generate expression matrix.

    Output
    - `{sample}_matrix.tsv.gz` The expression matrix of the well barcode, separated by tabs. 

    - `{sample}_count_detail.txt.gz` 4 columns: 
        - barcode  
        - gene ID  
        - UMI count  
        - read_count  

    - `{sample}_counts.txt` 5 columns:
        - well: Oligo id
        - Barcode: Well barcode sequence
        - readcount: read count of each well barcode
        - UMI: UMI count for each well barcode
        - gene: gene count for each well barcode

    - `{sample}_counts_raw.txt`, same as {sample}_counts.txt, including all wells without UMI filter.

    - `{sample}_downsample.txt` 2 columns：
        - percent: percentage of sampled reads
        - median_geneNum: median gene number per cell


    """

    def __init__(self, args, step):
        Step.__init__(self, args, step)
        if args.chemistry == 'customized':
            if args.whitelist:
                self.wellFile = args.whitelist
            else:
                print('Customized chemistry must provide whitelist, please check files.')
                sys.exit(1)
        else:
            self.wellFile = self.get_scope_bcwell(args.chemistry)
        self.bam = args.bam

        # set
        self.gtf_file = parse_genomeDir_rna(args.genomeDir)['gtf']
        self.id_name = utils.get_id_name_dict(self.gtf_file)
        self.well_name = self.bc2well(self.wellFile)
        self.umi_cutoff = args.UMI_cutoff
        if args.skip_umi_correct:
            self.umi_correct = False
        else:
            self.umi_correct = True

        # output files
        self.count_detail_file = f'{self.outdir}/{self.sample}_count_detail.txt'
        self.marked_count_file = f'{self.outdir}/{self.sample}_counts.txt'
        self.raw_count_file = f'{self.outdir}/{self.sample}_counts_raw.txt'
        self.matrix_file = f'{self.outdir}/{self.sample}_matrix.tsv.gz'
        self.stat_file = f'{self.outdir}/stat.txt'
        self.downsample_file = f'{self.outdir}/{self.sample}_downsample.txt'
        
    def run(self):
        self.bam2table()
        df = pd.read_table(self.count_detail_file, header=0)

        # df_sum
        Count.get_df_sums(df, self.well_name, self.raw_count_file, self.marked_count_file, self.stat_file, self.umi_cutoff)
        #Count.get_df_sum(df, self.marked_count_file, self.stat_file)

        # output matrix
        self.exp_matrix(df, self.matrix_file, self.well_name)

        # summary
        tbltxt = pd.read_csv(self.marked_count_file,header=0,sep="\t")
        keep_bc = tbltxt['Barcode'].to_list()
        tbltxt.drop(columns='Barcode',inplace=True)
        tbldiv = self.stat_table(tbltxt)

        # downsampling
        df_ds = df[df['Barcode'].isin(keep_bc)]
        res_dict = self.downsample(df_ds)

        # report
        self.report_prepare(tbldiv)
        self.add_content_item('metric', downsample_summary=res_dict)
        self.clean_up()


    def report_prepare(self,outable):
        self.add_data_item(accura_count=outable)
        df0 = pd.read_table(self.downsample_file, header=0)
        self.add_data_item(percentile=df0['percent'].tolist())
        self.add_data_item(MedianGeneNum=df0['median_geneNum'].tolist())


    @staticmethod
    def correct_umi(umi_dict, percent=0.1):
        """
        Correct umi_dict in place.

        Args:
            umi_dict: {umi_seq: umi_count}
            percent: if hamming_distance(low_seq, high_seq) == 1 and
                low_count / high_count < percent, merge low to high.

        Returns:
            n_corrected_umi: int
            n_corrected_read: int
        """
        n_corrected_umi = 0
        n_corrected_read = 0

        # sort by value(UMI count) first, then key(UMI sequence)
        umi_arr = sorted(
            umi_dict.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)
        while True:
            # break when only highest in umi_arr
            if len(umi_arr) == 1:
                break
            umi_low = umi_arr.pop()
            low_seq = umi_low[0]
            low_count = umi_low[1]

            for umi_kv in umi_arr:
                high_seq = umi_kv[0]
                high_count = umi_kv[1]
                if float(low_count / high_count) > percent:
                    break
                if utils.hamming_distance(low_seq, high_seq) == 1:
                    n_low = umi_dict[low_seq]
                    n_corrected_umi += 1
                    n_corrected_read += n_low
                    # merge
                    umi_dict[high_seq] += n_low
                    del (umi_dict[low_seq])
                    break
        return n_corrected_umi, n_corrected_read

    @utils.add_log
    def bam2table(self):
        """
        bam to detail table
        must be used on name_sorted bam
        """
        samfile = pysam.AlignmentFile(self.bam, "rb")
        with open(self.count_detail_file, 'wt') as fh1:
            fh1.write('\t'.join(['Barcode', 'geneID', 'UMI', 'count']) + '\n')

            def keyfunc(x):
                return x.query_name.split('_', 1)[0]
            for _, g in groupby(samfile, keyfunc):
                gene_umi_dict = defaultdict(lambda: defaultdict(int))
                for seg in g:
                    (barcode, umi) = seg.query_name.split('_')[:2]
                    if not seg.has_tag('XT'):
                        continue
                    gene_id = seg.get_tag('XT')
                    gene_umi_dict[gene_id][umi] += 1
                if self.umi_correct:
                    for gene_id in gene_umi_dict:
                        Count.correct_umi(gene_umi_dict[gene_id])

                # output
                for gene_id in gene_umi_dict:
                    for umi in gene_umi_dict[gene_id]:
                        fh1.write('%s\t%s\t%s\t%s\n' % (barcode, gene_id, umi,
                                                        gene_umi_dict[gene_id][umi]))
        samfile.close()



    def bc2well(self, bc_well_file):
        wells = {}
        with open(bc_well_file) as f:
            for i in f:
                ii= i.strip().split()
                if len(ii) == 1:
                    wells[ii[0]] = ii[0]
                else:
                    wells[ii[0]] = ii[1]
        return wells


    @staticmethod
    def get_df_sums(df, well_name, raw_count_file, marked_count_file, stat_file, umi_cutoff, col='UMI'):
    #def get_df_sum(df,  marked_count_file, stat_file, col='UMI'):

        def format_int(number):
            if isinstance (number,(int,float)):
                return str("{:,}".format(round(number)))
            else:
                return number

        df_sum_raw = df.groupby('Barcode').agg({
            'count': 'sum',
            'UMI': 'count',
            'geneID': 'nunique'
        })

        df_sum_raw.columns = ['readcount', 'UMI', 'gene']
        df_sum_raw = df_sum_raw.sort_values(col, ascending=False)
        id1 = df_sum_raw.index.to_series()
        name = id1.apply(lambda x: well_name[x])
        df_sum_raw.reset_index(inplace=True)
        df_sum_raw.index = name
        df_sum_raw.index.name = "well"

        df_sum = df_sum_raw[df_sum_raw['UMI']>=umi_cutoff]

        stat_header=['Median Reads per Well','Median UMI per Well','Median Genes per Well','Mean Reads per Well','Mean UMI per Well','Mean Genes per Well']
        stat_info=[format_int(df_sum['readcount'].median()), format_int(df_sum['UMI'].median()), format_int(df_sum['gene'].median()), format_int(df_sum['readcount'].mean()), format_int(df_sum['UMI'].mean()), format_int(df_sum['gene'].mean())]
        df_stat = pd.DataFrame({"item": stat_header, "count": stat_info})
        df_stat.to_csv(stat_file, sep=":", header=None, index=False)

        df_sum_raw = df_sum_raw.applymap(lambda x: format_int(x))
        df_sum = df_sum.applymap(lambda x: format_int(x))
        #df_sum = df_sum.sort_index(axis=0, ascending=True)
        df_sum_raw.to_csv(raw_count_file, sep='\t')
        df_sum.to_csv(marked_count_file, sep='\t')



    @utils.add_log
    def exp_matrix(self, df,  matrix_table_file, well_name):
    #def exp_matrix(self, df,  matrix_table_file):
        table = df.pivot_table(
                index='geneID', columns='Barcode', values='UMI',
                aggfunc=len).fillna(0).astype(int)

        # convert id to name; write table matrix
        id1 = table.index.to_series()
        name = id1.apply(lambda x: self.id_name[x])
        table.index = name
        table.index.name = ""
        bc = table.columns.to_series()
        sample = bc.apply(lambda x: well_name[x])
        table.columns= sample
        table = table.sort_index(axis=1, ascending=True)

        table.to_csv(
            matrix_table_file,
            sep="\t",
            compression='gzip')

    def stat_table(self,txt):
        marker_gene_table = txt.to_html(
                escape=False,
                index=False,
                table_id='accura_table',
                justify="center")

        return marker_gene_table

    def get_scope_bcwell(self,chemistry):
        import accuracode
        root_path = os.path.dirname(accuracode.__file__)
        bcwell_f = f'{root_path}/data/chemistry/{chemistry}/bcwell'
        return bcwell_f


    @staticmethod
    def sub_sample(fraction, df_cell, cell_read_index):
        cell_read = df_cell['count'].sum()
        frac_n_read = int(cell_read * fraction)
        subsample_read_index = cell_read_index[:frac_n_read]
        #index_dedup, counts = np.unique(subsample_read_index, return_counts=True)
        index_dedup = np.unique(subsample_read_index, return_counts=False)
        #n_count_once = np.sum(counts == 1)
        # total = UMI
        #umi_total = len(index_dedup)
        #umi_saturation = round((1 - n_count_once / umi_total) * 100, 2)
        #read_total = frac_n_read
        #read_saturation = round((1 - n_count_once / read_total) * 100, 2)

        # gene median
        df_cell_subsample = df_cell.loc[index_dedup, ]
        geneNum_median = float(df_cell_subsample.groupby(
            'Barcode').agg({'geneID': 'nunique'}).median())

        return geneNum_median


    @utils.add_log
    def downsample(self, df_cell):
        """saturation and median gene
        return fraction=1 saturation
        """
        cell_read_index = np.array(df_cell.index.repeat(df_cell['count']), dtype='int32')
        np.random.shuffle(cell_read_index)

        format_str = "%.2f\t%.2f\n"
        res_dict = {
            "fraction": [],
            "median_gene": []
        }
        with open(self.downsample_file, 'w') as fh:
            fh.write('percent\tmedian_geneNum\n')
            fh.write(format_str % (0, 0))
            for fraction in np.arange(0.1, 1.1, 0.1):
                geneNum_median = Count.sub_sample(
                    fraction, df_cell, cell_read_index)
                fh.write(format_str % (fraction, geneNum_median))
                #def format_float(x): return round(x / 100, 4)
                res_dict["fraction"].append(round(fraction, 1))
                res_dict["median_gene"].append(geneNum_median)

        return res_dict



@utils.add_log
def count(args):
    step_name = "count"
    runner = Count(args, step_name)
    runner.run()


def get_opts_count(parser, sub_program):
    parser.add_argument('--genomeDir', help='Required. Genome directory.')
    parser.add_argument('--chemistry', 
        help="""Predefined (pattern, barcode whitelist, linker whitelist) combinations. Can be one of:
- `accuracode96` Used for AccuraCode96 libraries.
- `accuracode384` Used for AccuraCode384 libraries.
- `customized` Used for user defined combinations. You need to provide `pattern`, `whitelist` and  `linker` at the
same time.""",
        choices=list(__PATTERN_DICT__.keys()),
        default='accuracode384'
    )

    parser.add_argument(
        '--whitelist',
        help='Cell barcode whitelist file path, one cell barcode per line.'
    )

    parser.add_argument(
        '--UMI_cutoff', default=500, type=int,
        help='UMI cutoff for output, default 500.'
    )

    parser.add_argument("--skip_umi_correct", help="Skip umi correction", action='store_true')

    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--bam', help='Required. BAM file from featureCounts.', required=True)

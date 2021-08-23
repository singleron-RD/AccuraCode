import os

__VERSION__ = "1.0.0"
__version__ = __VERSION__

ASSAY_DICT = {
    'accura': 'AccuraCode'
}

ROOT_PATH = os.path.dirname(__file__)

RELEASED_ASSAYS = ['accura']

HELP_DICT = {
    'gene_list': 'Required. Gene list file, one gene symbol per line. Only results of these genes are reported.',
    'genomeDir': 'Genome directory after running `mkref`.',
    'thread': 'Thread to use.',
    'debug': 'If this argument is used, accuracode may output addtional file for debugging.',
    'fasta': 'Required. Genome fasta file. Use relative path to `genomeDir`.',
    'outdir': 'Output directory.'
}

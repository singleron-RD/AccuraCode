import os

__VERSION__ = "1.1.0"
__version__ = __VERSION__

ASSAY_DICT = {
    'accura': 'AccuraCode_rna'
}

ROOT_PATH = os.path.dirname(__file__)

RELEASED_ASSAYS = ['accura']

HELP_DICT = {
    'gene_list': 'Required. Gene list file, one gene symbol per line. Only results of these genes are reported.',
    'genomeDir': 'Genome directory after running `mkref`.',
    'thread': 'Thread to use.',
    'debug': 'If this argument is used, accuracode may output addtional file for debugging.',
    'fasta': 'Required. Genome fasta file. Use relative path to `genomeDir`.',
    'outdir': 'Output directory.',
    'chemistry': 'You need to explicitly use `--chemistry accuracode96 or accuracode384 (default)`. `--chemistry customized` is used for user defined combinations that you need to provide `--pattern`, `--whitelist` and `--linker` at the same time.'
}

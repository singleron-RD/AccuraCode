## Arguments
`--genomeDir` Default='./'. Output directory.

`--thread` Default=6. Threads to use.

`--genome_name` Required, genome name.

`--dry_run` Only write config file and exit.

`--fasta` Required. Genome fasta file. Must be relative file path to genomeDir.

`--gtf` Required. Genome gtf file. Must be relative file path to genomeDir.

`--mt_gene_list` Mitochondria gene list file. Must be relative file path to genomeDir.
It is a plain text file with one gene per line. 
If not provided, will use `MT-` and `mt-` to determine mitochondria genes.


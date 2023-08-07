## Features
## Output files
## Arguments
`--mod` mod, sjm or shell.

`--mapfile` tsv file, 4 columns:
                1st col: LibName;
                2nd col: DataDir;
                3rd col: SampleName;
                4th col: optional;.

`--rm_files` remove redundant fq.gz and bam after running.

`--steps_run` Steps to run. Multiple Steps are separated by comma.

`--outdir` Output directory.

`--thread` Thread to use.

`--debug` If this argument is used, accuracode may output addtional file for debugging.

`--pattern` The pattern of R1 reads, e.g. `C9U12T18`. The number after the letter represents the number 
        of bases.  
- `C`: cell barcode  
- `L`: linker(common sequences, can be absent in pattern design)  
- `U`: UMI    
- `T`: poly T.

`--linker` Linker whitelist file path, one linker per line.

`--lowQual` Default 0. Bases in cell barcode and UMI whose phred value are lower than lowQual will be regarded as low-quality bases.

`--lowNum` The maximum allowed lowQual bases in cell barcode and UMI.

`--nopolyT` Outputs R1 reads without polyT.

`--noLinker` Outputs R1 reads without correct linker.

`--allowNoPolyT` Allow valid reads without polyT.

`--allowNoLinker` Allow valid reads without correct linker.

`--gzip` Output gzipped fastq files.

`--adapter_fasta` Addtional adapter fasta file.

`--minimum_length` Default `20`. Discard processed reads that are shorter than LENGTH.

`--nextseq_trim` Default `20`. Quality trimming of reads using two-color chemistry (NextSeq). 
Some Illumina instruments use a two-color chemistry to encode the four bases. 
This includes the NextSeq and the NovaSeq. 
In those instruments, a ‘dark cycle’ (with no detected color) encodes a G. 
However, dark cycles also occur when sequencing “falls off” the end of the fragment.
The read then contains a run of high-quality, but incorrect “G” calls at its 3’ end.

`--overlap` Default `10`. Since Cutadapt allows partial matches between the read and the adapter sequence,
short matches can occur by chance, leading to erroneously trimmed bases. 
For example, roughly 0.25 of all reads end with a base that is identical to the first base of the adapter. 
To reduce the number of falsely trimmed bases, the alignment algorithm requires that 
at least {overlap} bases match between adapter and read.

`--insert` Default `150`. Read2 insert length.

`--outFilterMatchNmin` Default `0`. Alignment will be output only if the number of matched bases 
is higher than or equal to this value.

`--out_unmapped` Output unmapped reads.

`--STAR_param` Other STAR parameters.

`--outFilterMultimapNmax` Default `1`. How many places are allowed to match a read at most.

`--starMem` Default `30`. Maximum memory that STAR can use.

`--gtf_type` Specify feature type in GTF annotation.

`--genomeDir` Required. Genome directory.

`--chemistry` Predefined (pattern, barcode whitelist, linker whitelist) combinations. Can be one of:
- `accuracode96` Used for AccuraCode96 libraries.
- `accuracode384` Used for AccuraCode384 libraries.
- `customized` Used for user defined combinations. You need to provide `pattern`, `whitelist`, `linker`(if designed) at the
same time.

`--whitelist` Cell barcode whitelist file path, one cell barcode per line or two columns(`barcode reportname`) per line.

`--UMI_cutoff` UMI cutoff for output, default 500.

`--skip_umi_correct` For big data to skip umi correct.


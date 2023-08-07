## Arguments
`--chemistry` Predefined (pattern, barcode whitelist, linker whitelist) combinations. Can be one of:  
- `accuracode96` Used for AccuraCode96 libraries.  
- `accuracode384` Used for AccuraCode384 libraries.  
- `customized` Used for user defined combinations. You need to provide `pattern`, `whitelist`, `linker`(if designed) at the 
same time.

`--pattern` The pattern of R1 reads, e.g. `C9U12T18`. The number after the letter represents the number 
        of bases.  
- `C`: cell barcode  
- `L`: linker(common sequences, can be absent in pattern design)  
- `U`: UMI    
- `T`: poly T.

`--whitelist` Cell barcode whitelist file path, one cell barcode per line.

`--linker` Linker whitelist file path, one linker per line.

`--lowQual` Default 0. Bases in cell barcode and UMI whose phred value are lower than lowQual will be regarded as low-quality bases.

`--lowNum` The maximum allowed lowQual bases in cell barcode and UMI.

`--nopolyT` Outputs R1 reads without polyT.

`--noLinker` Outputs R1 reads without correct linker.

`--allowNoPolyT` Allow valid reads without polyT.

`--allowNoLinker` Allow valid reads without correct linker.

`--gzip` Output gzipped fastq files.

`--fq1` R1 fastq file. Multiple files are separated by comma.

`--fq2` R2 fastq file. Multiple files are separated by comma.

`--outdir` Output diretory.

`--assay` Assay name.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, accuracode may output addtional file for debugging.


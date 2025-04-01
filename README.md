# DGR discovery
## Requirements
blast 2.15.0+
pandas 2.2.1+
biopython 1.83+
pyhmmer 0.10.9+
pyrodigal-gv 0.3.1+
rich 13.7.1+

## How to run code

    python DGRfinder.py -i [gzipped input contig file] -o [output dgr file .csv] -c [output contig file .fna] -cpu [num threads] -hmm resources/hmmprofiles/DGR_RTs_MSA

## Parameters
- Input file: 
	- Should be in fasta.gz format
- Outputs:
	- DGR file will contain information about identified DGRs
	- Contig file will contain whole contigs from which at least one DGR was identified
- hmm:
	- Precompiled HMM to identify DGR RTs

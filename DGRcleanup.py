from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline

import numpy as np
import csv, glob, math, datetime, os, time
import pandas as pd
from pathlib import Path

import DGRutils as utils

def cleanup_contigs(contig_file, min_length=2000):
	contig_path = Path(contig_file)
	contig_base_name = contig_path.stem
	contig_parent = str(contig_path.parent)
	output_file = contig_parent + '/' + contig_base_name + '_cleanup.fa'

	with open(output_file, 'w') as f:
		records = SeqIO.parse(contig_file, 'fasta')
		i = 0
		for item in records:
			contig_len = len(str(item.seq))
			if contig_len >= min_length:
				f.write('>Contig%i len:%i\n%s\n' % (i, contig_len, str(item.seq)))
				i += 1

from Bio import SeqIO
import csv, time, subprocess, os
import pandas as pd
from pathlib import Path
from rich.progress import Progress
from rich import print
import pyrodigal_gv
import pyhmmer
import concurrent.futures as cf
import gzip
import argparse

import DGRutils as utils
from dataclasses import dataclass


if __name__ == '__main__':
	# Argument parsing
	parser = argparse.ArgumentParser(description='Find DGRs in metagenomic assembly data')
	parser.add_argument('-i',dest='assembly',required=True,
	                    help='Path to assembly file in gzipped fasta format')
	parser.add_argument('-o',dest='outfile',required=True,
	                    help='Path to output file in .csv format')
	parser.add_argument('-hmm',dest='hmm',required=True,
	                    help='Path to HMM profile for RT searching')
	parser.add_argument('-c',dest='contigs_out',required=True,
	                    help='Path to output dgr-containing contigs to (fasta format)')
	parser.add_argument('-cpu',dest='threads',default=1,type=int,
						help='Number of cpus to use (default 1)')
	args = parser.parse_args()



def create_RT_hmmprofile(hmmfile=None):
	hmmfolder = 'hmmprofiles'
	#Look for hmmfile, if not there, then try default
	if hmmfile != None:
		if Path(hmmfile).exists():
			hmm = hmmfile
		else:
			if Path('hmmprofiles/DGR_RTs_MSA.fa').exists():
				hmm = 'hmmprofiles/DGR_RTs_MSA.fa'
			else:
				raise OSError(f"HMM Profile file: {hmmfile} does not exist.  Cannot find default file either.")
	else:
		if Path('hmmprofiles/DGR_RTs_MSA.fa').exists():
			hmm = 'hmmprofiles/DGR_RTs_MSA.fa'
		else:
			raise OSError('Default HMM profile hmmprofiles/DGR_RTs_MSA.fa does not exist.  Please download from github')

	hmm_files = ['hmmprofiles/DGR_RTs_MSA.h3m', 'hmmprofiles/DGR_RTs_MSA.h3i', 'hmmprofiles/DGR_RTs_MSA.h3f', 'hmmprofiles/DGR_RTs_MSA.h3p', 'hmmprofiles/DGR_RTs_MSA']
	for hf in hmm_files:
		if Path(hf).exists():
			subprocess.call(f"rm {hf}", shell=True)
	subprocess.call(f"hmmbuild hmmprofiles/DGR_RTs_MSA {hmm}", shell=True)
	subprocess.call('hmmpress hmmprofiles/DGR_RTs_MSA', shell=True)

def search_for_DGRs(contigs, outputfile, filtered_contig_file, hmm, threads):

	start_time = time.time()
	#Get the name of the rawdata file and remove the file extension
	rawdata_name = contigs.rsplit('/')[-1].rsplit('.f', 1)[0]
	#Create temp_directory in output_folder to store information
	temp_folder = f"temp/{rawdata_name}"
	#subprocess.call('rm -rf %s' % (temp_folder))
	#time.sleep(3)
	Path(temp_folder).mkdir(parents=True, exist_ok=True)

	print(f"Searching for DGRs in [bold red]{rawdata_name}[/bold red] using {threads} threads")

	#Find all ORFs in the rawdata and convert to protein sequences for input into hmmsearch using pyrodigal-gv in metagenome mode
	@dataclass
	class Protein:
		id: str
		contig: str
		start: int
		end: int
		nseq: str
		aseq: str

	### Load in input sequences
	with gzip.open(contigs, "rt") as handle:
		all_sequences = list(SeqIO.parse(handle, "fasta"))
		if len(all_sequences) < 1000:
			print(f"Detected less than 1000 contigs, assuming this is MAG\nWill not remove any contigs < 1000bps")
			sequences = all_sequences
		else:
			print(f"Detected assembly with {len(all_sequences)} contigs... assume metagenomic assembly\nWill remove any sequences < 1000 bps as unlikely to harbor DGR")
			sequences = [sequence for sequence in all_sequences if len(sequence.seq) >= 1000]
	### Calculate basic stats
	total_bp = 0
	contig_num = len(sequences)
	if len(sequences) > 1:
		for seq in sequences:
			total_bp += len(seq.seq)
		mbp = total_bp/1000000
		l50 = total_bp/2
		n50 = []
		running_total = 0
		if len(sequences) > 1:
			for seq in sequences:
				if running_total < l50:
					running_total += len(seq.seq)
				else:
					n50.append(len(seq.seq))
					running_total += len(seq.seq)
		else:
			for seq in sequences:
				n50 = len(seq)
		avg_contig_len = total_bp/contig_num
		print(f"After removing all contigs shorter than 1,000 bp (unlikely to harbor DGR)\nFound {contig_num} contigs\nAverage Length: {avg_contig_len} bp\nN50: {max(n50)} bp\nTotal bp: {mbp} Mb")
	### Find ORFs using multithreaded pyrodigal-gv and status bar
	orf_finder = pyrodigal_gv.ViralGeneFinder(meta=True, mask=True)
	protein_seqs = f"{temp_folder}/{rawdata_name}-protein_seqs.faa"
	def predict_genes(seq):
			return (seq.id, len(seq.seq), orf_finder.find_genes(bytes(seq.seq)))
	proteins = []
	with Progress() as progress:
		task = progress.add_task("Predicting ORFs with pyrodigal-gv...", total=total_bp)
		with cf.ThreadPoolExecutor(max_workers=threads) as tpe:
			for seq_id, length, preds in tpe.map(predict_genes, sequences):
				for pred in preds:
					nseq = pred.sequence()
					aseq = pred.translate()
					start = pred.begin
					end = pred.end
					protein_id = f"{seq_id}_{start}_{end}"
					proteins.append(Protein(protein_id, seq_id, start, end, nseq, aseq))
				progress.update(task, advance=length)
	with open(protein_seqs, "w") as prot_file:
		for protein in proteins:
			prot_file.write(f">{protein.id}\n{protein.aseq}\n")
	
	@dataclass
	class Result:
		query: str
		cog: str
		bitscore: int
	with pyhmmer.easel.SequenceFile(protein_seqs, digital=True, format="fasta", alphabet=pyhmmer.easel.Alphabet.amino()) as seq_file:
		seq_list = list(seq_file)

	orf_total = len(proteins)
	print(f"Found {orf_total} ORFs\nNow searching for putative RTs...")

	results = []
	with pyhmmer.plan7.HMMFile(hmm) as hmm_file:
		for hits in pyhmmer.hmmscan(seq_list, hmm_file, cpus=threads):
			cog = hits.query_name.decode()
			for hit in hits:
				if hit.included and hit.score > 20:
					results.append(Result(hit.name.decode(), cog, hit.score))

	prodigal_info = {}
	for protein in proteins:
		prodigal_info[protein.id] = [protein.contig, protein.start, protein.end]

	RT_hits = []
	for result in results:
		contig, start, end = prodigal_info[result.cog]
		aa_length = abs(end - start) / 3
		if aa_length >= 200 and aa_length <= 550:
			RT_hits.append([contig, start, end, result.cog])

	### Make list of contigs that contain RT hit
	contig_hits = []
	for RT_hit in RT_hits:
		contig_hits.append(RT_hit[0])
	
	#Extract contigs that contain RT hits and remove large formatted contigs file
	rt_hit_length = len(RT_hits)
	print(f"{rt_hit_length} putative RTs identified")
	
	print("Extracting contigs containing putative RTs for faster VR/TR searching...")
	filtered_contigs = "%s/%s-filtered_contigs.fa" % (temp_folder, rawdata_name)
	with gzip.open(contigs, "rt") as handle:
		records = list(SeqIO.parse(handle, "fasta"))
	sequences = (record for record in records if record.name in contig_hits)
	with open(filtered_contigs, "w") as output_handle:
		SeqIO.write(sequences, output_handle, "fasta")

	#Extract 20kb up and downstream of each RT in order to analyze for VR/TR pairs
	def blast(RT_hit):
		# Define extraction window
		extract_window = 20000
		# Define contig containing potential RT hit
		RT_contig = RT_hit[0]
		# Define start of extraction window from RT start position, cannot be less than 0
		RT_area_start = int(RT_hit[1]) - extract_window
		if RT_area_start < 0:
			RT_area_start = 0
		# Open formatted_contigs file as SeqIO object
		records = SeqIO.parse(filtered_contigs, 'fasta')
		contig = (record for record in records if record.name == RT_contig)
		for con in contig:
			seq = str(con.seq)
		#Define end of extraction window from length of contig after opening, cannot be past end of sequence
		RT_area_end = int(RT_hit[2]) + extract_window
		if RT_area_end > len(seq)-1:
			RT_area_end = len(seq)-1
		# Define length of sequence and number of chunks to be written to sliding window file
		window_size = 200
		step_size = 50
		RT_area_length = RT_area_end - RT_area_start
		num_chunks = ((RT_area_length-window_size)//step_size)+1
		# Create sliding window file with known RT_area parameters as defined above
		temp_out_base = f"{temp_folder}/{RT_contig}-{RT_hit[1]}-{RT_hit[2]}"
		sliding_window_file = f"{temp_out_base}-sliding_window.fa"
		with open(sliding_window_file, 'w') as w:
			for x in range(RT_area_start, RT_area_start+(num_chunks*step_size), step_size):
				y = x + window_size
				if y > RT_area_end:
					dna = seq[x:RT_area_end]
					y = RT_area_end
				else:
					dna = seq[x:x+window_size]
				w.write(f">{con.name}_{x}_to_{y}\n{dna}\n")
		w.close()

		### Make blastdb of 200bp segments surrounding RT
		sliding_window_blastdb = f"{temp_out_base}-sliding_window_blast_db"
		cmd = "makeblastdb -in %s -dbtype nucl -out %s" % (sliding_window_file, sliding_window_blastdb)
		subprocess.call(cmd, shell=True, stdout=open(os.devnull, 'wb'))
		
		### Define blast parameters
		blast_sliding_window_hits = f"{temp_out_base}-sliding_window-blast_output.txt"
		output_options = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'sseq', 'qseq']
		blast_out_str = ' '.join(output_options)
		word_size = 8
		reward = 1
		penalty = -1
		evalue = 1e-5
		gapopen = 6
		gapextend = 6
		perc_identity = 50
		num_threads = 1

		# Launch blastn command
		cmd = "blastn -query %s -db %s -out %s -outfmt '6 %s' -word_size %s -reward %s -penalty %s -evalue %s -gapopen %s -gapextend %s -perc_identity %s -num_threads %s" % (sliding_window_file, sliding_window_blastdb, blast_sliding_window_hits, blast_out_str, word_size, reward, penalty, evalue, gapopen, gapextend, perc_identity, num_threads)
		subprocess.call(cmd, shell=True)
	
	with Progress() as progress:
		task_len = len(RT_hits)
		task = progress.add_task("Searching for VR-TR pairs...", total=task_len)
		with cf.ThreadPoolExecutor(max_workers=threads) as executor:
			for RT_hit in executor.map(blast, RT_hits):
				progress.update(task, advance=1)

	print(f"Determining true VR/TR pairs for sample: {rawdata_name}")

	# Create list to store all raw_pair dicts
	filtered_pairs = []
	@dataclass
	class Pair:
		RT_hit: list
		vr_start: int
		vr_end: int
		vr_seq: str
		tr_start: int
		tr_end: int
		tr_seq: str

	def parse_blast_hits(RT_hit):
		raw_pairs = []
		#Define output options indices
		output_options = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'sseq', 'qseq']
		mismatch = output_options.index('mismatch')
		query_seq = output_options.index('qseq')
		sbjct_seq = output_options.index('sseq')
		query_name = output_options.index('qseqid')
		sbjct_name = output_options.index('sseqid')
		qstart = output_options.index('qstart')
		qend = output_options.index('qend')
		sstart = output_options.index('sstart')
		send = output_options.index('send')

		#Define blast output file
		temp_out_base = f"{temp_folder}/{RT_hit[0]}-{RT_hit[1]}-{RT_hit[2]}"
		blast_sliding_window_hits = f"{temp_out_base}-sliding_window-blast_output.txt"

		extract_window = 20000
		# Define start of extraction window from RT start position, cannot be less than 0
		RT_area_start = int(RT_hit[1]) - extract_window
		if RT_area_start < 0:
			RT_area_start = 0

		num_hits = 0
		with open(blast_sliding_window_hits, 'r') as blast_hits:
			blast_reader = csv.reader(blast_hits, delimiter='\t')
			for line in blast_reader:
				if int(line[mismatch]) > 0:
					pair = utils.determine_VR_TR_pair(line[query_seq], line[sbjct_seq])
					if pair != False:
						num_hits += 1
						seq_length = len(line[query_seq])
						if pair[0] == line[query_seq]:
							VR_name = line[query_name]
							TR_name = line[sbjct_name]
							vstart_pos = int(line[qstart])
							if int(line[qend]) < vstart_pos:
								vstart_pos = int(line[qend])
							tstart_pos = int(line[sstart])
							if int(line[send]) < tstart_pos:
								tstart_pos = int(line[send])
						else:
							VR_name = line[sbjct_name]
							TR_name = line[query_name]
							tstart_pos = int(line[qstart])
							if int(line[qend]) < tstart_pos:
								tstart_pos = int(line[qend])
							vstart_pos = int(line[sstart])
							if int(line[send]) < vstart_pos:
								vstart_pos = int(line[send])

						vstart = int(VR_name.split('_')[-3]) + vstart_pos-1
						tstart = int(TR_name.split('_')[-3]) + tstart_pos-1
					
						raw_pair_list = [vstart, vstart+seq_length, tstart, tstart+seq_length]
						raw_pairs.append(raw_pair_list)

		### Given list of putative TR-VR pairs, merge overlapping VR and TR intervals only if both the VR and TR intervals overlap, should keep unique VR-TR combos
		if num_hits>0:

			sorted_by_tstart = sorted(raw_pairs, key=lambda x: (x[2], x[0]))
			index=0
			merged = []
			for i in range(1, len(sorted_by_tstart)):
				if sorted_by_tstart[index][3] >= sorted_by_tstart[i][2] and sorted_by_tstart[index][1] >= sorted_by_tstart[i][0]:
					sorted_by_tstart[index][3] = max(sorted_by_tstart[index][3], sorted_by_tstart[i][3])
					sorted_by_tstart[index][1] = max(sorted_by_tstart[index][1], sorted_by_tstart[i][1])
				else:
					index = index + 1
					sorted_by_tstart[index] = sorted_by_tstart[i]
			for i in range(index+1):
				merged.append(sorted_by_tstart[i])
			sequences = SeqIO.parse(filtered_contigs, 'fasta')
			for sequence in sequences:
				if sequence.id == RT_hit[0]:
					for line in merged:
						vr_seq = sequence[line[0]:line[1]].seq
						tr_seq = sequence[line[2]:line[3]].seq
						filtered_pair = Pair(
							RT_hit=RT_hit, 
							vr_start=int(line[0]), 
							vr_end=int(line[1]), 
							vr_seq=vr_seq, 
							tr_start=int(line[2]), 
							tr_end=int(line[3]), 
							tr_seq=tr_seq)
						filtered_pairs.append(filtered_pair)
	
	with cf.ThreadPoolExecutor(max_workers=threads) as executor:
		executor.map(parse_blast_hits, RT_hits)

	# Determine if VR exists within an ORF and write outputs to Final dataclass
	@dataclass
	class Final:
		metagenome: str
		dgr_id: str
		contig: str
		rt_start: int
		rt_end: int
		rt_seq: str
		rt_prot_seq: str
		vr_start: int
		vr_end: int
		vr_seq: str
		tr_start: int
		tr_end: int
		tr_seq: str
		tg_start: int
		tg_end: int
		tg_seq: str
		vp_id: str
		vp_seq: str
	
	out = []
	contig_list = []
	def write_output(dgr):
		found = False
		#Check to see if VR exists within ORF from protein dataclass
		for protein in proteins:
			if protein.contig == dgr.RT_hit[0]:
				if (protein.start <= dgr.vr_start <= protein.end) or (protein.start <= dgr.vr_end <= protein.end):
					found = True
					vp_seq_id = protein.id
					vp_start = protein.start
					vp_end = protein.end
					rt_seq_id = dgr.RT_hit[3]
					### Extract gene sequences for RT and VP from protein dataclass
					for protein in proteins:
						if protein.id == rt_seq_id:
							rt_gene_seq = str(protein.nseq)
							rt_prot_seq = str(protein.aseq)
						if protein.id == vp_seq_id:
							vp_gene_seq = str(protein.nseq)
							vp_prot_seq = str(protein.aseq)

		if found == True:
			dgr_id = f"{rawdata_name}-{dgr.RT_hit[0]}-{dgr.RT_hit[1]}-{dgr.RT_hit[2]}"
			vp_id = f"{rawdata_name}-{dgr.RT_hit[0]}-{vp_start}-{vp_end}"
			final_class = Final(
				metagenome=rawdata_name, 
				dgr_id=dgr_id, 
				contig=dgr.RT_hit[0], 
				rt_start=dgr.RT_hit[1], 
				rt_end=dgr.RT_hit[2], 
				rt_seq=rt_gene_seq, 
				rt_prot_seq=rt_prot_seq, 
				vr_start=dgr.vr_start, 
				vr_end=dgr.vr_end, 
				vr_seq=dgr.vr_seq, 
				tr_start=dgr.tr_start, 
				tr_end=dgr.tr_end, 
				tr_seq=dgr.tr_seq, 
				tg_start=vp_start, 
				tg_end=vp_end, 
				tg_seq=vp_gene_seq, 
				vp_id=vp_id, 
				vp_seq=vp_prot_seq)
		out_row = pd.DataFrame([final_class])
		out.append(out_row)
		contig_list.append(final_class.contig)

	with cf.ThreadPoolExecutor(max_workers=threads) as executor:
		executor.map(write_output, filtered_pairs)
	
	output_path = f"{outputfile}"
	if len(out) >=1:
		final = pd.concat(out)
		final.to_csv(output_path, mode='w', index=False)
	if len(out) == 0:
		subprocess.call(f"touch {output_path}", shell=True)
	
	dgr_contig_list = list(set(contig_list))
	contigs = SeqIO.parse(filtered_contigs, 'fasta')
	dgr_contigs = [record for record in contigs if record.name in dgr_contig_list]
	to_add = f"{rawdata_name}-"
	with open(filtered_contig_file, 'w') as fc_writer:
		for record in dgr_contigs:
			record.id = (to_add + record.description)
			record.description = record.id
			SeqIO.write(record, fc_writer, 'fasta')
	
	subprocess.call(f"rm -r {temp_folder}/*", shell=True)
	subprocess.call(f"rm -r {temp_folder}", shell=True)
	
	end_time = time.time()
	total_time = end_time-start_time
	total_min = total_time/60
	print(f"{len(out)} VR-TR pairs found in {total_min} minutes")

search_for_DGRs(args.assembly, args.outfile, args.contigs_out, args.hmm, args.threads)

	
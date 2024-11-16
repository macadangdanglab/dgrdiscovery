from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline

import numpy as np
import csv, glob, math, datetime, os, time
import pandas as pd
from pathlib import Path

import DGRutils as utils

def determine_DGR_activity_with_aligners(rawdatafile, reference_genomefile, VR_file, TR_file, output_folder):
	#Create temp_directory in output_folder to store information
	temp_folder = '%s/temp' % (output_folder)
	os.system('rm -rf %s' % (temp_folder))
	time.sleep(3)
	Path(temp_folder).mkdir(parents=True, exist_ok=True)

	#Keep track of files that need to be deleted at the end of the analysis
	temp_files = []
	temp_blast_db = []

	#Format the reference genome file for easier searching later
	#formatted_ref_genomefile = '%s/formatted_ref_genome.fasta' % (temp_folder)
	formatted_ref_genomefile = reference_genome_name
	temp_files.append(formatted_ref_genomefile)
	utils.format_ref_genome_file(reference_genomefile, formatted_ref_genomefile)

	VR_start, VR_end, VR_contig_num, VR_seq = utils.extract_sequence(VR_file, formatted_ref_genomefile)
	TR_start, TR_end, TR_contig_num, TR_seq = utils.extract_sequence(TR_file, formatted_ref_genomefile)

	if len(VR_seq) != len(TR_seq):
		raise ValueError('VR and TR lengths are not equal')

	#Get the name of the reference genome file and remove the file extension
	reference_genome_name = '.'.join(reference_genomefile.split('/')[-1].split('.')[:-1])

	#Get the name of the rawdata file and remove the file extension
	rawdata_name = '.'.join(rawdatafile.split('/')[-1].split('.')[:-1])

	unique_name = '%s-%s' % (reference_genome_name, rawdata_name)

	#Extract 20kb around VR region to create reference genome for alignment
	#Define VR area +/- 100 bp in order to map to VR region
	vr_20kb_filename = '%s/VR-20kb.fasta' % (temp_folder)
	temp_files.append(vr_20kb_filename)
	with open(vr_20kb_filename, 'w') as vr_writer:
		with open(formatted_ref_genomefile, 'r') as ref_genome:
			rg_parser = SeqIO.parse(ref_genome, 'fasta')
			for contig in rg_parser:
				if contig.name == 'Contig%i' % (VR_contig_num):
					vr_area_start = VR_start - 10000
					if vr_area_start < 0:
						vr_area_start = 0
					vr_area_end = VR_end + 10000
					if vr_area_end > len(contig.seq):
						vr_area_end = len(contig.seq)
					vr_writer.write('>VR_area_20kb\n%s' % (str(contig.seq[vr_area_start:vr_area_end])))
					break

	
	#Create reference genome with BWA
	print('Creating index with BWA')
	os.system('bwa index %s' % (vr_20kb_filename))

	#Align to the 20kb region
	align_20kb_sam_file = '%s/%s-20kb_align.sam' % (temp_folder, unique_name)
	temp_files.append(align_20kb_sam_file)
	os.system('bwa mem -v 1 %s %s > %s' % (vr_20kb_filename, rawdatafile, align_20kb_sam_file))
	
	#Convert sam to bam and remove sam
	print('Converting from sam to bam')
	align_20kb_bam_file = '%s/%s-20kb_align.bam' % (temp_folder, unique_name)
	temp_files.append(align_20kb_bam_file)
	os.system('samtools view --threads 2 -b -o %s %s' % (align_20kb_bam_file, align_20kb_sam_file))
	#os.system('rm %s/%s-20kb_align.sam' % (temp_folder, unique_name))
	

	#Sort the bam file
	print('Sorting bam file')
	align_20kb_sorted_bam_file = '%s/%s-20kb_align_sorted.bam' % (temp_folder, unique_name)
	temp_files.append(align_20kb_sorted_bam_file)
	os.system('samtools sort -m 7G --threads 2 -o %s %s' % (align_20kb_sorted_bam_file, align_20kb_bam_file))

	#Filter the bam reads
	print('Filtering bam file')
	#os.system('java -jar /Users/brmac/tools/fgbio/fgbio-1.1.0.jar FilterBam -m 50 -i %s/%s-20kb_align.bam -o %s/%s-20kb_align_filtered.bam' % (temp_folder,unique_name,temp_folder,unique_name))
	align_20kb_filtered_bam_file = '%s/%s-20kb_align_filtered.bam' % (temp_folder, unique_name)
	temp_files.append(align_20kb_filtered_bam_file)
	os.system('java -jar /home/brm4/fgbio/fgbio-1.1.0.jar FilterBam -m 50 -i %s -o %s' % (align_20kb_sorted_bam_file, align_20kb_filtered_bam_file))

	#Convert back to sam
	print('Converting filtered bam back to sam')
	align_20kb_filtered_sam_file = '%s/%s-20kb_align_filtered.sam' % (temp_folder, unique_name)
	temp_files.append(align_20kb_filtered_sam_file)
	os.system('samtools view --threads 2 -o %s %s' % (align_20kb_filtered_sam_file, align_20kb_filtered_bam_file))

	#Build index with BBMap
	os.system('bbmap.sh ref=%s' % (vr_20kb_filename))

	#Index using BBMap
	global_align_file = '%s/%s-20kb_global_align.sam' % (output_folder, unique_name)
	#temp_files.append(global_align_file)
	os.system('bbmap.sh vslow minid=0 indelfilter=2 inslenfilter=3 dellenfilter=3 in=%s out=%s' % (align_20kb_filtered_sam_file, global_align_file))
	
	#Extract all reads that mapped to the VR region, which starts at 10,000 (so extract +/- 100 bp from 10,000)
	vr_region_aligners = '%s/VR_BWA_BBmap-%s.fasta' % (output_folder, unique_name)
	with open('%s' % (global_align_file)) as f:
		with open(vr_region_aligners, 'w') as vr_writer:
			reader = csv.reader(f, delimiter='\t')
			for line in reader:
				if line[0][0] != '@':
					read_start = int(line[3])
					if read_start > 9900 and read_start < (10000 + len(VR_seq) + 100):
						vr_writer.write('>%s\n%s\n' % (line[0], line[9]))

	
	#Align VRs usingn blast
	aligned_VR_sequences = '%s/aligned_VR_sequences_BWA_BBMap-%s.fasta' % (output_folder, unique_name)
	vr_blast_database = '%s/vrblastdb' % (temp_folder)
	vr_blast_output = '%s/vr_aligning_blastout.xml' % (temp_folder)

	temp_blast_db.append(vr_blast_database)
	temp_files.append(vr_blast_output)

	cline = NcbimakeblastdbCommandline(dbtype='nucl', input_file=VR_file, out=vr_blast_database)
	cline()

	cline = NcbiblastnCommandline(out=vr_blast_output, db=vr_blast_database, query=vr_region_aligners, outfmt=5, word_size=8, reward=1, penalty=-1, evalue=1e-4, gapopen=2, gapextend=1, perc_identity=50)
	cline()

	amiss, tmiss = 0, 0
	for i in range(len(VR_seq)):
		if VR_seq[i] != TR_seq[i]:
			if TR_seq[i] == 'A':
				amiss += 1
			elif TR_seq[i] == 'T':
				tmiss += 1
	reverse = False
	if tmiss > amiss:
		reverse = True

	with open(vr_blast_output, 'r') as f:
		results = NCBIXML.parse(f)
		with open(aligned_VR_sequences, 'w') as aligned_VR_writer:
			vr_oriented = VR_seq
			if reverse:
				vr_oriented = str(Seq(VR_seq).reverse_complement())
			aligned_VR_writer.write('>%s\n%s\n' % ('VR', vr_oriented))
			try:
				for result in results:
					if len(result.alignments) > 0:
						lowest = 1
						best = [0, 0]
						for anum, alignment in enumerate(result.alignments):
							for hnum, hsp in enumerate(alignment.hsps):
								if hsp.expect < lowest:
									lowest = hsp.expect
									best = [anum, hnum]
						alignment = result.alignments[best[0]]
						hsp = alignment.hsps[best[1]]
						seq_out = ''
						if hsp.sbjct_start < hsp.sbjct_end:
							start = hsp.sbjct_start - 1
						else:
							start = hsp.sbjct_end - 1
						i = 0
						for i in range(start):
							seq_out += '-'
						if hsp.sbjct_start < hsp.sbjct_end:
							seq_out += hsp.query
						else:
							seq_out += str(Seq(hsp.query).reverse_complement())
						end = len(VR_seq) - start - len(hsp.query)
						i = 0
						for i in range(end):
							seq_out += '-'
						if reverse:
							seq_out = str(Seq(seq_out).reverse_complement())
						aligned_VR_writer.write('>%s\n%s\n' % (result.query, seq_out))
			except:
				print('XML file empty')
				utils.cleanup(temp_files, temp_blast_db)

	utils.cleanup(temp_files, temp_blast_db)

def process_determine_DGR_activity(rawdatafile, reference_genomefile, VR_file, TR_file, output_folder, rawdatafile2, vt_folder):
	if vt_folder:
		if VR_file[-1] != '/':
			VR_file += '/'
		if TR_file[-1] != '/':
			TR_file += '/'
		VRs = sorted(glob.glob('%s*-VR.f*' % (VR_file)))
		TRs = sorted(glob.glob('%s*-TR.f*' % (TR_file)))
		if len(VRs) == 0:
			print('No DGRs found for this genome. Returning')
			return
	else:
		VRs = [VR_file]
		TRs = [TR_file]
	determine_DGR_activity_from_metagenome(rawdatafile, reference_genomefile, VRs, TRs, output_folder, rawdatafile2)
	temp_folder = '%s/temp' % (output_folder)
	os.system('rm -rf %s' % (temp_folder))


def determine_DGR_activity_from_metagenome(rawdatafile, reference_genomefile, VRs, TRs, output_folder, rawdatafile2):
	#Create temp_directory in output_folder to store information
	temp_folder = '%s/temp' % (output_folder)
	Path(temp_folder).mkdir(parents=True, exist_ok=True)

	#Keep track of files that need to be deleted at the end of the analysis
	temp_files = []
	temp_blast_db = []

	#Format the reference genome file for easier searching later
	formatted_ref_genomefile = '%s/formatted_ref_genome.fasta' % (temp_folder)
	temp_files.append(formatted_ref_genomefile)
	utils.format_ref_genome_file(reference_genomefile, formatted_ref_genomefile)

	#Get the name of the reference genome file and remove the file extension
	reference_genome_name = '.'.join(reference_genomefile.split('/')[-1].split('.')[:-1])
	reference_genome_name = reference_genome_name.replace('_contigs', '')

	#Get the name of the rawdata file and remove the file extension
	rawdata_name = rawdatafile.split('/')[-1].split('.')[0]

	#Create error list for later debugging:
	errors = []
	error_file = '%s/DGR_analysis_%s_errors.txt' % (output_folder, rawdata_name)

	#Determine if rawdatafile is in .gz format
	if rawdatafile.split('.')[-1] == 'gz':
		print('Uncompressing raw data')
		os.system('cp %s %s/%s.fastq.gz' % (rawdatafile, temp_folder, rawdata_name))
		os.system('unpigz -p 4 %s/%s.fastq.gz' % (temp_folder, rawdata_name))
		temp_files.append('%s/%s.fastq' % (temp_folder, rawdata_name))
		if rawdatafile2 is not None:
			os.system('cp %s %s/%s_2.fastq.gz' % (rawdatafile, temp_folder, rawdata_name))
			os.system('unpigz -p 4 %s/%s_2.fastq.gz' % (temp_folder, rawdata_name))
			os.system('cat %s/%s_2.fastq >> %s/%s.fastq' % (temp_folder, rawdata_name, temp_folder, rawdata_name))
			os.system('rm %s/%s_2.fastq' % (temp_folder, rawdata_name))
		rawdatafile = '%s/%s.fastq' % (temp_folder, rawdata_name)


	#If the rawdata file is in fastq, it needs to be converted to fasta for blast, delete that file at the end
	#Otherwise just use the fastafile supplied and it does not need ot be deleted
	if utils.is_fastq_file(rawdatafile):
		print('Converting rawdata from FASTQ to FASTA')
		rawdata_fasta = '%s/%s.fa' % (temp_folder, rawdata_name)
		temp_files.append(rawdata_fasta)
		with open(rawdata_fasta, 'w') as f:
			data = SeqIO.parse(rawdatafile, 'fastq')
			i = 1
			for seq in data:
				f.write('>Sequence%i\n%s\n' % (i, str(seq.seq)))
				i += 1
	else:
		rawdata_fasta = rawdatafile

	
	vr_100_filename = '%s/VR-100bp.fasta' % (temp_folder)
	temp_files.append(vr_100_filename)
	#Create empty file
	with open(vr_100_filename, 'w') as vr_100:
		pass

	for VR_file, TR_file in zip(VRs, TRs):
		if Path(VR_file).stat().st_size > 0:
			VR_start, VR_end, VR_contig_num, VR_seq = utils.extract_sequence(VR_file, formatted_ref_genomefile)
			TR_start, TR_end, TR_contig_num, TR_seq = utils.extract_sequence(TR_file, formatted_ref_genomefile)

			if len(VR_seq) != len(TR_seq):
				raise ValueError('VR and TR lengths are not equal')

			unique_name = '%s-Contig%s_%s_%s' % (reference_genome_name, VR_contig_num, VR_start, VR_end)

			print('Creating VR area files for %s' % (unique_name))

			#Define VR area +/- 100 bp in order to map to VR region
			with open(vr_100_filename, 'a') as vr_100:
				with open(formatted_ref_genomefile, 'r') as ref_genome:
					rg_parser = SeqIO.parse(ref_genome, 'fasta')
					for contig in rg_parser:
						if contig.name == 'Contig%i' % (VR_contig_num):
							vr_area_start = VR_start - 100
							if vr_area_start < 0:
								vr_area_start = 0
							vr_area_end = VR_end + 100
							if vr_area_end > len(contig.seq):
								vr_area_end = len(contig.seq)
							vr_100.write('>VR_area_100-%s\n%s\n' % (unique_name, str(contig.seq[vr_area_start:vr_area_end])))
							break
		else:
			errors.append('%s in %s did not contain a VR' % (VR_file, rawdata_name))

	#Store all sequences from raw data that may potentially match to the VR area +/- 100 bp
	sequences_matched_to_vr_100_area = '%s/%s-seqs_match_VR100.fasta' % (temp_folder, rawdata_name)
	vr_100_blast_database = '%s/vr100_blastdb' % (temp_folder)
	blastoutput_from_vr100_blast = '%s/blastoutput_vr100.txt' % (temp_folder)
	
	#Delete these files at the end during cleanup
	temp_files.append(sequences_matched_to_vr_100_area)
	temp_blast_db.append(vr_100_blast_database)
	temp_files.append(blastoutput_from_vr100_blast)

	#Create blast database with the VR+/-100 area
	cline = NcbimakeblastdbCommandline(dbtype='nucl', input_file=vr_100_filename, out=vr_100_blast_database)
	cline()

	print("Finding potential raw data sequences that match to the surround VR area")

	#Setup blast output options
	#output_options = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'sseq', 'qseq']
	output_options = ['qseqid', 'length', 'qlen']
	blast_out_str = ' '.join(output_options)

	#Blast all rawdata to the VR+/-100 blast database
	cline = NcbiblastnCommandline(out=blastoutput_from_vr100_blast, db=vr_100_blast_database, query=rawdata_fasta, outfmt='6 %s' % (blast_out_str), word_size=8, reward=1, penalty=-1, evalue=1e-4, gapopen=6, gapextend=6, perc_identity=80, task='blastn', dust='no')
	cline()	

	
	#Find sequences that had at least 50% of their sequence align (gets rid of partially aligned sequences)
	potential_seqs = []
	with open(blastoutput_from_vr100_blast, 'r') as f:
		reader = csv.reader(f, delimiter='\t')
		for query_id, align_len, query_len in reader:
			#Check to see if at least 50% of the rawdata sequence aligned
			if int(align_len) >= (0.8 * int(query_len)):
				potential_seqs.append(query_id)
	
	if len(potential_seqs) == 0:
		#utils.cleanup(temp_files, temp_blast_db)
		print('There were no sequences that matched to VR')
		return


	#Now read in the rawdata file and create new data file is only the sequences that aligned to VR area +/- 100 bp
	#Output: sequences_matched_to_vr_100_area now contains all rawdata reads that potentially match to VR
	with open(sequences_matched_to_vr_100_area, 'w') as vr_100_writer:
		data = SeqIO.parse(rawdata_fasta, 'fasta')
		for sequence in data:
			if sequence.name == potential_seqs[0]:
				vr_100_writer.write('>%s\n%s\n' % (sequence.name, str(sequence.seq)))
				if len(potential_seqs) > 1:
					potential_seqs.pop(0)
				else:
					break

	

	print("Determining which candidate sequences match best to VR")

	#Go through the rawdata reads that potentially map to VR and determine if there is a better match somewhere else by using the entire ref genome
	#First create a new blast database with the entire reference genome
	entire_ref_genome_blast_database = '%s/ref_gen_blastdb' % (temp_folder)
	temp_blast_db.append(entire_ref_genome_blast_database)

	cline = NcbimakeblastdbCommandline(dbtype='nucl', input_file=formatted_ref_genomefile, out=entire_ref_genome_blast_database)
	time.sleep(3)
	cline()

	#Next create files for output, need to parse XML files to determine is best match is somewhere else
	best_alignments_blastoutput = '%s/best_alignments_blastoutput.xml' % (temp_folder)
	temp_files.append(best_alignments_blastoutput)

	#Blast the candidate reads to the entire genome to find alignment, using a more stringent alignment.
	#If read doesn't match anywhere else or matches best to VR region, then will use the read for downstream analysis
	#Input file: sequences_matched_to_vr_100_area
	#To do: find the optimum search settings for more stringency
	cline = NcbiblastnCommandline(out=best_alignments_blastoutput, db=entire_ref_genome_blast_database, query=sequences_matched_to_vr_100_area, outfmt=5, word_size=20, reward=1, penalty=-2, evalue=1e-4, gapopen=6, gapextend=2, perc_identity=80)
	cline()

	#Now parse the XML file and find the best match (or no match)
	#Output includes all transcripts that match to VR (but some may also match to TR and thus need to be filtered)
	vr_tr_transcripts = '%s/vr_tr_transcripts.fasta' % (temp_folder)
	reads_that_matched_better_somewhere_else = '%s/reads_that_matched_better_somewhere_else.fasta' % (temp_folder)
	temp_files.append(vr_tr_transcripts)
	#temp_files.append(reads_that_matched_better_somewhere_else)

	for VR_file, TR_file in zip(VRs, TRs):
		if Path(VR_file).stat().st_size > 0:
			VR_start, VR_end, VR_contig_num, VR_seq = utils.extract_sequence(VR_file, formatted_ref_genomefile)
			TR_start, TR_end, TR_contig_num, TR_seq = utils.extract_sequence(TR_file, formatted_ref_genomefile)

			unique_name = '%s-Contig%s_%s_%s' % (reference_genome_name, VR_contig_num, VR_start, VR_end)

			print('Processing blast output for %s' % (unique_name))

			vr_tr_names = []
			with open(best_alignments_blastoutput, 'r') as blastoutput:
				parser = NCBIXML.parse(blastoutput)
				for result in parser:
					#For each sequence, find the best hit
					if len(result.alignments) > 0:
						lowest = 1
						best = [0, 0]
						for anum, alignment in enumerate(result.alignments):
							for hnum, hsp in enumerate(alignment.hsps):
								if hsp.expect < lowest:
									lowest = hsp.expect
									best = [anum, hnum]
						alignment = result.alignments[best[0]]
						hsp = alignment.hsps[best[1]]
						if alignment.hit_def == 'Contig%i' % (VR_contig_num):
							within = False
							if hsp.sbjct_start >= VR_start and hsp.sbjct_start <= VR_end:
								within = True
							if hsp.sbjct_end >= VR_start and hsp.sbjct_end <= VR_end:
								within = True
							if hsp.sbjct_start <= VR_start and hsp.sbjct_end >= VR_end:
								within = True
							if within:
								vr_tr_names.append(result.query)

			if len(vr_tr_names) == 0:
				#utils.cleanup(temp_files, temp_blast_db)
				print('There were no sequences that matched to VR after searching the entire genome')
			else:
				with open(vr_tr_transcripts, 'w') as vr_tr_writer:
					sequences = SeqIO.parse(sequences_matched_to_vr_100_area, 'fasta')
					for sequence in sequences:
						if sequence.name == vr_tr_names[0]:
							vr_tr_writer.write('>%s\n%s\n' % (sequence.name, str(sequence.seq)))
							if len(vr_tr_names) > 1:
								vr_tr_names.pop(0)
							else:
								break

				#Remove any potential TR sequences
				#Extract the TR area +/- 100 bp, then see if any of the tr_vr_sequences align perfectly to TR
				#To do: allow for 1-2 mismatches
				sequences = SeqIO.parse(formatted_ref_genomefile, 'fasta')
				for sequence in sequences:
					if sequence.name == 'Contig%i' % (TR_contig_num):
						start = TR_start - 100
						if start < 0:
							start = 0
						end = TR_end + 100
						if end > len(sequence.seq):
							end = len(sequence.seq)
						tr_area = sequence.seq[start:end]
						tr_area_rev = tr_area.reverse_complement()

				vr_sequences = '%s/VR_sequences-%s.fasta' % (output_folder, unique_name)
				num_trs = 0
				with open(vr_sequences, 'w') as vr_writer:
					sequences = SeqIO.parse(vr_tr_transcripts, 'fasta')
					for sequence in sequences:
						if sequence.seq in tr_area or sequence.seq in tr_area_rev:
							num_trs += 1
						else:
							vr_writer.write('>%s\n%s\n' % (sequence.name, str(sequence.seq)))

				#Align VR sequences and orient the VR/TR pair
				aligned_VR_sequences = '%s/aligned_VR_sequences-%s.fa' % (output_folder, unique_name)
				vr_blast_database = '%s/vrblastdb' % (temp_folder)
				vr_blast_output = '%s/vr_aligning_blastout.xml' % (temp_folder)

				temp_blast_db.append(vr_blast_database)
				temp_files.append(vr_blast_output)

				cline = NcbimakeblastdbCommandline(dbtype='nucl', input_file=VR_file, out=vr_blast_database)
				cline()

				cline = NcbiblastnCommandline(out=vr_blast_output, db=vr_blast_database, query=vr_sequences, outfmt=5, word_size=8, reward=1, penalty=-1, evalue=1e-4, gapopen=2, gapextend=1, perc_identity=50)
				cline()

				amiss, tmiss = 0, 0
				for i in range(len(VR_seq)):
					if VR_seq[i] != TR_seq[i]:
						if TR_seq[i] == 'A':
							amiss += 1
						elif TR_seq[i] == 'T':
							tmiss += 1
				reverse = False
				if tmiss > amiss:
					reverse = True

				with open(vr_blast_output, 'r') as f:
					if Path(vr_blast_output).stat().st_size > 0:
						results = NCBIXML.parse(f)
						with open(aligned_VR_sequences, 'w') as aligned_VR_writer:
							vr_oriented = VR_seq
							if reverse:
								vr_oriented = str(Seq(VR_seq).reverse_complement())
							aligned_VR_writer.write('>%s\n%s\n' % ('VR', vr_oriented))
							for result in results:
								if len(result.alignments) > 0:
									lowest = 1
									best = [0, 0]
									for anum, alignment in enumerate(result.alignments):
										for hnum, hsp in enumerate(alignment.hsps):
											if hsp.expect < lowest:
												lowest = hsp.expect
												best = [anum, hnum]
									alignment = result.alignments[best[0]]
									hsp = alignment.hsps[best[1]]
									seq_out = ''
									if hsp.sbjct_start < hsp.sbjct_end:
										start = hsp.sbjct_start - 1
									else:
										start = hsp.sbjct_end - 1
									i = 0
									for i in range(start):
										seq_out += '-'
									if hsp.sbjct_start < hsp.sbjct_end:
										seq_out += hsp.query
									else:
										seq_out += str(Seq(hsp.query).reverse_complement())
									end = len(VR_seq) - start - len(hsp.query)
									i = 0
									for i in range(end):
										seq_out += '-'
									if reverse:
										seq_out = str(Seq(seq_out).reverse_complement())
									aligned_VR_writer.write('>%s\n%s\n' % (result.query, seq_out))

	#utils.cleanup(temp_files, temp_blast_db)
	#if len(errors) > 0:
	#	utils.print_errors(error_file, errors)
	with open('%s/done.txt' % (output_folder), 'w') as w:
		w.write('Done')



if __name__ == '__main__':
	print('Please run DGR.py')
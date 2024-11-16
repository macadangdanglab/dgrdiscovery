from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline

import numpy as np
import csv, glob, math, datetime, os, time
import pandas as pd
from pathlib import Path

import DGRutils as utils

#Performs quality filtering, removing human reads, and removing viral reads
def process_raw_data(r1, r2, output_folder):
	if r2 is None:
		print('No r2 found. Proceeding with single ended analysis using r1 only')

	#Create temp_directory in output_folder to store information
	temp_folder = '%s/temp' % (output_folder)
	#os.system('rm -rf %s' % (temp_folder))
	#time.sleep(3)
	Path(temp_folder).mkdir(parents=True, exist_ok=True)

	print('Creating subfolders within output folder: 01_QC, 02_Contigs, 03_Assembly')
	QC_folder = '%s/01_QC' % (output_folder)
	Contigs_folder = '%s/02_Contigs' % (output_folder)
	Assembly_folder = '%s/03_Assembly' % (output_folder)
	Path(QC_folder).mkdir(parents=True, exist_ok=True)
	Path(Assembly_folder).mkdir(parents=True, exist_ok=True)

	#Keep track of files that need to be deleted at the end of the analysis
	temp_files = []
	
	sample_name = r1.split('/')[-1].split('.')[0].split('_')[0]

	
	#Create samples.txt for illumina-utils quality filtering
	samples_file = '%s/samples-%s.txt' % (temp_folder, sample_name)
	temp_files.append(samples_file)
	with open(samples_file, 'w') as w:
		writer = csv.writer(w, delimiter='\t')
		writer.writerow(['Sample', 'r1', 'r2'])
		writer.writerow([sample_name, r1, r2])

	#Create ini file for iu
	os.system('iu-gen-configs %s -o %s' % (samples_file, temp_folder))

	#Run quality filtering
	os.system('iu-filter-quality-minoche --ignore-deflines %s/%s.ini' % (temp_folder, sample_name))

	temp_files.append('%s/%s-QUALITY_PASSED_R1.fastq' % (temp_folder, sample_name))
	temp_files.append('%s/%s-QUALITY_PASSED_R2.fastq' % (temp_folder, sample_name))
	

	#Filter out reads that map to human genome
	'''
	if r2 is not None:
		os.system('kneaddata --input %s --input %s -db genomes/hg37dec_v0.1 --trimmomatic /home/brm4/Trimmomatic-0.39/ --output %s -t 4 --trimmomatic-options "SLIDINGWINDOW:4:20 MINLEN:50" --bowtie2-options "--very-sensitive --dovetail" --remove-intermediate-output' % (r1, r2, temp_folder))
		os.system('cp %s/%s_1_kneaddata_paired_1.fastq %s/%s_1_kneaddata_paired_1.fastq' % (temp_folder, sample_name, QC_folder, sample_name))
		os.system('cp %s/%s_1_kneaddata_paired_2.fastq %s/%s_1_kneaddata_paired_2.fastq' % (temp_folder, sample_name, QC_folder, sample_name))
		kd_files = glob.glob('%s/*kneaddata*' % (temp_folder))
		for f in kd_files:
			temp_files.append(f)
	else:
		pass
		#os.system('kneaddata --input%s --reference-db genomes/hg37dec_v0.1 --trimmomatic /home/brm4/Trimmomatic-0.39/ --output %s' % (r1, QC_folder))
	'''
	#Filter out reads that map to human genome using bowtie2
	if r2 is not None:
		#os.system('bowtie2 --very-sensitive --dovetail --quiet --threads 4 -x genomes/hg37dec_v0.1 -1 %s/%s-QUALITY_PASSED_R1.fastq -2 %s/%s-QUALITY_PASSED_R2.fastq --un-conc %s/%s_clean.fastq --al-conc %s/%s_contam.fastq' % (temp_folder, sample_name, temp_folder, sample_name, QC_folder, sample_name, temp_folder, sample_name))
		aligned_file = '%s/%s-aligned.sam' % (temp_folder, sample_name)
		temp_files.append(aligned_file)
		unmapped_bam = '%s/%s-unmapped.bam' % (temp_folder, sample_name)
		temp_files.append(unmapped_bam)
		sorted_unmapped_bam = '%s/%s-sorted_unmapped.bam' % (temp_folder, sample_name)
		temp_files.append(sorted_unmapped_bam)
		unmapped_fq1 = '%s/%s-unmapped1.fastq' % (QC_folder, sample_name)
		unmapped_fq2 = '%s/%s-unmapped2.fastq' % (QC_folder, sample_name)
		contam_file = '%s/%s_contam.fastq' % (temp_folder, sample_name)
		temp_files.append(contam_file)
		print('Aligning with bowtie2...')
		os.system('bowtie2 --very-sensitive --dovetail --threads 4 -x genomes/hg37dec_v0.1 -1 %s/%s-QUALITY_PASSED_R1.fastq -2 %s/%s-QUALITY_PASSED_R2.fastq -S %s' % (temp_folder, sample_name, temp_folder, sample_name, aligned_file))
		print('Converting sam to bam')
		os.system('samtools view --threads 4 -bS -f 12 %s > %s' % (aligned_file, unmapped_bam))
		print('Sorting bam')
		os.system('samtools sort -m 4G --threads 4 -n -o %s %s' % (sorted_unmapped_bam, unmapped_bam))
		print('Converting bam to fastq')
		os.system('samtools fastq --threads 4 -1 %s -2 %s -0 /dev/null -s /dev/null -n %s' % (unmapped_fq1, unmapped_fq2, sorted_unmapped_bam))
	else:
		pass
	

	unmapped_fq1 = '%s/%s-unmapped1.fastq' % (QC_folder, sample_name)
	unmapped_fq2 = '%s/%s-unmapped2.fastq' % (QC_folder, sample_name)
	
	#Use megaHIT to asssemble contigs
	if r2 is not None:
		print('Assembling contigs with megahit')
		os.system('megahit -1 %s -2 %s --min-contig-len 500 -m 0.85 --presets meta-sensitive -o %s -t 4' % (unmapped_fq1, unmapped_fq2, Contigs_folder))
	else:
		pass
	
	#Remove intermediate contigs
	files = glob.glob('%s/intermediate_contigs/*' %(Contigs_folder))
	for f in files:
		temp_files.append(f)

	#Rename_contigs:
	os.system('mv %s/final.contigs.fa %s/%s_contigs.fa' % (Contigs_folder, Contigs_folder, sample_name))
	
	#Map reads back to contigs using bowtie2
	#Create bowtie2 database from the contigs file
	os.system('bowtie2-build %s/%s_contigs.fa %s/%s_db' % (Contigs_folder, sample_name, Contigs_folder, sample_name))
	bowtie2_files = glob.glob('%s/*.bt2' % (Contigs_folder))
	for f in bowtie2_files:
		temp_files.append(f)

	#Align to contigs
	print('Running bowtie2 to align reads to contigs')
	aligned_to_contigs_sam = '%s/%s-aligned.sam' % (Assembly_folder, sample_name)
	temp_files.append(aligned_to_contigs_sam)
	os.system('bowtie2 --sensitive-local --threads 4 -x %s/%s_db -1 %s -2 %s -S %s' % (Contigs_folder, sample_name, unmapped_fq1, unmapped_fq2, aligned_to_contigs_sam))

	print('Converting sam to bam')
	aligned_to_contigs_bam = '%s/%s-aligned.bam' % (Assembly_folder, sample_name)
	temp_files.append(aligned_to_contigs_bam)
	os.system('samtools view -F 4 --threads 4 -bS -o %s %s' % (aligned_to_contigs_bam, aligned_to_contigs_sam))

	print('Sorting bam file')
	aligned_to_contigs_bam_sorted = '%s/%s-aligned_sorted.bam' % (Assembly_folder, sample_name)
	#temp_files.append(aligned_to_contigs_bam_sorted)
	os.system('samtools sort -m 7G --threads 2 -o %s %s' % (aligned_to_contigs_bam_sorted, aligned_to_contigs_bam))

	utils.cleanup_files(temp_files)
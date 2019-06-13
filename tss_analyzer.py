#!/usr/bin/env python3

import sys
import argparse
import pysam
import HTSeq
import Readers
import os
from operator import itemgetter
from scipy.stats import poisson 

"""
Analyzes a BED file with scores that indicate read density at TSS positions (as provided by dbtss).
The script than analyzes regions surrounding annotated TSSs of known genes using a peak finder
algorithm based loosely on sliding window possion model (Strnenac D et al, BMC Genomics, 2013).
Calculates per gene sum(nr of reads at each position within a peak * nr of subsequent T/C or A/G at that position) / total nr of reads within all peaks
"""


def main() :
	
	parser = argparse.ArgumentParser( description='Calclulates the distribution of T to C mismatches in a bam file for a given gene' )
	parser.add_argument( 'gtf_file' )
	parser.add_argument( 'genome' )
	parser.add_argument( 'bed_file', help="bed or bam file, stdin (bed-file) if 'STDIN' is indicated instead of a file name (in that case an output file name needs to be specified")
	parser.add_argument( '--gene', help="Reports TSS for 1 gene only on STDOUT.")
	parser.add_argument( '--mode', default = 'top', help="top: calculate top-score (default), TSS: define TSS position + sequence, nt_freq: calculate nt frequencies in first 10 nt multiplied by transcription frequency from site, nt_freq_position: calculate nt frequencies per position in first 10 nt multiplied by transcription frequency from site" )
	parser.add_argument( '--outfile', help="output file name" )
	options = parser.parse_args()

	gtf_file = options.gtf_file
	bed_file = options.bed_file
	genome_file = options.genome
	gene = options.gene
	mode = options.mode
	outfile = options.outfile
	
	
	# Load BED file into GA
	# Type is: BED, BAM, or STDIN
	filetype = 'BED'
	if bed_file[-3:] == 'bam' : filetype = 'BAM'
	elif bed_file == 'STDIN' : filetype = 'STDIN'
	ga = initialize_ga( bed_file, filetype )
	
	# Load GTF file
	gtf = Readers.GTF_Reader( gtf_file )

	# Load genome fasta
	genome = pysam.FastaFile( genome_file )

	# determine the output file basename
	if outfile : outfile_basename = outfile
	else : outfile_basename = bed_file

	# determine which analysis to do.
	if gene : 
		process_gene( gene, gtf, ga, genome )
	
	elif mode == 'TSS':
		output_file = os.path.join(os.getcwd(), os.path.basename(outfile_basename)) + ".tss"
		process_tss_max( gtf, ga, genome, output_file )
		
	else :
		output_file = os.path.join(os.getcwd(), os.path.basename(outfile_basename)) + ".topscore"
		
		if mode == 'nt_freq':
			output_file = output_file + '_nt_freq'
		if mode == 'nt_freq_position':
			output_file = output_file + '_nt_freq_position'
		process_tss_multi( gtf, ga, genome, output_file, mode )

def initialize_ga( filename, filetype ) :
	"""
	Initialize an HTseq genomic array with TSS data from BED or BAM files.
	"""
	if filetype == "BED" : ga = load_bedfile_to_ga( filename )
	elif filetype == "BAM" : ga = load_bamfile_to_ga( filename )
	elif filetype == "STDIN" : ga = load_stdin_bed_to_ga()
	else :
		print("Can't load TSS data!")
		exit()
	
	return ga


def process_gene( gene, gtf, ga, genome ) :
	"""
	Report TSS for a single gene to STDOUT.
	"""

	for g in gtf.features() :
		if g.gene_id() != gene : continue
		
		print("%s (%s)" % (gene, g.transcript_id()))
		qiv = get_tss_query_interval( g, upstream=1000, downstream=1000 )

		try :
			qiv_scores = ga[ qiv ]
		except IndexError :
			print( qiv )
	
		# get the read densities from the BED file for the requested region
		qiv_coverage_vector = chromvector_to_array( qiv_scores )

		# search for peaks
		peak_regions = find_peak_regions( qiv_coverage_vector, window=50, adjacent=200 )

		# retrieve tss sequences within all peaks
		for piv in peak_regions : 
			
			coverage = qiv_coverage_vector[ piv[0]:piv[1] ]
			genomic_iv = (qiv.chrom, qiv.start + piv[0], qiv.start + piv[1], g.strand())
			sequences = get_tss_sequences( genomic_iv, coverage, genome )

			print("%s:%d-%d,%s" % genomic_iv)
			for s in sequences : 
				print("%s\t%d" % (s[1],s[0]))


def process_tss_multi( gtf, ga, genome, output_file, mode ) :
	"""
	Calcuate a TOPscore and NT frequencies that use all reads in a TSS peak.
	"""

	scored = []	
	ntcounter = NtFreqCounter()
	
	for gene in gtf.features() :
		
		# returns per gene a list of (score, sequence until 15 nt downstream) for all positions within all peaks with score > 0
		sequences = retrieve_gene_tss( gene, ga, genome )

		if mode == 'top':
			if len(sequences) > 0 :
				#max_tss = max( sequences, key= lambda x: x[0] )
				#ntcounter.record( max_tss[1] )
				countsum = sum( [ s[0] for s in sequences ] )
				sum_cutoff = 0
				base_cutoff = 0
				if countsum > sum_cutoff :
					for s in sequences : 
						if s[0] > base_cutoff :
							ntcounter.record( s[1], count=(s[0]/float(countsum)) )
				
				scores = score_tss_sequences( sequences )
				scored.append((
					gene.transcript_id(), 
					gene.gene_id(), 
					scores['numreads'], 
					scores['avg_ct_len'], 
					scores['avg_ag_len'] ))

		if mode == 'nt_freq':
			scores = score_tss_sequences_nt( sequences )
			scored.append((
				gene.transcript_id(), 
				gene.gene_id(), 
				scores['numreads'], 
				scores['C'],
				scores['T'],
				scores['G'], 
				scores['A'], 
				scores['CT'], 
				scores['CG'], 
				scores['CA'], 
				scores['GT'], 
				scores['GA'], 
				scores['TA'] ))

		if mode == 'nt_freq_position':
			scores = score_per_nt_position( sequences )
			scored.append((
				gene.transcript_id(), 
				gene.gene_id(), 
				scores['numreads'], 
				','.join([str(n) for n in scores['C']]), 
				','.join([str(n) for n in scores['T']]),
				','.join([str(n) for n in scores['G']]),
				','.join([str(n) for n in scores['A']])
				))
	
	header_start = 	('transcript','gene','numreads')
	if mode == 'top':		
		header = "\t".join( header_start + ('ctlen','aglen') )
		
	if mode == 'nt_freq':
		header = "\t".join( header_start + ('C','T','G','A','CT','CG','CA','GT','GA','TA') )
		
	if mode == 'nt_freq_position':
		header = "\t".join( header_start + ('C','T','G','A') )
		
	write_file( output_file, header, scored )
	ntcounter.write( output_file + ".freqs" ) 
		

def process_tss_max( gtf, ga, genome, output_file ):
	"""
	Identical to process_tss_multi, but only considers the single position within the TSS region that has the maximum number of reads.
	Output: file with per transcript position of position with highest score in all TSS peak regions along with first 15 nt.
	"""
	
	of = open(output_file, 'w')
	header = ('transcript','gene','chrom','position','strand','numreads','numreads_peakregion','peak_width','sequence')
	of.write("\t".join(header) + "\n")

	ntcounter = NtFreqCounter()

	for gene in gtf.features() :
		
		# define region start of gene +/- 1000 nt
		qiv = get_tss_query_interval( gene, upstream=1000, downstream=1000 )

		if qiv.chrom == 'chrM' : 
			continue
	
		# get the read densities from the BED file for the requested region
		try :
			qiv_scores = ga[ qiv ]
		except IndexError :
			print( qiv )
			continue

		#convert scores in interval to array (list)
		qiv_coverage_vector = chromvector_to_array( qiv_scores )
	

		# search for peaks (window of 50 nt is compared to adjacent 200 nt windows, 
		# if significantly enriched -> peak, is then trimmed of zeros and merged with overlapping peaks, yields list of peaks(start,end)
		peak_regions = find_peak_regions( qiv_coverage_vector, window=50, adjacent=200 )
		
		
		# returns list of (position, score, numreads in peak) for all positions in all peaks that have this max score (calculated over all peaks). 
		tss = find_TSS( peak_regions, qiv_coverage_vector )
		
		for t in tss:
			try:
				first_nts = genome.fetch( qiv.chrom, qiv.start + t[0], qiv.start + t[0]+15 ).upper() 
				if qiv.strand == '-':
					first_nts = genome.fetch( qiv.chrom, qiv.start + t[0]-14, qiv.start + t[0]+1 ).upper()
					first_nts = reverse_complement(first_nts)
				
				row = (gene.transcript_id(), gene.gene_id(), qiv.chrom, str(qiv.start + t[0]+1), 
					qiv.strand, str(t[1]), str(t[2]), str(t[3]), first_nts)
				of.write("\t".join(row) + "\n") 

				ntcounter.record( first_nts, count=t[2] )

	
			except KeyError :
				print('%s: sequence not available' %(gene.transcript_id()))		
	
	of.close()		
	ntcounter.write( output_file + ".freq" )


def retrieve_gene_tss( gene, ga, genome ) :
	
	"""
	Find peaks within a TSS and return a list of counts, sequences.
	"""
	
	# define region start of gene +/- 1000 nt
	qiv = get_tss_query_interval( gene, upstream=1000, downstream=1000 )

	if qiv.chrom == 'chrM' : 
		return []
	
	# get the read densities from the BED file for the requested region
	try :
		qiv_scores = ga[ qiv ]
	except IndexError :
		print( qiv )
		return []

	#convert scores in interval to array (list)
	qiv_coverage_vector = chromvector_to_array( qiv_scores )
	

	# search for peaks (window of 50 nt is compared to adjacent 200 nt windows, 
	# if significantly enriched -> peak, is then trimmed of zeros and merged with overlapping peaks, yields list of peaks(start,end)
	peak_regions = find_peak_regions( qiv_coverage_vector, window=50, adjacent=200 )

	# retrieve tss sequences within all peaks
	sequences = []
	for piv in peak_regions : 
			
		coverage = qiv_coverage_vector[ piv[0]:piv[1] ]
		genomic_iv = (qiv.chrom, qiv.start + piv[0], qiv.start + piv[1], gene.strand())
		sequences.extend( get_tss_sequences( genomic_iv, coverage, genome ))
		#returns a list of (score, sequence from position to 15 nt downstream) for all positions with score > 0 within a peak
		
	return sequences


def write_file( output_file, header, data ) :
	
	with open(output_file, "w") as ofh : 
		
		ofh.write( header + "\n")

		for item in sorted( data, key=itemgetter(3) , reverse=True) : 
			row = [ str(r) for r in item ] 
			ofh.write( "\t".join( row ) + "\n" )
	

def reverse_complement( seq ) :
	rc = {'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C', 'N' : 'N'}
	rcseq = "".join( [rc[nt] for nt in seq[::-1]])
	return rcseq

def load_bedfile_to_ga( bed_file ) :
	
	ga = HTSeq.GenomicArrayOfSets( "auto", stranded=True )
	with open(bed_file, "r") as fh :
		
		for line in fh : 

			row = line.strip().split("\t")

			# field values based on the bed files from DBTSS. BED files are 0 based, no adjusting necessary.
			try :	
				chrom = row[0]
				pos = int(row[1])
				strand = row[5]
				score = float(row[4])
			except IndexError as e :
				print(e)
				print(row)
				continue
			
			try :
				ga[ HTSeq.GenomicInterval( chrom, pos-1, pos, strand )] = score
			except ValueError as e:
				print("Error loading GA:")
				print(row)
				continue
	return ga

def load_bamfile_to_ga( bam_file ) :
	
	ga = HTSeq.GenomicArray( "auto", stranded=True )
	bam_reader = HTSeq.BAM_Reader( bam_file )

	plusone_mismatch_reads = 0
	plusone_mismatch = {'A':0, 'T':0, 'C':0, 'G':0}
	total_reads = 0
	
	for aln in bam_reader:

		if aln.aligned:

			total_reads += 1
			ivstart = aln.iv.start_d

			#if has_plusone_mismatch( aln ) :
		#		#print(aln.read.seq, aln.iv.strand, aln.optional_field("MD"))
		#		plusone_mismatch_reads += 1
		#		plusone_mismatch[ aln.read.seq.decode('UTF-8')[0] ] += 1
		#		if aln.iv.strand == '+' : ivstart += 1
		#		else : ivstart -= 1

			iv = HTSeq.GenomicInterval( aln.iv.chrom, ivstart, ivstart+1, aln.iv.strand )
			ga[iv] += 1

		# FOR REMOVING MULTIMAPPERS
		#if "NH" in aln.optional_fields:
		#		if aln.optional_field("NH") == 1:
	#				ga[iv] += 1
#			else:
#				if "XS" not in aln.optional_fields:
#					ga[iv] += 1
	
	print("+1 mismatches: %d/%d" % (plusone_mismatch_reads, total_reads) )
	print( plusone_mismatch) 
	return ga

def has_plusone_mismatch( aln ) : 
	"""
	Check whether an alignment has a mismatch in the +1 position, which generally indicates an artifact of RT
	"""
	mdstring = aln.optional_field('MD')
	readseq = aln.read.seq.decode('UTF-8')
	strand = aln.iv.strand
	offset = 0
	
	# test whether read has a mismatched +1 G
	# note that htseq reports the origial read sequence for both plus and minus strand reads
	# but that the MD string is always the plus strand alignment.
	if (( strand=='+' and mdstring[0]=='0')  
		or ( strand=='-' and mdstring[-1]=='0' and not(mdstring[-2].isdigit()) )) :
		#print( readseq, strand, mdstring )	
		return True
	else : 
		return False

def load_stdin_bed_to_ga() :
	
	ga = HTSeq.GenomicArrayOfSets( "auto", stranded=True )
		
	for line in sys.stdin : 

		row = line.strip().split("\t")

		# field values based on the bed files from FANTOM. BED files are 0 based, no adjusting necessary.
		try :	
			chrom = row[0]
			pos = int(row[1])
			strand = row[5]
			score = float(row[4])
		except IndexError as e :
			print(e)
			print(row)
			continue

		ga[ HTSeq.GenomicInterval( chrom, pos, pos+1, strand )] = score
		
	return ga

def get_tss_query_interval( gene, upstream=0, downstream=0 ) :
	"""
	Return an HTSeq.GenomicInverval object for the region to search for the TSS.
	Parameters on that remain to be defined.
	"""

	exons = gene.exons()
	strand = gene.strand()

	if strand == '+' : 
		start = exons[0].start - upstream 
		end = exons[0].start + downstream 
	else :
		end = exons[-1].end + downstream 
		start = exons[-1].end - upstream 
	
	return (HTSeq.GenomicInterval( exons[0].chr, start, end, strand))

def chromvector_to_array( chromvector ) :
	"""
	Transform an HTSeq chromvector object to a numeric array.
	"""

	v = []
	for s in chromvector.steps() :
		ivlen = s[0].end - s[0].start
		if isinstance( s[1], float ) : v.extend( [s[1]]*ivlen )
		else : v.extend( [0]*ivlen )	

	return v

def find_peak_regions( iv, window, adjacent ) :

	"""
	Searches a vector of scores for peaks according to a simple poisson enrichment score.
	In short, the expected score within the main window (test) is compared to the score density 
	in the adjoining windows (b1 = 200 nt upstream, b2 = 200 nt downstream). 
	Significance of enrichment is calcuated via Poisson.
	If region is significantly enriched, it's trimmed for bounding zeros and returned.
	"""

	winsize = window + 2*adjacent
	peaks = []

	for winstart in range( 0, len(iv) - winsize, int(window/2) ) :
		
		b1 = sum( iv[ winstart:winstart+adjacent ] )
		b2 = sum( iv[ winstart+window+adjacent : winstart+window+adjacent*2 ] )
		test = sum( iv[ winstart+adjacent : winstart+adjacent+window ] )

		if test == 0 : continue
		
		# poisson test raises nan error if lambda == 0
		if b1 == 0 : b1 = 1
		if b2 == 0 : b2 = 1

		p1 = poisson.sf( test, b1/4.0 )
		p2 = poisson.sf( test, b2/4.0 )

		cutoff = 0.000000001
		if (p1 < cutoff) & (p2 < cutoff) :
			
			start = winstart + adjacent
			end = winstart + adjacent + window
			while iv[start] == 0 : 
				start += 1
			while iv[end] == 0 : 
				end -= 1
			peaks.append( (start, end+1) )

	# merge overlapping peaks
	merged = []
	if len(peaks) > 0 :
		prev = peaks[0]
		for i in range(1, len(peaks)) : 
			if prev[1] > peaks[i][0] : prev = (prev[0], peaks[i][1])
			else : 
				merged.append( prev )
				prev = peaks[i]
		merged.append( prev )

	return merged


def find_TSS( peaks, iv ):
	
	"""
	Identifies position(s) with maximum value over all peaks 
	Returns: list of (position, position_counts, interval_counts, interval_width) for all positions with counts=max
	"""
	
	tss = []
	
	for p in peaks:
		counts = sum( iv[p[0]:p[1]] )
		max_peak_score = max( iv[p[0]:p[1]] )
		for i in range(p[0], p[1]):
			if iv[i] == max_peak_score:
				tss.append( (i, iv[i], counts, p[1]-p[0]) )
	
	if len(tss) > 0:
		max_score = max ( [tss[i][1] for i in range(0,len(tss))] )
		tss_max = list(set([t for t in tss if t[1] == max_score]))
	else:
		tss_max = []
		
	return tss_max

def get_tss_sequences( iv, coverage, genome, seqlen=15 ) :
	
	"""
	Args: genomic interval, read coverage vector across interval, pysam.Faidx object for genome.
	returns: array of ( count, seq ) for all sequences of length seqlen.
	"""

	try : 
		seq = genome.fetch( iv[0], iv[1]-seqlen, iv[2]+seqlen ).upper()
	except KeyError :
		return []
	
	if iv[3] == '-' : 
		seq = reverse_complement(seq)
		coverage = coverage[::-1]

	sequences = [ (coverage[i], seq[i+seqlen:i+seqlen*2]) for i in range(0,len(coverage)) if coverage[i] > 0]
	
	return sequences
	
def score_tss_sequences( sequences ) :
	"""
	Calculate the average unbroken C/T and A/G lengths for all TSS sequences included in the list of sequences. The weighted average
	of C/T lengths is TOPscore. Other calculations of TSS features can also be included here in the future.
	Args: list of ( count, sequence ) pairs for a TSS peak
	Returns: dict of score_name -> score
	"""

	scores = dict() 
	scores['numreads'] = sum( [s[0] for s in sequences] )

	if scores['numreads'] > 0 : 
		scores['avg_ct_len'] = sum( [ s[0]*calculate_ct_len(s[1]) for s in sequences ] ) / float(scores['numreads'])
		scores['avg_ag_len'] = sum( [ s[0]*calculate_ag_len(s[1]) for s in sequences ] ) / float(scores['numreads'])
	else : 
		scores['avg_ct_len'] = 0.0
		scores['avg_ag_len'] = 0.0

	return scores
	
def score_tss_sequences_nt( sequences ):
	
	scores = {'C':0,'G':0,'T':0,'A':0,'CT':0,'CG':0,'CA':0,'GT':0,'GA':0,'TA':0}
	scores['numreads'] = sum( [s[0] for s in sequences] )
	
	if scores['numreads'] > 0:
		for s in sequences:
			for k in calculate_nt_freq( s[1] ).keys():
				scores[k] += s[0]*calculate_nt_freq( s[1] )[k] / float(scores['numreads'])
	
	return scores 

def score_per_nt_position( sequences ):
	"""
	Calculate the frequency of nucleotides at each position in a list of numreads,sequence pairs.
	Args: List of (numreads, seq)
	Returns: hash of A,T,C,G -> [ freq per position]
	"""

	scores = {'C':10*[0],'G':10*[0],'T':10*[0],'A':10*[0]}
	scores['numreads'] = sum( [s[0] for s in sequences] )	
	
	if scores['numreads'] > 0:
		for s in sequences:
			for i in range(10):
				scores[s[1][i]][i] += s[0]/scores['numreads']
	
	return scores


def calculate_ct_len( seq ) :

	i=0
	while i < len(seq) :
		if (seq[i] == 'C') or (seq[i] == 'T') : i += 1
		else : break
	
	return i
	
def calculate_ag_len( seq ) :
	i=0
	while i < len(seq) :
		if (seq[i] == 'A') or (seq[i] == 'G') : i += 1
		else : break
	return i

def calculate_nt_freq( seq ) :
	
	freq = {'C':0,'G':0,'T':0,'A':0,'CT':0,'CG':0,'CA':0,'GT':0,'GA':0,'TA':0}
	for i in range(0,10):
		for k in freq.keys():
			if seq[i] in k:
				freq[k] += 1
	
	return freq
	

class NtFreqCounter :
	
	"""
	Simple class for recording nt frequencies at each position in a set of sequences.
	"""

	def __init__( self ) :
		self.ntlen = 10
		self.count = 0
		self.freqs = [None]*self.ntlen
		for i in range(0, self.ntlen) : self.freqs[i] = {'A':0.1, 'C': 0.1, 'T': 0.1, 'G': 0.1, 'N': 0.1 }

	def record( self, seq, count=1 ) :
		
		for i in range(0, len( self.freqs )) :
			if i < len(seq) :
				self.freqs[i][ seq[i] ] += count

		self.count += count

	def write( self, filename ) :
		
		with open( filename, "w" ) as ofh :

			# print header
			nts = ('A', 'T', 'C', 'G', 'N')
			ofh.write( "Pos\t" + "\t".join( nts ) + "\n" )

			for i in range(0, len(self.freqs) ) :
				ntcounts = [ self.freqs[i][nt] for nt in nts ]
				ntfreqs = [ str(nt/sum(ntcounts)) for nt in ntcounts ]
				ofh.write( str(i) + "\t" + "\t".join(ntfreqs) + "\n" )
	

if __name__ == "__main__" : 
	main()

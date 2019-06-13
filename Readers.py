import HTSeq
import copy
import sys
import re
from collections import namedtuple
from operator import attrgetter
import os

Interval = namedtuple( 'Interval', 'chr start end strand frame type attrs' )

class Feature :

	"""
	Describes a genomic feature consisting of a list of genomic intervals and methods for accessing feature attributes.

	Parameters:
	intervals : a list of Interval named tuples, most often generated from a GFF file.

	"""

	# instantiate with a list of Interval objects
	def __init__(self, intervals) : 
		self._intervals = intervals

	# access basic attributes
	def get_attrs(self) : 
		return self._intervals[0].attrs	

	def transcript_id(self) : 
		return self.get_attrs()['transcript_id']

	def gene_id(self) : 
		attrs = self.get_attrs()
		if 'gene_id' in attrs : return attrs['gene_id']
		else : return None

	def strand(self) : 
		return self._intervals[0].strand
		
	# retrieve intervals of a specific type (eg. exon, CDS, etc.). Default is none.
	def get(self, feature_type=None ) :
		e = [f for f in self._intervals if f.type == feature_type]
		return e

	def exons(self) : 
		return self.get('exon')
	
	def introns(self) :

		exons = self.exons()
		introns = []

		for i in range( 1, len(exons) ) :

			iv = Interval( chr = exons[ i ].chr,
							type = 'intron',
							start = exons[ i-1 ].end,
							end = exons[ i ].start,
							strand = exons[ i ].strand,
							frame = ".",
							attrs = exons[ i ].attrs
							)

			if ( iv.end - iv.start ) > 0 : 
				introns.append( iv )

		return introns

	def cds(self) : 
		return self.get('CDS')

	def five_utr(self) :
		strand = self.strand()
		if strand == '+' : return self.get_start_utr()
		else : return self.get_end_utr()
			
	def three_utr(self) : 
		strand = self.strand()
		if strand == '+' : return self.get_end_utr()
		else : return self.get_start_utr()
	
	# returns the utr at the left end of the transcript
	def get_start_utr(self) :
		exons = sorted(self.exons(), key=attrgetter('start'))

		cds_exons = self.cds()
		if len( cds_exons ) == 0 : return []

		cds_start = min([x.start for x in cds_exons])
			
		utr = []
			
		for exon in exons :
			iv_start, iv_end = (0,0)
			
			if exon.end < cds_start : 
				iv_start = exon.start
				iv_end = exon.end
			else : 
				iv_start = exon.start
				iv_end = cds_start
			
			if ( iv_end - iv_start ) > 0 : 
				iv = Interval( chr = exon.chr,
							   type = 'utr',
							   start = iv_start,
							   end = iv_end,
							   strand = exon.strand,
							   frame = exon.frame,
							   attrs = exon.attrs
							   )
				utr.append( iv )
			
			if exon.end >= cds_start : break
			
		return utr
		
	# returns the utr at the right end of the transript
	def get_end_utr(self) :

		# exons are returned sorted by increasing position
		exons = sorted(self.exons(), key=attrgetter('start'), reverse=True)

		cds_exons = self.cds()
		if len(cds_exons) == 0 : return []

		cds_end = max([x.end for x in cds_exons])
		utr = []
			
		for exon in exons :
			iv_start, iv_end = (0,0)
			
			if exon.start > cds_end : 
				iv_start = exon.start
				iv_end =- exon.end
			else : 
				iv_start = cds_end
				iv_end = exon.end

			if ( iv_end - iv_start ) > 0 : 

				iv = Interval( chr = exon.chr,
							   type = 'utr',
							   start = iv_start,
							   end = iv_end,
							   strand = exon.strand,
							   frame = exon.frame,
							   attrs = exon.attrs)

				utr.append( iv )			  
			
			if exon[1] <= cds_end : break
			
		return utr

			
	# methods for accessing transcript coordinates
	# arguments are closed end genomic coordinates
	# returns a list of half-open intervals
	def genomic_fragment(self, start, end) :
		exons = self.exons()
		if start > end : 
			s = end 
			e = start
		else : 
			s = start
			e = end 

		# get chr and strand information
		strand = self.strand()
		chrom = self.chrom()
			
		intervals = []
		i = 0
		curr_f = 0
		exon = exons[i]

		# check bounds, and truncate start, end if necessary
		if s < exons[0][1] : s = exons[0][1]
		if e >= exons[-1][2] : e = exons[-1][2] - 1

		# search for the first exon
		while True :
			
			# if the interval start is right of the exon bounds, go to the next exon
			# else start recording intrevals
			if (exon[1] <= s) and (exon[2] > s) : 
					curr_f = s
					break
			
			i += 1
			if i < len(exons) :
				exon = exons[i]
			else :
				return None

		# search for the last exon, appending exons along the way
		while True :
			if exon[2] <= e :
				intervals.append((chrom, curr_f, exon[2], strand))
			else :
				if exon[1] <= e :
					intervals.append((chrom, curr_f, e+1, strand))
					break
			i += 1
			if i < len(exons) :
				exon = exons[i]
				curr_f = exon[1]
			else :
				break
			
		return intervals
					
		
	def tx_fragment(self, start, end) :
		# if start and end are outside of bounds (start < 1, end > tx_length)
		# extend 1st and last exon, respectively
		
		genomic_start = self.tx_to_genomic(start)
		genomic_end = self.tx_to_genomic(end)

		return self.genomic_fragment(genomic_start, genomic_end)

	def tx_to_genomic(self, tx_pos) :

		# force the position to be in bounds
		tx_length = self.tx_length()
		if tx_pos < 1 : tx_pos = 1
		if tx_pos > tx_length : tx_pos = tx_length
		
		# make position relative to left end
		if self.strand() == '-' :
			tx_pos = tx_length - tx_pos + 1
		
		# find genomic position
		exons = self.exons()
		if tx_pos < 1 : 
			return exons[0][1] + tx_pos -1
		elif tx_pos > tx_length : 
			return exons[-1][2] + (tx_length - tx_pos) - 1
		else :
			running_total = 0
			for e in self.exons() :
				exon_length = e[2]-e[1]
				running_total += exon_length
				diff = tx_pos - running_total
				if diff <= 0 :
					return e[2]+diff-1

		return None
				

	# note that this returns 1-based coordinates
	def genomic_to_tx(self, genomic_pos) :

		exons = self.exons()
		position = 0
		
		# if the position is out of bounds, return none.
		if genomic_pos < exons[0][1] : return None
		if genomic_pos >= exons[-1][2] : return None

		for e in exons :
			if genomic_pos > e[2] :
				position += e[2]-e[1]
			else :
				position += genomic_pos-e[1]+1
				break

		# calcuate cds length
		tx_len = self.tx_length()

		if self.strand() == '+' : return position
		else : return tx_len-position+1


	def tx_cds_coordinates(self) :
		cds = self.cds()
		start = min([x[1] for x in cds])
		end = max([x[2] for x in cds])
		
		tstart = self.genomic_to_tx(start)
		tend = self.genomic_to_tx(end-1)
		
		if tend > tstart : return (tstart, tend)
		else : return (tend, tstart)
			
	# Methods for generating coverage vectors

	def coverage(self, samfile) : 
		return self.coverage_over_intervals(self.exons(), samfile)

	def coverage_over_intervals(self, intervals, samfile) : 
		int_length = sum([s[2]-s[1] for s in intervals])

		coverage = [0]*int_length
		position_offset = 0
	
		for interval in intervals : 
			alignedreads = samfile.fetch(reference=interval[0], start=interval[1]-1, end=(interval[2]))

			for read in alignedreads : 
			
				position = read.pos
			
			#aend is in 1-based coordinates -- this appears to be a pysam error
			#or, it's half-open, but that's not documented.
				if read.is_reverse : position = read.aend - 1

				if position >= interval[1] and position < interval[2] :
					position = position - interval[1] + position_offset
					coverage[position] += 1

			position_offset += interval[2] - interval[1]
		
		if intervals[0][3] == '+' : return coverage
		else : return list(reversed(coverage))

# READERS

############################################################################
# GTF_Reader: Class for parsing GTF files. 
# NOTE: returns intervals in 0-based coordinates per Python convention, NOT the 1-based coordinates in original GTF files.
############################################################################
class GTF_Reader : 
	
	def __init__(self, filename, id_tag = 'transcript_id') : 
		self._filename = filename
		self._file = open( filename, "r" )

		#process first row of gff
		line = self._file.readline().strip()
		self._current_iv = self.parse_row( line ) 

		#set which identifier from the attributes field to use as the id
		self._idtag = id_tag

	def records(self) : 
		record = self.get_next_record()
		while( record ) : 
			yield record
			record = self.get_next_record()

	def features(self) : 
		record = self.get_next_record()
		while( record ) :
			yield Feature( record )
			record = self.get_next_record()

	# return next record as a list of intervals
	# return None if no more records exist
	def get_next_record(self) :
		
		current_iv = self._current_iv
		current_gene = None 
		record = None

		# if current_iv is None, there are no more records to return.
		if current_iv : 
			record = [ current_iv ]
			current_gene = current_iv.attrs[ self._idtag ]
		else : 
			return None

		while( True ) : 
	
			next_line = self._file.readline()

			if next_line :	
				next_iv = self.parse_row( next_line.strip() )
				next_gene = next_iv.attrs[ self._idtag ]
				if next_gene == current_gene :
					record.append( next_iv )
				else : 
					self._current_iv = next_iv
					break
		
			else :
				self._current_iv = None
				break

		return record
					

	def parse_row(self, line) : 
			
		# limit 
		valid_chromosome = re.compile("chr[\dMXY]+$") 
		re_tag = re.compile("\s*(\w+) \"(.+?)\"")
				
		row = line.strip().split("\t")
		# parse tags from field 8 (includes transcript_id, etc.
		
		tags = dict()

		for t in row[8].split(";") : 
			match = re.match(re_tag, t)
			if match : 
				tags[ match.group(1) ] = match.group(2)
				
		#create an interval record
		#also, GTF files are 1-based. Make 0-based.
		iv = Interval( chr = row[0],
			start = int( row[3] ) - 1,
			end = int( row[4] ),
			strand = row[6],
			frame = row[7],
			type = row[2],
			attrs = tags)	
		return iv



class Super_GFF_Reader( HTSeq.GFF_Reader ) :

	""" 
	Class that inherits from and is a wrapper for HTseq.GFF_Reader. Returns
	HTseq.GenomicFeature objects grouped by transcript_id (or gene_id) instead
	of line by line. Additionally allows for retrieval of only the longest
	transcript for a given gene_id or a union of all intervals reported for a
	given gene id.  
	"""
	
	def __init__( self, filename, end_included=True ) :
		
		super().__init__( filename, end_included )


	def fetch_by_transcript( self, ivtype=None, id='transcript_id' ) :

		"""
		Returns intervals grouped by transcript matching a specified type (eg. exon, CDS, etc.)
		Parameters:
		ivtype: type of interval to return. Includes CDS, 5UTR, 3UTR, exons, introns.
		id: 'gene_id' if used with union2gene gtf file
		"""

		interval_type_methods = { 
			'exon' : self.get_exons,
			'intron' : self.get_introns,
			'5UTR' : self.get_5utr,
			'3UTR' : self.get_3utr,
			'CDS' : self.get_cds,
			'union2gene' : self.get_union2gene
		}

		current_record= [] 
		current_id = None

		for row in self :

			row_id = row.attr[id]

			if current_id == row_id :
				current_record.append( row )
			else :

				if current_id :

					if ivtype : 
						current_record = interval_type_methods[ ivtype ]( current_record )

					if current_record :
						yield current_record

				current_id = row_id
				current_record = [ row ]

		if ivtype :
			current_record = interval_type_methods[ ivtype ]( current_record )
		if current_record :
			yield current_record
			

	def fetch_by_gene( self, ivtype=None, method='union' ) :
		
		"""
		Return intervals from a GTF file grouped by gene.
		
		Parameters:
		
		ivtype: class of intervals to return. Includes CDS, 5UTR, 3UTR, exons, introns.
		method: method for collapsing transcripts into single entries for each gene.
		"""

		merge_method = {
			'union' : self.get_union, 
			'longest' : self.get_longest 
		}

		current_records = []
		current_id = None

		for record in self.fetch_by_transcript( ivtype=ivtype) :
			
			row_id = record[0].attr['gene_id']

			if row_id == current_id :
				
				current_records.append( record )

			else : 
				
				if current_id : 
					yield merge_method[ method ]( current_records )

				current_records = [ record ]
				current_id = row_id
				
		yield current_records

	"""
	Interval type methods
	"""
	
	def get_exons( self, records ) :
		return self.get_type( records, 'exon' )

	def get_introns( self, records ) :
		introns = self.get_joining_intervals( self.get_exons(records) )
		for i in introns : i.type = 'intron'
		return introns

	def get_5utr( self, records ) :
		
		strand = records[0].iv.strand
		cds =  self.get_cds( records )
		if len(cds) == 0 : return None
		
		utr = []
		if strand == '+' :
			cds_start = min( [iv.iv.start for iv in cds] )
			utr = self.clip_leading_intervals( self.get_exons(records), cds_start )

		else :
			cds_end = max( iv.iv.end for iv in cds )
			utr = self.clip_trailing_intervals( self.get_exons(records), cds_end )
		
		for u in utr : u.type = '5UTR'
		return utr

	def get_3utr( self, records ) :

		strand = records[0].iv.strand
		cds = self.get_cds( records )
		if len(cds) == 0 : return None
	
		utr = []
		if strand == '+' :
			cds_end = max( [iv.iv.end for iv in cds] )
			utr = self.clip_trailing_intervals( self.get_exons(records), cds_end )

		else :
			cds_start = min( [iv.iv.start for iv in cds] )
			utr = self.clip_leading_intervals( self.get_exons(records), cds_start )

		for u in utr : u.type = '3UTR'
		return utr

	def get_cds( self, records ) :
		return self.get_type( records, 'CDS' )
		
	def get_union2gene( self, records ) :
		return self.get_type( records, 'union2gene' )

	def get_type( self, records, ivtype ) :
		return [ r for r in records if r.type == ivtype ]

	"""
	Interval joining methods
	"""

	def get_longest( self, records ) :
	
		max_len = 0
		max_record = None

		for record in records :

			rec_len = sum( [ r.iv.end - r.iv.start for r in record ] )

			if rec_len > max_len : 
				max_len = rec_len
				max_record = record

		return max_record

	def get_union( self, records ) :

		"""
		Return a list of intervals reflecting the union of all transcripts for a given gene.
		"""

		# collapses records into a single list rathern than list of lists.
		all_records = sum( records, [] )

		keyfun = lambda x : x.iv.start
		ivs = sorted( all_records, key=keyfun )
		
		union = [ [ ivs[0].iv.start, ivs[0].iv.end ] ]
		
		for iv in ivs :
			
			if iv.iv.start <= union[-1][1] :
				union[-1][1] = iv.iv.end
			else : 
				union.append( [iv.iv.start, iv.iv.end] )

		genomic_ivs = []
		for iv in union :
			genomic_iv = copy.deepcopy( ivs[0] )
			genomic_iv.iv.start = iv[0]
			genomic_iv.iv.end = iv[1]
			genomic_ivs.append( genomic_iv )

		return genomic_ivs


	def get_joining_intervals( self, intervals ) :

		"""
		Calcluate a new set of GenomicIntervals derived from those in between the intervals provided.
		"""

		keyfun = lambda x : x.iv.start
		ivs = sorted( intervals, key=keyfun )
		between = []

		for i in range(0, len(ivs)-1 ) :
			gi = copy.deepcopy( intervals[0] )
			gi.iv.end = ivs[i+1].iv.start
			gi.iv.start = ivs[i].iv.end

			between.append( gi )
		
		return between


	def clip_leading_intervals( self, intervals, max_position ) :

		"""
		Return intervals preceding a defined position (eg. CDS start) in genomic coordinates.

		Parameters: 

		intervals: list of HTSeq.GenomicFeature objects (note that the interval info is attr 'iv'
		max_position: genomic coordinate of position to clip at.

		"""

		keyfun = lambda x : x.iv.start
		ivs = sorted( intervals, key=keyfun )
		new = []

		for iv in ivs : 
			
			if iv.iv.end < max_position : 
				new.append( iv )
			else :
				
				if iv.iv.start < max_position :
					new_iv = copy.deepcopy( iv )
					new_iv.iv.end = max_position
					new.append( new_iv )
				
				break

		return new


	def clip_trailing_intervals( self, intervals, min_position ) :
	
		"""
		Return intervals trailing a defined position (eg. CDS stop) in genomic coordinates
		"""

		keyfun = lambda x : x.iv.end
		ivs = sorted( intervals, key=keyfun, reverse=True )
		new = []

		for iv in ivs :
			
			if iv.iv.start > min_position :
				new.append( iv )
			else :
				
				if iv.iv.end > min_position :
					new_iv = copy.deepcopy( iv )
					new_iv.iv.start = min_position
					new.append( new_iv )

				break

		return new


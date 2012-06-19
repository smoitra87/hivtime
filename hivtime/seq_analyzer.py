"""
	visualizes time series patient sequences
"""

import sys,os
from operator import itemgetter
import copy,glob

import numpy as np
import pylab as pl
import sqlite3 as lite
import Bio
from Bio import AlignIO
from Bio import SeqIO
from Bio import Motif
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.Align import AlignInfo
import markup
import mechanize
import urlparse

import utility as util

#----------------------------------------------------------------------
# Constants and definitions
dbpath = 'hivtime.db' 
datadir = util.get_datadir()
JALVIEW_DIR = "/home/subhodeep/Jalview"
LOSDB_URL = "http://www.hiv.lanl.gov/components/sequence/HIV/search/search.html"

dbname = 'hivtime.db'
dbpath = os.path.join(util.get_datadir(),dbname)
#tbls = ['p24.tbl','full_genome.tbl']
tbls = ['p24.tbl']
tbls = [os.path.join(util.get_datadir(),tbl) for tbl in tbls]

# If FETCH_ALN is True, then the alignment is downloaded anyway
FETCH_ALN = True
EXEC_JALVIEW = True

#----------------------------------------------------------------------
# Scripts

def get_tblname_from_fpath(fpath) : 
	""" Get the names of the tables """
	fname = os.path.split(fpath)[1]
	tblname = os.path.splitext(fname)[0]
	return tblname

class SeqAnalyzer(object)  :
	""" Performs analysis of sequences"""
	def __init__(self,cur,con,tblname) : 
		self.cur = cur
		self.con = con
		self.name = tblname
		self.header = map(itemgetter(1),self._get_header_meta())
		self.pcodes = self._get_patientcodes()
		self.time_fields = ('SamplingYear','DaysfromfirstSample',\
			'FiebigStage','DaysfromInfection','DaysfromSeroconversion')
		self.pcodes_time = None
		self.pat_obj = None
		self.page = None

	def analyze(self) : 
		""" Analyze all the sequences """
		print("Currently Analyzing %s"%(self.name))
	
		# Set up reporting
		page = markup.page()
		self.page = page
		title = "Index for Visualizations"
		page.h2(title)
		page.init(title=title)

		#--------------------------------------------------------------	
		# Find patients who have enough time data
		self.pcodes_time = self._find_time_pids()

		
		#--------------------------------------------------------------
		# For every patient create an alignment object and process 

#		avl_files = glob.glob(os.path.join(datadir,"*gag.fasta"))
#		avl_files = [os.path.split(f)[1] for f in avl_files]
#		pcodes_avl = [f.split('_gag.fasta')[0] for f in avl_files]
#
#		for pcode in pcodes_avl : 
#			self._process_pcode(pcode)

		for pcode in self.pcodes_time.keys() : 
			self._process_pcode(pcode)

		# Fill in the index elements
		page.hr()
		page.h3("Full Capsid Sequence Alignment")
		for pcode in self.pcodes_time : 
			page.a("%s p24 Alignment"%(pcode),\
				href="../data/images/%s_seq.html"%(pcode))	
			page.br()

		page.hr()
		page.h3("Full Capsid Changes Visualization Alignment")
		for pcode in self.pcodes_time : 
			page.a("%s p24 Visualization"%(pcode),\
				href="../data/images/%s_viz.html"%(pcode))	
			page.br()

		with open("index_viz.html","w") as fout :
			fout.write(page.__str__())

	
	def _process_pcode(self,pcode) : 
		""" Process every patient code"""
		self.pat_obj = PatientAlignment(pcode)
		self.pat_obj.write_aln()
		if EXEC_JALVIEW :
			self.pat_obj.print_jalview()
 
		
	def _get_header_meta(self) : 
		""" Get the header metadata """	
		cur.execute("PRAGMA table_info(%s)"%self.name)
		return cur.fetchall()

	def _get_patientcodes(self)  :
		""" Get all the patientcodes"""
		cur.execute("SELECT DISTINCT patientcode from %s\
			WHERE patientcode is NOT NULL"%(self.name))
		return map(itemgetter(0),cur.fetchall())
		

	def _find_time_pids(self) : 
		"""
Read the db to get patient ids of patients that have reasonable amount
of time series information

# Criteria to use ??
Patients who have atleast 4 time points any field of the 5
discriminating fields

What to do if multiple time fields have discriminating information ?
Pick timeforsero > timeinfec > dayssample > fiebig > sampling year
		"""
		cur = self.cur
		# REMOVEME : Going to work with n_filter = 11 items first for 
		# debug stage
		

		# n_filter=4 retreives about 60 patients
		n_filter = 4
		cmd_str = "DROP TABLE IF EXISTS dummy1"
		cur.execute(cmd_str)
		field_codes = 'ABCDE'
		field_code_dict = dict(zip(field_codes,self.time_fields))
		field_code_dict2 = dict(zip(self.time_fields,field_codes))

		cmd_str1 = ",".join(field_codes)
		cmd_str2 = ", ".join(["COUNT(DISTINCT %s) AS %s"%(f[0],f[1])\
			for f in zip(self.time_fields,field_codes)])
		cmd_str3 = " OR ".join(["%s>="%(f)+"%d" for f in field_codes])	
		cmd_str3 = cmd_str3%((n_filter,)*len(field_codes))
		cmd_str = "CREATE TABLE dummy1 AS SELECT patientcode," + \
			cmd_str1 + " FROM (SELECT patientcode," + cmd_str2 + \
			" FROM %s GROUP BY patientcode) WHERE "%(self.name) + \
			cmd_str3 + " ;" 

		# Example command from above construct
#		cmd_str = "CREATE TABLE dummy1 AS SELECT patientcode,A,B,C,D,E \
#FROM (SELECT \
# patientcode,COUNT(DISTINCT daysFROMfirstsample) AS A, \
#COUNT(DISTINCT daysFROMinfection) AS B, COUNT(DISTINCT \
#daysFROMseroconversion) AS C,COUNT(DISTINCT samplingyear) AS D,\
#COUNT(DISTINCT fiebigstage) AS E FROM p24 GROUP BY patientcode) \
#where A>=%d OR B>=%d OR C>=%d OR D>=%d OR E>=%d ;"%\
#			((n_filter,)*5)
		
		cur.execute(cmd_str)
		cmd_str = "SELECT count(patientcode) from dummy1"
		cur.execute(cmd_str)
		count_pats = cur.fetchone()[0]
		print("The number of patients satisfying criteria %d = %d"%\
				(n_filter,count_pats))
			
		# Find the fields which have the excess time points as well
		# Patient codes satisfying the filter constraints

		# Dicts of patientcodes and their correpsonding fields
		# This is returned from func 
		pcodes_time = {}

		#Priority of time fields
		priority_time_fields = ('DaysfromSeroconversion',\
			'DaysfromInfection','DaysfromfirstSample','SamplingYear',\
			'FiebigStage')

		priority_field_codes = [field_code_dict2[f] for f in \
			priority_time_fields]
	
		# Loop through all the fields that have the information
		for field in priority_field_codes : 
			cmd_str = "SELECT patientcode,%s FROM dummy1  where %s \
				>= %d"%(field,field,n_filter)
			cur.execute(cmd_str)
			records = cur.fetchall()
			pcodes = map(itemgetter(0),records)
			for pcode in pcodes :
				if not pcodes_time.has_key(pcode) : 
					pcodes_time[pcode] = field_code_dict[field]

		return pcodes_time

class PatientAlignment(object) :
	"""Creates an alignment of patient sequences """
	def __init__(self,pat_code) :
		self.pcode = pat_code
		self.FLAG_FETCH_ALN = FETCH_ALN
		self.aln = None
		self.fname = self.pcode+'_gag.fasta'
		self.fpath = os.path.join(datadir,self.fname)
		self._p24_start = 'PIVQN'
		self._p24_end = 'KARVL'
		self._load_aln()
		self.consensus = self._get_consensus()
		self._create_viz_aln()

	def _load_aln(self) :
		""" Load the alignment 
Try to see if the fasta file is already present in the data directory
otherwise force a download of the alignment using browser automation
"""
		if self.FLAG_FETCH_ALN: 
			self._fetch_aln()
		else :
			if not self._check_aln_exists() : 
				print("Did not find %s"%self.fpath)
				self.aln = self._fetch_aln()
	
		# Fix the fasta file first
		self._fix_fasta()

		# Read in the fasta file as a biopython alignment object
		self.aln = AlignIO.read(open(self.fpath),"fasta")
		print "Gag Alignment Lenth is %d"% \
			self.aln.get_alignment_length()

		# Find the start and end of the p24 region
		m = Motif.Motif(alphabet=IUPAC.protein)
		m.add_instance(Seq(self._p24_start,m.alphabet))
		m.add_instance(Seq(self._p24_end,m.alphabet))

		hxb2_seq = self.aln[0].seq # hxb2 seq is always on top
		
		search_res = list(m.search_instances(hxb2_seq))
		if len(search_res) != 2 : 
			raise ValueError("Could not find the motifs for pat"\
				%self.pcode)

		# Cut the p24 region
		p24_startpos = search_res[0][0]
		p24_endpos = search_res[1][0]
		self.aln = self.aln[:, p24_startpos:p24_endpos+\
			len(self._p24_end)]
		print "Capsid alignment length is %d"%\
			(self.aln.get_alignment_length())

	def _get_consensus(self) :
		""" Find the consensus sequence for the alignment """
		aln_info = AlignInfo.SummaryInfo(self.aln)
		return aln_info.gap_consensus(threshold=0.7)

	def _create_viz_aln(self) : 
		""" Create an alignment for visualization purposes showing only those positions which are changing """
		aln_viz = copy.deepcopy(self.aln)
		
		def f1(zip_tuple) :
			""" Compare two aa. If they are the same then replace by .
else return aa_aln
 """		
			aa_aln,aa_cons = zip_tuple
			if aa_aln == aa_cons :
				return '.'
			else : 
				return aa_aln
						
		for rec in aln_viz : 
			seq_str = ''.join(map(f1,zip(rec,self.consensus)))		
			rec.seq = Seq(seq_str)

		fname = self.pcode + '_viz.fasta'
		fpath = os.path.join(datadir,fname)
		with open(fpath,"w") as fout : 
			AlignIO.write([aln_viz],fout,"fasta")	

		return aln_viz

	def write_aln(self,suffix='p24') :
		""" Write out the current alignment """
		fname = self.pcode + '_' + suffix + ".fasta"
		fpath = os.path.join(datadir,fname)
		with open(fpath,'w') as fout :
			AlignIO.write([self.aln],fout,"fasta")

	def print_jalview(self) : 
		""" Print out an html page containing images"""
		cur_dir = os.getcwd()
	
		# Issue command to execute Jalview
		os.chdir(JALVIEW_DIR)
		cmd_str = './Jalview -nodisplay -open ~/projects/hivtime/data/%s_viz.fasta -colour zappo -png ~/projects/hivtime/data/images/%s_viz.png -imgMap ~/projects/hivtime/data/images/%s_viz.html'%((self.pcode,)*3)
		os.system(cmd_str)
	
		cmd_str = './Jalview -nodisplay -open ~/projects/hivtime/data/%s_p24.fasta -colour zappo -png ~/projects/hivtime/data/images/%s_seq.png -imgMap ~/projects/hivtime/data/images/%s_seq.html'%((self.pcode,)*3)
		os.system(cmd_str)
		
		os.chdir(cur_dir)
	
	def _fix_fasta(self) :
		""" Read in the fasta file and make some fixes 

Fix 1 : Some of the traling characters gaps are missing. Add them. 


"""
		records = list(SeqIO.parse(self.fpath,"fasta"))
		max_len = max([len(rec) for rec in records])
		for rec in records : 
			if len(rec) < max_len : 
				rec.seq = rec.seq + '-'*(max_len-len(rec))

		SeqIO.write(records,self.fpath,format="fasta")

	def _fetch_aln(self) : 
		""" Download the alignment using a browser session """
		print("Downloading %s via Browser session..."\
				%(self.pcode))

		# Start the broser
		br = mechanize.Browser()	
		
		# Set the broser options
		br.set_handle_robots(False)

		# Open the search page
		br.open(LOSDB_URL)
		br.select_form(nr=1)
		br["value SEQ_SAMple SSAM_postInfect_days 4"] = '*'
		br["value SEQ_SAMple SSAM_postSeroConv_days 4"] = '*'
		br["value PATient PAT_code 1 exact"] = self.pcode	
		br['Genomic Region'] = ['P24']
		br.submit()
		br.select_form(nr=1)
		br.form["INCLUDE_HXB2"] = ['1']
		br.form["translate"] = ['FALSE_AA']
		br.submit(name='save')
		link = br.find_link(text_regex="gag.fasta")
		full_url = urlparse.urljoin(link.base_url,link.url)
		br.retrieve(full_url,os.path.join(datadir,\
			self.pcode+'_gag.fasta'))
		br.close()

	def _check_aln_exists(self) :
		""" Check if the patient file exists """
		return os.path.isfile(self.fpath)
	
	def _cut_aln(self) : 
		""" Cut the alignment according to p24 dimensions """
		pass

#-----------------------------------------------------------------------
# Main Script
#-----------------------------------------------------------------------

with lite.connect(dbpath) as con : 
	# Check SQL version for sanity check
	cur = con.cursor()
	cur.execute("SELECT SQLITE_VERSION()")
	print "SQL Version %s"%cur.fetchone()

	# Extract tablenames from full path	
	tblnames = map(get_tblname_from_fpath,tbls)
		
	for tblname in tblnames : 
		seq_analyzer =  SeqAnalyzer(cur,con,tblname);
		seq_analyzer.analyze()

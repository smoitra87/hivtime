"""
	visualizes time series patient sequences
"""

import sys,os,pdb
from operator import itemgetter

import numpy as np
import pylab as pl
import sqlite3 as lite
import Bio

import utility as util

#----------------------------------------------------------------------
# Constants and definitions
dbpath = 'hivtime.db' 
datadir = util.get_datadir()

dbname = 'hivtime.db'
dbpath = os.path.join(util.get_datadir(),dbname)
#tbls = ['p24.tbl','full_genome.tbl']
tbls = ['p24.tbl']
tbls = [os.path.join(util.get_datadir(),tbl) for tbl in tbls]

# If FETCH_ALN is True, then the alignment is downloaded anyway
FETCH_ALN = False

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

	def analyze(self) : 
		""" Analyze all the sequences """
		print("Currently Analyzing %s"%(self.name))
		
		# Find patients who have enough time data
		self.pcodes_time = self._find_time_pids() 
		
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
	""" """
	def __init__(self,pat_code) :
		self.pcode = pat_code
		self.FLAG_FETCH_ALN = FETCH_ALN
		self.aln = None

	def _load_aln(self) :
		""" Load the alignment 
Try to see if the fasta file is already present in the data directory
otherwise force a download of the alignment using browser automation
"""
		if self.FLAG_FETCH_ALN(self) : 
			self._fetch_aln()
		else :
			if not self._check_aln_exists() : 
				self.aln = self._fetch_aln()
		
		self.aln = self._read_aln() 

	def _fetch_aln(self) : 
		""" Download the alignment using a browser session """
		pass


	def _check_aln_exists(self) :
		""" Check if the patient file exists """
		fname = self.pcode+'_gag.fasta'
		return os.path.isfile(fname)
	
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

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
tbls = ['full_genome.tbl']
tbls = [os.path.join(util.get_datadir(),tbl) for tbl in tbls]



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

	def analyze(self) : 
		""" Analyze all the sequences """
		print("Currently Analyzing %s"%(self.name))
		
		# Find patients who have enough time data
		self.find_time_pids() 


	def _get_header_meta(self) : 
		""" Get the header metadata """	
		cur.execute("PRAGMA table_info(%s)"%self.name)
		return cur.fetchall()

	def _get_patientcodes(self)  :
		""" Get all the patientcodes"""
		cur.execute("SELECT DISTINCT patientcode from %s\
			WHERE patientcode is NOT NULL"%(self.name))
		return map(itemgetter(0),cur.fetchall())
		

	def find_time_pids(self) : 
		"""
Read the db to get patient ids of patients that have reasonable amount
of time series information

# Criteria to use ??
Patients who have atleast 4 time points any field of the 5 discriminating fields
		"""
		cur = self.cur
		n_filter = 2
		for iter_filter in range(n_filter,n_filter+15) : 
			cmd_str = "DROP TABLE IF EXISTS dummy1"
			cur.execute(cmd_str)
			cmd_str = "CREATE TABLE dummy1 AS SELECT patientcode,A,B,C,D,E \
	FROM (SELECT \
	 patientcode,COUNT(DISTINCT daysFROMfirstsample) AS A, \
	COUNT(DISTINCT daysFROMinfection) AS B, COUNT(DISTINCT \
	daysFROMseroconversion) AS C,COUNT(DISTINCT samplingyear) AS D,\
	COUNT(DISTINCT fiebigstage) AS E FROM p24 GROUP BY patientcode) \
	where A>=%d OR B>=%d OR C>=%d OR D>=%d OR E>=%d ;"%\
				((iter_filter,)*5)
			cur.execute(cmd_str)
			cmd_str = "SELECT count(patientcode) from dummy1"
			cur.execute(cmd_str)
			count_pats = cur.fetchone()[0]
			print("The number of patients satisfying %d criteria=%d"%\
					(iter_filter,count_pats))
		


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

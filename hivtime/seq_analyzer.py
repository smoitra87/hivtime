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
		"""
		cmd_str = "SELECT * from "

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

""" 
Performs Data Analysis for HIV Time varying sequences
	
	@author: Subhodeep Moitra (smoitra@cs.cmu.edu)
"""

import sqlite3 as lite
import os,sys
import utility as util
from utility import get_datadir

#----------------------------------------------------------------------
# Params --------------------------------------------------------------
#----------------------------------------------------------------------

# Filenames and paths
# Only one db but multiple tables
dbname = 'hivtime.db'
dbpath = os.path.join(get_datadir(),dbname)
tbls = ['p24.tbl','full_genome.tbl']
#tbls = ['p24.tbl']
tbls = [os.path.join(get_datadir(),tbl) for tbl in tbls]

# Boolean params
TABLE_CREATE=False
DATA_ANALYZE=True

#----------------------------------------------------------------------
def parse_table(fpath) : 
	""" Parse the table """
	with open(fpath,'r') as fin :
		tbldat = fin.readlines()
	
	# Skip preheader lines 
	while tbldat[0].strip() == '' : 
		tbldat = tbldat[1:]
	while tbldat[0].startswith('Number of records') : 
		tbldat = tbldat[1:]

	# Get the header and the body 
	header,records = tbldat[0],tbldat[1:]
	header = header.strip().split('\t')
	header = [h.replace(' ','') for h in header]
	records = [record.strip().split('\t') for record in records]

	return header,records
	

def create_table(con,cur,fpath) :
	# Get the name of the table from the file name
	fname = os.path.split(fpath)[1]
	tblname = os.path.splitext(fname)[0]
	cur.execute('DROP TABLE IF EXISTS %s'%(tblname))	
	header,records = parse_table(fpath)
	code_tbl = dict(util.code_tbl)

	# Create table with header
	cmd = 'CREATE TABLE %s('%(tblname)	
	header_t = ['\''+h+'\' '+code_tbl[h] for h in header]
	cmd += ','.join(header_t) + ')'
	cur.execute(cmd)

	# Populate the table
	for record in records :
		values = ','.join(['\''+r+'\'' for r in record])
		cur.execute("INSERT INTO %s VALUES (%s)"%(tblname,values))
		

def get_tblname_from_fpath(fpath) : 
	""" Get the names of the tables """
	fname = os.path.split(fpath)[1]
	tblname = os.path.splitext(fname)[0]
	return tblname


class TableAnalyzer(object) : 
	""" Performs analysis on tables """	
	def __init__(self,cur,tblname): 
		self.cur = cur # The database cursor handles
		self.name = tblname
		self.header_meta = self._get_header_meta()
	
	def _get_header_meta(self) : 
		""" Get the header metadata """	
		cur.execute("PRAGMA table_info(%s)"%self.name)
		return cur.fetchall()

	def analyze(self)	: 
		""" Perform the analysis """
		print "Analyzing %s"%(self.name)
		
#----------------------------------------------------------------------
#  Main Script
#----------------------------------------------------------------------

# Connect to the database
with lite.connect(dbpath) as con : 
	# Check SQL version for sanity check
	cur = con.cursor()
	cur.execute("SELECT SQLITE_VERSION()")
	print "SQL Version %s"%cur.fetchone()
	
	# Overwrite previous tables id TABLE_CREATE is set to true
	if TABLE_CREATE : 
		for tblname in tbls : 
			create_table(con,cur,tblname)

	if DATA_ANALYZE : 
		tblnames = map(get_tblname_from_fpath,tbls)
		for tblname in tblnames : 
			tbl_analyzer = TableAnalyzer(cur,tblname)
			tbl_analyzer.analyze()
	
					

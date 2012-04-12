""" 
Performs Data Analysis for HIV Time varying sequences
	
	@author: Subhodeep Moitra (smoitra@cs.cmu.edu)
"""

import sqlite3 as lite
import os,sys
import utility as util
from utility import get_datadir
import markup
from operator import itemgetter

import numpy as np
import pylab as pl



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
PLOT_AGAIN=False


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
	header = [h.replace('#','_n_') for h in header]
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
	header_t = [h+' '+code_tbl[h] for h in header]
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


def make_table(page,rows,header=None) : 
	""" Make an HTML table for a markup page """
	page.table(border=1)
	page.tr()
	page.th(header)
	page.tr.close()
	for row in rows :
		page.tr()
		page.td(row)
		page.tr.close()
	page.table.close()


class TableAnalyzer(object) : 
	""" Performs analysis on tables """	
	def __init__(self,cur,tblname): 
		self.cur = cur # The database cursor handles
		self.name = tblname
		self.header_meta = self._get_header_meta()
		self.header = [h[1] for h in self.header_meta ]
		self.page = None
	
	def _get_header_meta(self) : 
		""" Get the header metadata """	
		cur.execute("PRAGMA table_info(%s)"%self.name)
		return cur.fetchall()

	def analyze(self)	: 
		""" Perform the analysis """
		print "Analyzing %s"%(self.name)
		code_tbl = dict(util.code_tbl)
	
		# Set up reporting
		page = markup.page()
		self.page = page
		title="Analysis for %s"%(self.name)
		page.h2(title)
		page.init(title=title)

		#----------------------------------------------------------		
		page.h4("How many non null elements in every Field?")	
		page.ul(class_='mylist')
		for h in self.header :
			cur.execute("SELECT COUNT(*) from %s where %s !=''"%\
				(self.name,h))
			val = cur.fetchone()
			page.li(h+' : '+str(val[0]),class_='myitem')
		page.ul.close()
		
		#----------------------------------------------------------
		page.h4("How many non null distinct elements in every field?")
		page.ul(class_='mylist')
		for h in self.header :
			cur.execute(""" SELECT COUNT(*) from (
					SELECT DISTINCT %s from %s where %s !='')
				"""%(h,self.name,h))
			val = cur.fetchone()
			page.li(h+' : '+str(val[0]),class_='myitem')
		page.ul.close()


		#---------------------------------------------------------
		page.h4("Analysis of all Numeric fields")

		for h in self.header :
			if code_tbl[h] not in ('INT','FLOAT') : 
				continue
			else : 

				page.h5("Numerical analysis of %s"%h)
				cur.execute("SELECT %s from %s where %s != ''"%\
					(h,self.name,h))
				records = cur.fetchall()
				vals = map(itemgetter(0),records)
				
				# Special cases
				if h == u'DaysfromSeroconversion' :
					# All unicode characters are cast as -1
					page.p("Warning..!!! early, late and \
						pre-seroconversion are special values")
			
					nearly,nlate,npre=0,0,0
					for ii,val in enumerate(vals) : 
						if val == u'early' : 
							vals[ii] = 30
							nearly += 1
						if val == u'late' :
							vals[ii] = 1000
							nlate += 1
						if val == u'pre-seroconversion' :
							vals[ii] = -10
							npre += 1
					page.p("Using early=30,late=1000,pre-sero=-10")
					page.p("n_early=%d,nlate=%d,npre=%d"%\
						(nearly,nlate,npre))

				vals = np.array(map(int,vals))
				page.p("Number of elements in field : %s"%(len(vals)))			
	
				if len(vals) == 0  : 
					page.p("No elements found..!")
					continue
	
				page.ul(Class="mylist")
				page.li("Max %s : %f"%(h,vals.max()),class_='myitem')
				page.li("Min %s : %f"%(h,vals.min()),class_='myitem')
				page.li("Avg %s : %f"%(h,vals.mean()),class_='myitem')
				page.li("Std %s : %f"%(h,vals.std()),class_='myitem')
				page.ul.close()

				# Plot histograms
				figpath = "figs/"+self.name+"_"+h+".png"
				if PLOT_AGAIN :
					pl.figure()
					pl.hist(vals,40,facecolor='green')
					pl.savefig(figpath)	
		
				page.img(width=400,height=300,alt="HistPlots",\
					src=figpath)

		#------------------------------------------------------------
		# Write out results
		with open(self.name+"_report.html",'w') as fout : 
			fout.write(page.__str__())
		
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

		# create the index page 
		page = markup.page()
		page.init(title="Index page")		
		page.h2("Index")
		for tblname in tblnames : 
			page.a(tblname+" Report",href=tblname+"_report.html")
			page.br()
		with open("index.html",'w')	 as fout : 
			fout.write(page.__str__())

	
		
		for tblname in tblnames : 
			tbl_analyzer = TableAnalyzer(cur,tblname)
			tbl_analyzer.analyze()

	
					

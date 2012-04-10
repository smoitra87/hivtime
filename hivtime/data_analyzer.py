""" 
Performs Data Analysis for HIV Time varying sequences
	
	@author: Subhodeep Moitra (smoitra@cs.cmu.edu)
"""

import sqlite3 as lite
import os,sys
from utility import get_datadir



# Only one db but multiple tables
dbname = 'hivtime.db'

class DBHandler(object) : 
	""" Deals with database """
	def __init__(self,dbname) :
		""" Intiializer """
		self.dbname = dbname
		datadir = get_datadir()
		self.dbpath = os.path.join(datadir,dbname)

	def _check_db(self) :
		""" Check if the db exists """
		return os.path.exists(self.dbpath) 

	def _create_db(self) : 
		""" Creates and populates db from flat files """
		pass
	
		

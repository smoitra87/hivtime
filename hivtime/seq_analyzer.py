"""
	visualizes time series patient sequences

"""


import sys,os,pdb

import numpy as np
import pylab as pl
import sqlite3 as lite
import Bio

import utility as util



class SeqAnalyzer(object)  :
	""" Performs analysis of sequences"""
	def __init__(self) : 
		self.db_name = ""

	def analyze(self) : 
		""" Analyze all the sequences """
		pass

	def find_pids(self) : 
		""" Read the db"""
		pass	


if __name__ == '__main__' : 
	seq_analyzer =  SeqAnalyzer();
	seq_analyzer.analyze()

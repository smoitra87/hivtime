""" 
This script provides utility functions
"""

import os,sys

code_tbl = [
	('#','INT'),
	('PatientCode','TEXT'),
	('PatientId','TEXT'),
	('Accession','TEXT'),
	('Name','TEXT'),
	('Subtype','TEXT'),
	('Country','TEXT'),
	('SamplingYear','INT'),
	('DaysfromfirstSample','INT'),
	('FiebigStage','TEXT'),
	('Daysfromtreatmentend','INT'),
	('#ofpatienttimepoints','INT'),
	('Daysfromtreatmentstart','INT'),
	('DaysfromInfection','INT'),
	('DaysfromSeroconversion','INT'),
	('PubmedID','TEXT'),
	('Start','INT'),
	('Stop','INT'),
	('SequenceLength','INT'),
	('Organism','TEXT')
]

def get_basedir() :
	""" Get the base directory of the project"""
	paths = os.path.split
	pathr = os.path.realpath
	pathd = os.path.dirname
	pathj = os.path.join
	pathe = os.path.exists
	scriptname = sys.argv[0]
	scriptdir = pathr(pathd(scriptname))
	basedir = paths(scriptdir)[0]
	return basedir

def get_homedir() :
	""" Get the home directory of the user"""
	homedir = os.getenv('HOME')
	return homedir

def get_datadir() : 
	""" Get the Data directory of the project  """
	basedir = get_basedir()
	datadir = os.path.join(basedir,'data')
	return datadir

basedir = get_basedir()



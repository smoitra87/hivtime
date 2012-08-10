""" 
This script provides utility functions
"""

import os,sys

code_tbl = [
	('_n_','INT'),
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
	('_n_ofpatienttimepoints','INT'),
	('Daysfromtreatmentstart','INT'),
	('DaysfromInfection','INT'),
	('DaysfromSeroconversion','INT'),
	('PubmedID','TEXT'),
	('Start','INT'),
	('Stop','INT'),
	('SequenceLength','INT'),
	('Organism','TEXT')
]

"""
 Maps from aa type to a tuple containing polarity and charge information 
 (0/1,-1/0/1) - NonPolar/Polar, Negative,Neutral,Positive
"""
aa_prop_tbl =  {
	'A' : (0,0),
	'G' : (0,0),
	'L' : (0,0),
	'I' : (0,0),
	'M' : (0,0),
	'F' : (0,0),
	'P' : (0,0),
	'W' : (0,0),
	'V' : (0,0),
	'N' : (1,0),
	'C' : (1,0),
	'Q' : (1,0),
	'H' : (1,0),
	'S' : (1,0),
	'T' : (1,0),
	'Y' : (1,0),
	'R' : (1,1),
	'K' : (1,1),
	'D' : (1,-1),
	'E' : (1,-1),
}

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



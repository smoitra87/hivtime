"""
	visualizes time series patient sequences
"""

import sys,os
from operator import itemgetter
import copy,glob
import time
import itertools


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
import pymol
from pymol import cmd

import utility as util

#----------------------------------------------------------------------
# Constants and definitions
dbpath = 'hivtime.db' 
datadir = util.get_datadir()
JALVIEW_DIR = "/home/subhodeep/Jalview"
LOSDB_URL = "http://www.hiv.lanl.gov/components/sequence/HIV/search/search.html"

dbname = 'hivtime.db'
dbpath = os.path.join(util.get_datadir(),dbname)
pdbname="3H4E.pdb"
pdbpath = os.path.join(datadir,pdbname)
edgefname="31-Aug-2011-edges.dat"
edgepath = os.path.join(datadir,edgefname)
#tbls = ['p24.tbl','full_genome.tbl']
tbls = ['p24.tbl']
tbls = [os.path.join(util.get_datadir(),tbl) for tbl in tbls]
ent_top_cutoff = 0.03 # used for annotation top entropy values
ent_stat_cutoff = 0.4 # used for speeding up stat calculations

# If FETCH_ALN is True, then the alignment is downloaded anyway
FETCH_ALN = False
EXEC_JALVIEW = False
DO_PROCESS = False
DO_HOTSPOT = False
DO_STAT = True

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

		for pcode in self.pcodes_time.keys() : 
			if DO_PROCESS :
				self._process_pcode(pcode)

		#--------------------------------------------------------------
		# Create a hotspot visualization based on the entropy scores of
		# all sequences pooled together

		if DO_HOTSPOT : 
			print("Beginning hotspot analysis")
			pat_ent_aa = []
			for pcode in self.pcodes_time.keys() :
				p24fpath = os.path.join(datadir,pcode+'_p24.fasta')
				if not os.path.exists(p24fpath) : 
					raise IOError('File %s not found'%p24fpath)
				aln = AlignIO.read(open(p24fpath),"fasta")
				ent_aa = self._entropy_aa(aln)
				pat_ent_aa.append(ent_aa)
			
			# Calculate the mean of the entropy values
			# WARNING..!This is a not a rigorous technique
			avg_ent = np.mean(pat_ent_aa,axis=0)
			
			# Visualize the entropy values and plot it
			pl.figure()
			pl.plot(avg_ent)
			pl.title("Plot of the average entropy across patients")
			pl.xlabel('AA position')
			pl.ylabel('Entropy Values')

			# Get the canonical sequence to get the amino acids
			self._canon_seq = AlignIO.read(open(p24fpath),"fasta")[0]
			self._canon_seq = self._canon_seq.seq

			# Place text on peaks
			idx_peaks = np.where(avg_ent > ent_top_cutoff)[0]			
			for idx in idx_peaks :
				pl.text(idx,avg_ent[idx]+0.002,str(idx+1)+\
					self._canon_seq[idx])
		
			pl.axhline(y=ent_top_cutoff,color="black",linewidth=2,\
				linestyle="dashed")
	
			pl.savefig("figs/spikeplot.png")
			pl.close()

			# Create Pymol visualization of hotspot
			self._draw_hotspot_pymol(avg_ent)

			# Create hotspot page
			self._create_hotspot_page()
	
		if DO_STAT : 
			print("Beginning Statistical analysis")
	
			# Exploratory analysis for finding distribution of entropy
			# per column per patient
			#self._explore_stats1()
			
			for pcode in self.pcodes_time.keys(): 
 				self._calc_stats(pcode)
					
	
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

	def _explore_stats1(self) : 
		""" Do exploratory statistical analysis for entropy per patient
		"""
		ent_pat_aa = []
		for pcode in self.pcodes_time.keys(): 
			p24fpath = os.path.join(datadir,pcode+'_p24.fasta')
			if not os.path.exists(p24fpath) : 
				raise IOError('File %s not found'%p24fpath)
			aln = AlignIO.read(open(p24fpath),"fasta")
			ent_aa = self._entropy_aa(aln)
			ent_pat_aa.extend(ent_aa)

		pl.hist(ent_pat_aa,bins=100,facecolor='green',log=True)
		pl.savefig("figs/explore_stat1.png")
		pl.close()
		

	def _calc_stats(self,pcode) : 
		""" 
Calculate statistics for each patient alignment
Questions answered by stat analysis : 

Q1. Does a residue position change a lot ?
	A1. Hotspot analysis. 

Q2. Does a residue position change in a patient with time ?
	A2. Autocorrelation/MI/Hyp tests

Q3. Is a residue position coupled with another residue position ? 
	A3. CHAVI gremlin analysis

Q4. is a residue position coupled with another residue per patient 
and with time ? 
	A4. Paired aa counts
		"""
		
		analysis = {}
		analysis['pcode'] = pcode;

		p24fpath = os.path.join(datadir,pcode+'_p24.fasta')
		if not os.path.exists(p24fpath) : 
			raise IOError('File %s not found'%p24fpath)
		aln = AlignIO.read(open(p24fpath),"fasta")
		aln = aln[1:] # Remove reference row
		ent_aa = np.array(self._entropy_aa(aln))

		# find indices of all columsn greater than cutoff
		idx_cols = np.where(ent_aa > ent_stat_cutoff)[0]
		analysis['candidate_cols'] = idx_cols;
	
		#print pcode," : ", idx_cols

		# Calculate columnwise statistics		
		times = [int(rec.id.split('.')[0]) for rec in aln ]
		uniq_times = sorted(set(times))		
		if len(uniq_times) > 1 : 
			scale_factor = len(uniq_times)-1 
		else : 
			scale_factor = 1
		
		# List of order means for each candidate column in patient
		patient_col_mean = []

		for col in idx_cols:
			col_aa = aln[:,col]
		 	uniq_aa = set(col_aa) # Unique aa in column
			col_aa_mean = {}
			for aa in uniq_aa : 
				# Find row ids or aa types
				row_idx = np.where(map(lambda s:s==aa,col_aa))[0]
				# Find corresponding times for aa types
				times_idx = filter(lambda x : x[0] in row_idx,\
					enumerate(times))
				times_idx = map(itemgetter(1),times_idx)
				uniq_times_aa = sorted(set(times_idx))
				# Get rank-ordering of occurence of aa type
				ranks_aa =	filter(lambda x:x[1] in uniq_times_aa,\
					enumerate(uniq_times))
				ranks_aa = map(itemgetter(0),ranks_aa)
				#scale ranks to be b/w 0 and 1
				ranks_aa_scaled = \
					map(lambda x: (0.0+x)/scale_factor,ranks_aa)
				aa_order_mean = np.mean(ranks_aa_scaled)
				col_aa_mean[aa] = aa_order_mean
			patient_col_mean.append(col_aa_mean)	
		
		if len(patient_col_mean) > 0  :
			self._plot_patient_order_means(patient_col_mean)
			
		# Calculate column pair statistics 
		for col1, col2 in itertools.combinations(idx_cols,2) : 
			pass
	def _plot_patient_order_means(self,patient_col_mean):
		""" Make a bar plot like figure of order means for every
		column in a  patient"""
		pl.figure()
		pl.axis([0,len(patient_col_mean),0,1])
		for i in range(len(patient_col_mean))[1:] :
			pl.axvline(x=i)
		pl.show()
		pl.close()	
	
	def _draw_hotspot_pymol(self,avg_ent) :
		""" Draws hotspot values on the Capsid structure"""
		pymol.finish_launching()
		cmd.load(pdbpath)	
		cmd.bg_color('white')
		cmd.hide('all')
		cmd.show('cartoon','3H4E and chain A')
		cmd.set_view(\
		    '-0.912057340, -0.295883119, -0.283901036,\
			0.399889141, -0.795006156, -0.456119299,\
			-0.090745762, -0.529538393, 0.843419492,\
			0.000087761, -0.000581503, -183.590560913,\
			8.321296692, -66.967765808, -33.590000153,\
			-36.799331665, 403.937011719, -20.000000000' )	
		
		# set the b-factors
		from pymol import stored
		stored.newB = []
		stored.newB.extend(avg_ent)
		
		# clear out the old B Factors
		cmd.alter('3H4E','b=0.0');

		# update the B Factors with new properties
		cmd.alter('3H4E and n. CA','b=stored.newB.pop(0)');
		cmd.spectrum("b",selection="3H4E and n. CA")	

		# Label the high entropy residues
		idx_peaks = np.where(avg_ent>ent_top_cutoff)[0]+1
		for resi in idx_peaks : 
			cmd.label('3H4E//A/'+str(resi)+'/CA','"%s%s"%(resn,resi)')

		cmd.ray()
		sys.stderr.write("Sleeping for 1 sec before saving image")
		time.sleep(1)
		figpath="figs/hot_capsid.png"
		cmd.png(figpath,width=800,height=600)
		
		#--------------------------------------------------------------
		# Now draw edges on the structure and save that image too
		cutoff=2

		with open(edgepath,'r') as fin :
		    data = fin.readlines();
		edges = [line.strip().split() for line in data]
		
		for n,e in enumerate(edges) :
		    t = cmd.dist('dist'+str(n),'3H4E//A/'+e[0]+'/CA',\
		'3H4E//A/'+e[1]+'/CA')
		
		    if t < cutoff :
		        cmd.delete('dist'+str(n))
		    else :
		        cmd.color('magenta','dist'+str(n))
		        cmd.hide('labels','dist'+str(n))
		        cmd.show('spheres','3H4E//A/'+e[0]+'/CA')
		        cmd.show('spheres','3H4E//A/'+e[1]+'/CA')
		        #cmd.color('yellow','3H4E//A/'+e[0]+'/CA')
		        #cmd.color('yellow','3H4E//A/'+e[1]+'/CA')
		
		cmd.set('dash_radius','0.2');
		cmd.set('dash_gap','0.0');
		cmd.set('sphere_scale','0.3')
		#cmd.set('sphere_scale','0.8','retinal');
		cmd.ray()
		figpath="figs/hot_capsid_gremlin.png"
		cmd.png(figpath,width=800,height=600)

		# Save the PyMol session
		sesspath=os.path.join(datadir,"hot_capsid.pse")
		cmd.save(sesspath)
		cmd.delete('all')
		cmd.quit()

	def _create_hotspot_page(self) : 
		""" Creates the Hotspot HTML page"""
		_page_prev = self.page # save the prev page
		self.page = markup.page()
		page = self.page
		title ="Hotspot visualization page"
		page.h2(title)
		page.init(title=title)
		
		page.h4("Hotspot Spike Plot")
		figpath="figs/spikeplot.png"
		page.img(width=800,height=600,alt="SpikePlot",src=figpath)

		page.hr()
		page.h4("Hotspots laid out on structure")
		page.a("Hotspot Pymol Session",href="../data/hot_capsid.pse")
		page.br()

		figpath="figs/hot_capsid.png"
		page.img(width=800,height=600,alt="Hotspots on Capsid",\
			src=figpath)

		page.br()
		page.h4("Hotspots laid out on structure with gremlin edges")
		figpath="figs/hot_capsid_gremlin.png"
		page.img(width=800,height=600,alt="Hotspots Capsid + gremlin",\
			src=figpath)

		with open("hotspot_viz.html","w") as fout : 
			fout.write(page.__str__())

		# Write out the hotspot page
		self.page = _page_prev # restore the previous page


	def _entropy_aa(self,aln) : 
		""" Returns a list containing the entropy in each aa position
		"""
		# Remove the first sequence as it is the reference sequence
		aln = aln[1:]	
		
		# Calculate the entropy for each position
		info_cols = [] # Non conserved columns
		for col_id in range(aln.get_alignment_length()) :
			col = aln.get_column(col_id)
			ftab = {} # Freq table
			for c in col :
				ftab[c] = ftab.get(c,0) + 1
			counts = ftab.values()
			sum_c = sum(counts)
			p = map(lambda x : (0.0+x)/sum_c,counts) # probabilities
			e = -reduce(lambda x,y:x+y,map(lambda t:t[0]*t[1],\
				zip(p,np.log(p))))
			info_cols.append(e)
		return info_cols
	

	
	def _process_pcode(self,pcode) : 
		""" Process every patient code"""
		self.pat_obj = PatientAlignment(pcode,self)
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
	def __init__(self,pat_code,seqobj) :
		self.pcode = pat_code
		self.seqobj = seqobj # Pointer to parent seqeunce analyzer obj
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
Fix 2 : Examine the accession ids from the fasta filed with those 
stored in the SQL table. There might be discrepancies because patient
codes are not unique. We should have used patientids instead. 
Fix 3 : Make sure that the sequences are sorted in ascending time
Fix 4 : Time information is inserted as first field in record name

"""
		records = list(SeqIO.parse(self.fpath,"fasta"))

		# Fix 1
		max_len = max([len(rec) for rec in records])
		for rec in records : 
			if len(rec) < max_len : 
				rec.seq = rec.seq + '-'*(max_len-len(rec))

		# Fix 2
		acc_names = [self._get_accname(rec) for rec in records[1:]]
		# query the sqlite database to get the list of accession names
		# for the patient
		cur.execute("SELECT DISTINCT accession from %s\
			WHERE patientcode='%s'"%(self.seqobj.name,self.pcode))
		acc_names_tbl = map(itemgetter(0),cur.fetchall())
		# sequences in fasta alignment but not in db search
		diff = set(acc_names).difference(acc_names_tbl)

		records = [rec for rec in records if\
			 self._get_accname(rec) not in diff]

		# Fix 3  - Sequences are in ascending time
		timefield = self.seqobj.pcodes_time[self.pcode]
		rec_hxb = records[0]
		records= records[1:]
	
		# Get the accessions and time info for this patient		
		cur.execute("SELECT DISTINCT accession,%s FROM %s WHERE \
			patientcode='%s'"%(timefield,self.seqobj.name,self.pcode))
		acc_time = dict(cur.fetchall()) # dict key-acc , val-time
		accs = map(self._get_accname,records)
		times = map(lambda acc : dict.__getitem__(acc_time,acc),accs)
		sorted_tups = sorted(zip(records,times),key=itemgetter(1))
		for rec,time in sorted_tups : 
			rec.name = str(time) +"."+ rec.name
			rec.id = str(time) + "." + rec.id
		records = map(itemgetter(0),sorted_tups)
		records.insert(0,rec_hxb)
		
	
		SeqIO.write(records,self.fpath,format="fasta")

	def _get_accname(self,rec) : 
		""" Gets accession from seq record"""
		return rec.name.split('.')[-1]

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

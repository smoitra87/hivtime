"""
	visualizes time series patient sequences
"""

import sys,os
from operator import itemgetter
import copy,glob
import time
import itertools
from pdb import set_trace as stop

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
from scipy.stats import mannwhitneyu

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
ent_stat_cutoff = 0.6 # used for speeding up stat calculations
sig_lev = 0.05

# If FETCH_ALN is True, then the alignment is downloaded anyway
FETCH_ALN = False
EXEC_JALVIEW = False
DO_PROCESS =False
DO_HOTSPOT = False
DO_STAT = True
DO_STAT_PYMOL = False

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
			analyses = {}		
			for pcode in self.pcodes_time.keys(): 
 				analyses[pcode] = self._calc_stats(pcode)
			self._gen_report_stats(analyses)
	
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
		analysis['plot_res'] = False
		analysis['plot_pair'] = False


		p24fpath = os.path.join(datadir,pcode+'_p24.fasta')
		if not os.path.exists(p24fpath) : 
			raise IOError('File %s not found'%p24fpath)
		aln = AlignIO.read(open(p24fpath),"fasta")
		aln = aln[1:] # Remove reference row
		ent_aa = np.array(self._entropy_aa(aln))

		# find indices of all columsn greater than cutoff
		idx_cols = np.where(ent_aa > ent_stat_cutoff)[0]
		analysis['idx_cols'] = idx_cols;
	
		#print pcode," : ", idx_cols

		# Calculate columnwise statistics		
		times = [int(rec.id.split('.')[0]) for rec in aln ]
		uniq_times = sorted(set(times))		
		if len(uniq_times) > 1 : 
			scale_factor = len(uniq_times)-1 
		else : 
			scale_factor = 1
		
		#--------------------------------------------------------------	
		# List of order means for each candidate column in patient
		patient_col_mean = []
		patient_col_rank = []


		for col in idx_cols:
			col_aa = aln[:,col]
		 	uniq_aa = set(col_aa) # Unique aa in column
			col_aa_mean = {}
			col_aa_rank = {}
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
				col_aa_rank[aa] = ranks_aa
				#scale ranks to be b/w 0 and 1
				ranks_aa_scaled = \
					map(lambda x: (0.0+x)/scale_factor,ranks_aa)
				aa_order_mean = np.mean(ranks_aa_scaled)
				col_aa_mean[aa] = aa_order_mean
			patient_col_mean.append(col_aa_mean)	
			patient_col_rank.append(col_aa_rank)
			# Run Mann whitney U test
		analysis['col_order_mean'] = patient_col_mean	
		analysis['col_rank'] = patient_col_rank

		results_test = self._test_wilcox(patient_col_rank)
		analysis['result_res']= results_test

		if len(patient_col_mean) > 0  :
			analysis['plot_res'] = True
#			self._plot_patient_order_means(patient_col_mean,pcode,\
#				idx_cols,mode=1)

		#--------------------------------------------------------------
		# Order means per column pair
		patient_col_mean = []
		patient_col_rank = []

		for col1, col2 in itertools.combinations(idx_cols,2) : 
			col1_aa = aln[:,col1]
			col2_aa = aln[:,col2]
			col_aa = [s1+s2 for s1,s2 in zip(col1_aa,col2_aa)]
		 	uniq_aa = set(col_aa) # Unique aa in column
			col_aa_mean = {}
			col_aa_rank = {}
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
				col_aa_rank[aa] = ranks_aa
				#scale ranks to be b/w 0 and 1
				ranks_aa_scaled = \
					map(lambda x: (0.0+x)/scale_factor,ranks_aa)
				aa_order_mean = np.mean(ranks_aa_scaled)
				col_aa_mean[aa] = aa_order_mean
			patient_col_mean.append(col_aa_mean)	
			patient_col_rank.append(col_aa_rank)
		analysis['pair_order_mean'] = patient_col_mean
		analysis['pair_rank'] = patient_col_rank

		results_test = self._test_wilcox(patient_col_rank)
		analysis['result_pair']	 = results_test

	
		if len(patient_col_mean) > 0  :	
			analysis['plot_pair'] = True
#			self._plot_patient_order_means(patient_col_mean,pcode,\
#				idx_cols,mode=2)

		return analysis

	def _test_wilcox(self,pat_col_rank) : 
		""" Runs the Mann-Whitney U / Wilcoxon Rank-sum test
for every pair of aa in each candidate column of the alignment. 
Returns the summary of the results
			"""
		result_all = [] ;
		for col in pat_col_rank : 
			result_col = {}
			for aa1,aa2 in itertools.combinations(col.keys(),2) :
				if col[aa1] == col[aa2] : 
					continue
				pval = mannwhitneyu(col[aa1],col[aa2])[1]
				result_col[str(aa1)+"<"+str(aa2)] = pval	
			result_all.append(result_col)
		return result_all


	def _gen_report_stats(self,analyses) :
		""" Generates reports for the statistics calculations"""

		self._page_prev = self.page # save the prev page
		self.page = markup.page()
		page = self.page
		title ="Patient Sequence Statistics page"
		page.h2(title)
		page.init(title=title)
	

		#-----------------------------------------------------------
		# Generate patientwise time dependent residues/columns
		page.h4("Patientwise Time Dependent Residues ")
		pstr = ["""
The figures below display the time preference for each residue type for select columns on a per patient basis. The columns per patient were selected based an entropy cutoff. Time was ordered and shrunk down to a scale of 0 to 1 with 0 being early stage and 1 being late stage. Normally if a residue type is perfectly conserved or does not show any dependence on time then the residue type is likely to have a mean time value(y-value) of around 0.5. However if the residue is significantly above or below the 0.5 dashed blue line then it shows that the residue has a prefernce for occuring in the early or late stages. NOTE :  That this scheme cannot differentiate between residue types that are in middle stage, since they are likely to have a mean y-value of 0.5 as well. 
		"""]
		page.p(pstr)

		page.p("<b>NOTE:</b> Current Significance Level is \\alpha=0.05")

		#--------------------------------------------------------------
		# Order means per column images
		self._write_divs1(analyses)

		#--------------------------------------------------------------
		# Order means per column pairs images
		page.hr()
		page.br()
		page.h4("Patientwise time dependent residue pairs")

		self._write_divs2(analyses)

		#--------------------------------------------------------------
		# Visualization of major statistically significant res and edges
		page.hr()
		page.br()
		self._viz_sig(analyses)
		page.h4("Statistically significant residues visualized on structure")
		figpath="figs/sig_res.png"
		page.img(width=640,height=480,alt="sig_res",\
				src=figpath)
		
		page.h4("Statistically significant edges visualized on structure")
		figpath="figs/sig_res_pair.png"
		page.img(width=640,height=480,alt="sig_res_pair",\
				src=figpath)


		with open("patient_stats.html","w") as fout : 
			fout.write(page.__str__())
		self.page = self._page_prev

	def _viz_sig(self,analyses) : 
		""" Visualize and add image for stat sig residues and pairs. Note that
the residues across all patients are aggregated and then displayed onto the 
protein structure. Normally a residue occurs only in one patient. But in case
there are multiple patients that have the same significant residue then the 
name of the patient is displayed as well
		"""
		# Get the list of the significant residues
		self.sig_res = self._get_sig_elems(analyses,'res')
		# Get the list of the stat significant pairs
		self.sig_pair = self._get_sig_elems(analyses,'pair')

		# Create a Pymol instance and plot them 
		if DO_STAT_PYMOL : 
			self._draw_sig_pymol()			


	def _draw_sig_pymol(self) : 
		""" Draws the significant residues on the protein structure"""
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
		cmd.color('wheat')
		cmd.color("marine","chain A and i. "+\
			"+".join(map(str,np.array(self.sig_res.keys())+1)))	
		cmd.show("spheres","color marine and name CA")
		cmd.set("sphere_transparency","0.3","color marine")
		cmd.set("sphere_scale","0.7","color marine")
		for resi in np.array(self.sig_res.keys())+1 : 
			cmd.label('3H4E//A/'+str(resi)+'/CA',\
				'"%s%s "%(resn,resi)')
		cmd.set("label_size","-2")
		cmd.set('label_position','(0,2,3)')
		cmd.ray()
		sys.stderr.write("Sleeping for 1 sec before saving image")
		time.sleep(1)
		figpath="figs/sig_res.png"
		cmd.png(figpath,width=800,height=600)
	
		#assert False
		for n,e in enumerate(self.sig_pair.keys()) : 
			t = cmd.dist('dist'+str(n),'3H4E//A/'+str(e[0])+'/CA',\
				'3H4E//A/'+str(e[1])+'/CA')
			cmd.color('red','dist'+str(n))
			cmd.show('spheres','3H4E//A'+'/'+str(e[0])+'/CA')
			cmd.show('spheres','3H4E//A'+'/'+str(e[1])+'/CA')
			cmd.color('red','3H4E//A'+'/'+str(e[0])+'/CA')
			cmd.color('red','3H4E//A'+'/'+str(e[1])+'/CA')
			cmd.label('3H4E//A/'+str(e[0])+'/CA','"%s%s "%(resn,resi)')
			cmd.label('3H4E//A/'+str(e[1])+'/CA','"%s%s "%(resn,resi)')
			cmd.hide('labels','dist'+str(n))

		cmd.set("sphere_transparency","0.3","color red")
		cmd.set("sphere_scale","0.7","color red")
		cmd.set('dash_radius','0.2');
		cmd.set('dash_gap','0.0');
		cmd.ray()
		sys.stderr.write("Sleeping for 1 sec before saving image")
		time.sleep(1)
		figpath="figs/sig_res_pair.png"
		cmd.png(figpath,width=800,height=600)
		

	def _get_sig_elems(self,analyses,type_res_pair) :
		""" Helper function that gets a list of significant residues or pairs
based on the argument that is passed to the function
WARNING : The idx returned by the function is one less than the residue number

		"""
		sig_type = {} # init empty dict
		for pcode in filter(lambda x:analyses[x]['plot_'+type_res_pair],\
			self.pcodes_time.keys()) : 
			idx_cols = analyses[pcode]['idx_cols']	
			result = analyses[pcode]['result_'+type_res_pair]
			result = map(lambda x: \
				filter(lambda y: y[1]<sig_lev,x.items()),result)
			result = filter(lambda x: len(x)>0,enumerate(result))
		
			if type_res_pair  == 'pair' : 
				idx_cols=list(itertools.combinations(np.array(idx_cols),2))
	
			for idx,result_elem in result :
				if len(result_elem) == 0 :
					continue 
				if not sig_type.has_key(idx_cols[idx]) : 
					sig_type[idx_cols[idx]] = {}
				sig_type.get(idx_cols[idx])[pcode] = result_elem

		return sig_type

	def _prep_div_tbl1(self,result) : 
		""" Helper function that prepares string for div """
		ret = []
		for col in result :
			ret.append(\
				"</br>".join([h.replace('<','-')+":"+"%.4f"%p for\
				 h,p in col.items() if p < sig_lev]))
		return ret

	def _write_divs1(self,analyses) : 
		""" Writes out divs and images for columns"""
		page = self.page
		ncell_row = 3 # Number of cells per row

		_pcodes =  map(itemgetter(0),\
			filter(lambda x:x[1]['plot_res'],analyses.items()))
		pcodes = filter(lambda x: x in _pcodes,self.pcodes_time.keys())

		page.div(id="container",\
			style="display:table;border: 1px solid black")
		for i in range(len(pcodes))[::ncell_row] :
			page.div(id="row",style="display:table-row")
			for pcode in pcodes[i:i+ncell_row] : 
				page.div(id="cell",\
					style="display:table-cell;padding:1em;"+\
					"border: 1px solid black")	
				figpath="figs/"+pcode+"_res_order_means.png"
				# Saving image
				page.img(width=400,height=300,alt=pcode+"order_means",\
				src=figpath)
				# Printing p-values
				page.table(style="margin:auto;border:1px solid black")	
				page.tr()
				page.td(["Res"+str(r) for r in \
					np.array(analyses[pcode]['idx_cols'])+1],\
						style="margin:auto;border:1px solid black;"+\
						"font-size:70%;padding:0.5em")
				page.tr.close()
				page.tr()
				page.td(self._prep_div_tbl1(\
					analyses[pcode]['result_res']),\
					style="margin:auto;border:1px solid black;"+\
						"font-size:50%;padding:0.5em")
				page.tr.close()
				page.table.close()

				page.div.close()
			page.div.close()
		page.div.close()

	def _write_divs2(self,analyses) : 
		""" Writes out divs and images for column pairs"""
		page = self.page
		ncell_row = 3 # Number of cells per row

		_pcodes =  map(itemgetter(0),\
			filter(lambda x:x[1]['plot_pair'],analyses.items()))
		pcodes = filter(lambda x: x in _pcodes,self.pcodes_time.keys())

		page.div(id="container",\
			style="display:table;border: 1px solid black")
		for i in range(len(pcodes))[::ncell_row] :
			page.div(id="row",style="display:table-row")
			for pcode in pcodes[i:i+ncell_row] : 
				page.div(id="cell",\
					style="display:table-cell;padding:1em;"+\
					"border: 1px solid black")	
				figpath="figs/"+pcode+"_pair_order_means.png"
				page.img(width=400,height=300,alt=pcode+"order_means",\
				src=figpath)
				# Printing p-values
				page.table(style="margin:auto;border:1px solid black")	
				page.tr()
				page.td(["Res"+str(r) for r in \
					np.array(analyses[pcode]['idx_cols'])+1],\
						style="margin:auto;border:1px solid black;"+\
						"font-size:70%;padding:0.5em")
				page.tr.close()
				page.tr()
				page.td(self._prep_div_tbl1(\
					analyses[pcode]['result_pair']),\
					style="margin:auto;border:1px solid black;"+\
						"font-size:50%;padding:0.5em")
				page.tr.close()
				page.table.close()
				page.div.close()
			page.div.close()
		page.div.close()

	def _plot_patient_order_means(self,patient_col_mean,pcode,idx_cols,
		mode):
		""" Make a bar plot like figure of order means for every
		column in a  patient"""
		pl.figure()
		d=0.2
		pl.axis([0-d,len(patient_col_mean)+d,1+d,0-d])
		for i in range(len(patient_col_mean))[1:] :
			pl.axvline(x=i,color='black')
		pl.title("Order Means; Pat - %s"%(pcode))
		pl.xlabel("Residue Number")
		pl.ylabel("Order Scale")
		
		# Choose different plots for cols and pairs
		if mode == 1 : 
			_xticks = (np.arange(len(idx_cols))+0.5,\
				map(str,np.array(idx_cols)+1))
			figname= pcode+"_res_order_means.png"
		if mode == 2 :
			pairs = list(itertools.combinations(np.array(idx_cols)+1,2))
			_xticks = (np.arange(len(pairs))+0.5,\
				[str(c1)+"-"+str(c2) for c1,c2 in pairs])
			figname= pcode+"_pair_order_means.png"
		
		pl.xticks(*_xticks)
		pl.axhline(y=0.5,linestyle="dashed",linewidth=2)		
		
		for ii,aa_col_mean in enumerate(patient_col_mean)  :
			sorted_tups = sorted(aa_col_mean.items(),key=itemgetter(1))
			sorted_aa = map(itemgetter(0),sorted_tups)
			for jj,aa in enumerate(sorted_aa) :
				d = 0.1
				xpos = ii+(jj+1)*1.0/(len(aa_col_mean)+1)
				ypos = aa_col_mean[aa]
				pl.text(xpos,ypos,aa, \
					bbox=dict(facecolor='red', alpha=0.5))
			
		figpath="figs/"+figname
		pl.savefig(figpath)
		#pl.show()
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
		self.fixname = self.pcode+'_fix.fasta'
		self.fixpath = os.path.join(datadir,self.fixname)
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
		self.aln = AlignIO.read(open(self.fixpath),"fasta")
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
		
		# Filter all records where no time info is None
		tups = filter(lambda x : x[1] is not None,zip(records,times))
		sorted_tups = sorted(tups,key=itemgetter(1))

		for rec,time in sorted_tups : 
			rec.name = str(time) +"."+ rec.name
			rec.id = str(time) + "." + rec.id
		records = map(itemgetter(0),sorted_tups)
		records.insert(0,rec_hxb)
		
		SeqIO.write(records,self.fixpath,format="fasta")

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

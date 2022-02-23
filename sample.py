import json
import numpy as np
import pandas as pd

from decorator_equality import class_equality_attributes

from gene import Gene
from snv import SNV
from amplification import Amplification
from chromosome import Chromosome
from position import Position
from segment import Segment
from mutation_rate import Mutation_rate
import matplotlib.pyplot as plt



@class_equality_attributes
class Sample():

	def __init__(self,config,name,save = False):
		self.name = name
		self.config = config
		path_preprocessed_data = config['folder_preprocessed_data']/(name+'.json')
		if path_preprocessed_data.exists():
			self.load_preprocessed_data(path_preprocessed_data)
		else:
			self.preprocess_data()
			if save:
				self.save()

	def preprocess_data(self,threshold_amplification_no_WGD = 5):
		self.clinical_data = self.config['summary_table'].loc[self.config['summary_table']['samplename'] == self.name].iloc[0].to_dict()
		
		for key in self.clinical_data.keys(): #used since json can not be saved with np.bool_ types
			if isinstance(self.clinical_data[key], np.bool_):
				self.clinical_data[key] = bool(self.clinical_data[key])
		self.segments = self.load_cnas()
		SNVs, self.snv_type = self.load_SNVs() # TODO : what to do if gray list
		self.match_snvs_to_segments(SNVs)

		self.major_cn_mode = self.get_major_cn_mode(self.segments)

		self.threshold_amplification = threshold_amplification_no_WGD if self.major_cn_mode <=1 else 2*threshold_amplification_no_WGD
		self.amplifications = self.get_amplifications(self.segments,self.threshold_amplification)
		self.mutation_rate = self.get_mutation_rate(self.segments)
		self.subclonal_structure = pd.read_csv(self.config["folder_subclonal_structure"]/ (self.name+self.config['suffix_subclonal_structure']), sep='\t')

	def load_cnas(self):
		df = pd.read_csv(self.config['folder_cna'] / (self.name+self.config['suffix_cnas']), sep='\t')
		N_segments = df.shape[0]
		list_segments = []

		for i in range(N_segments):
			list_segments.append(Segment(df.iloc[i],self.config,self.clinical_data))
			if i>0:
				if list_segments[i].start <= list_segments[i-1].end:
					raise Exception('CNAs not ordered in input file')

		return list_segments

	def load_SNVs(self):
		vcf_path = self.config['folder_snv_white'] /(self.name+self.config['suffix_snvs'])
		if vcf_path.exists():
			return self.parse_vcf(vcf_path), 'white'

		vcf_path = self.config['folder_snv_gray'] /(self.name+self.config['suffix_snvs'])
		if vcf_path.exists():	
			return self.parse_vcf(vcf_path),'gray'
		else:
			return [], 'missing'

	def parse_vcf(self,vcf_path):
		df = pd.read_csv(vcf_path, sep='\t', comment = '#', header = None, names = self.config['vcf_column_names'])
		SNVs = []
		for index, row in df.iterrows():
			chromosome = Chromosome(row['CHROM'])
			pos = Position(chromosome,row['POS'],self.config)
			info = row['INFO'].split(';')
			t_ref_count = 0
			for s in info:
				if s.startswith(self.config["vcf_prefix_counts"]["ref_count"]):
					t_ref_count = int(s[len(self.config["vcf_prefix_counts"]["ref_count"]):])
					break
			t_alt_count = 0
			for s in info:
				if s.startswith(self.config["vcf_prefix_counts"]["alt_count"]):
					t_alt_count = int(s[len(self.config["vcf_prefix_counts"]["alt_count"]):])
					break
			SNVs.append(SNV(chromosome,
				pos,
				t_ref_count,
				t_alt_count,
				row['REF'],
				row['ALT']))
		SNVs = sorted(SNVs, key=lambda SNV: SNV.pos)
		return SNVs
		
	def match_snvs_to_segments(self,SNVs):
		i = 0 #segment index 
		for SNV in SNVs:
			while SNV.pos > self.segments[i].end:
				i+=1
				if i == len(self.segments):
					#last SNVs not matched 
					return

			if SNV.pos >= self.segments[i].start:
				self.segments[i].add_SNV(SNV)

	def get_major_cn_mode(self,segments,max_cn = 5):
		# returns the most common major copy number (where all copy numbers higher than max_cn are mapped back to max_cn)
		n_bases_cn = [0 for i in range(max_cn+1)]
		for segment in segments:
			cn = min(max_cn,segment.major_cn)
			n_bases_cn[cn] += segment.get_length()
		return n_bases_cn.index(max(n_bases_cn))

	def get_amplifications(self,segments,threshold_amplification):
		potential_amplifications = dict()
		for chromosome in (list(range(1,23))+['X', 'Y']):
			potential_amplifications[str(chromosome)] = dict()
			for arm in ['p','q']:
				potential_amplifications[str(chromosome)][arm] = Amplification(Chromosome(chromosome),arm,[],threshold_amplification,self.clinical_data,self.config)

		for segment in segments:
			if segment.major_cn >=threshold_amplification:
				if segment.start.arm != segment.end.arm:
					print('Warning, amplified segment that spans both arms', self.name, segment.start, segment.end)
				segment.match_genes()
				potential_amplifications[segment.chromosome.chromosome][segment.start.arm].add_segment(segment)

		amplifications = []
		for chromosome in (list(range(1,23))+['X', 'Y']):
			for arm in ['p','q']:
				if len(potential_amplifications[str(chromosome)][arm].segments) != 0:
					potential_amplifications[str(chromosome)][arm].set_oncogenes()
					amplifications.append(potential_amplifications[str(chromosome)][arm])
		return amplifications

	def get_mutation_rate(self,segments):
		# TODO: maybe compute here the different subclonal snvs?
		n_1_0 = 0
		l_1_0 = 0
		n_1_1 = 0
		l_1_1 = 0
		for segment in segments:
			if segment.major_cn == 1:
				if segment.minor_cn == 0:
					n_1_0 += len(segment.SNVs)
					l_1_0 += segment.get_length()
				else:
					n_1_1 += len(segment.SNVs)
					l_1_1 += segment.get_length()
		return Mutation_rate(n_1_0,l_1_0,n_1_1,l_1_1)

	def plot_CN_profile(self,add_snvs = False, chromosome = None, start_pos = 0, max_bases = None,total_cn=False):
		if chromosome is not None:
			if not isinstance(chromosome, Chromosome):
				chromosome = Chromosome(chromosome)
			segments = []
			for s in self.segments:
				if s.chromosome == chromosome:
					segments.append(s)
		else:
			segments = self.segments
		x = np.zeros(2*len(segments)) 
		major = np.zeros(2*len(segments))
		minor = np.zeros(2*len(segments))
		for i,segment in enumerate(segments):
			if chromosome is None:
				x[2*i] = segment.start.absolute_position
				x[2*i+1] = segment.end.absolute_position
			else:
				x[2*i] = segment.start.position
				x[2*i+1] = segment.end.position
			major[2*i:2*(i+1)] = segment.major_cn
			minor[2*i:2*(i+1)] = segment.minor_cn
		fig, ax =  plt.subplots()
		if total_cn:
			plt.plot(x,major+minor,label="Total")
		else:
			plt.plot(x,major,label="Major")
		plt.plot(x,minor,label="Minor")
		if add_snvs:
			x_snvs = []
			mcn = []
			for s in segments:
				for snv in s.SNVs:
					if (snv.alt_count + snv.ref_count) == 0:
						continue
					if chromosome is None:
						x_snvs.append(snv.pos.absolute_position)
					else:
						x_snvs.append(snv.pos.position)
					vaf = snv.alt_count/(snv.alt_count + snv.ref_count)
					tot_cn = s.major_cn+s.minor_cn
					rho = self.clinical_data['purity']
					mcn.append(vaf*(tot_cn + ((1-rho)/rho)*s.get_ploidy_healthy()))
			plt.plot(x_snvs,mcn,marker='.', color='k', linestyle='None',label='Mutations')

		plt.axhline(y=0,color='k', linestyle='--',  linewidth=0.4)
		plt.axhline(y=1,color='k', linestyle='--',  linewidth=0.4)
		if total_cn:
			plt.axhline(y=2,color='k', linestyle='--',  linewidth=0.4)
		plt.axhline(y=self.threshold_amplification,color='k', linestyle='--',  linewidth=0.4)
		if chromosome is None:
			for i in range(len(self.config["cum_length_chr"])):
				plt.axvline(x=self.config["cum_length_chr"][i],color='k', linestyle='--',  linewidth=0.4)
			ax.set_xticks(self.config["cum_length_chr"])
			ax.set_xticklabels('')
			ax.set_xticks([(self.config["cum_length_chr"][i]+self.config["cum_length_chr"][i+1])/2 for i in range(len(self.config["cum_length_chr"])-1)], minor=True)
			ax.set_xticklabels(list(range(1,23))+['X','Y'], minor=True)
			for tick in ax.xaxis.get_minor_ticks():
				tick.tick1line.set_markersize(0)
				tick.tick2line.set_markersize(0)
				tick.label1.set_horizontalalignment('center')
			ax.set_xlabel('Chromosome')
		else:
			if start_pos != 0 or max_bases is not None:
				if max_bases is None:
					ax.set_xlim(start_pos,1.1*x[-1])
				else:
					ax.set_xlim(start_pos,start_pos+max_bases)
			ax.set_xlabel('Position')
		ax.set_ylabel('Copy number')
		plt.legend()
		plt.show()

	def load_preprocessed_data(self,path_preprocessed_data):
		with open(self.config['folder_preprocessed_data']/ (self.name+'.json'), 'r') as fp:
			data = json.load(fp)
		self.clinical_data = data['clinical_data']
		self.snv_type = data['SNVs_calls_type']
		self.segments = [Segment(segment,self.config,self.clinical_data) for segment in data['segments']]
		self.major_cn_mode = data['major_cn_mode']
		self.threshold_amplification = data['threshold_amplification']
		self.subclonal_structure = pd.DataFrame.from_dict(data['subclonal_structure'])

		self.amplifications =  []
		for amplification in data['amplifications']:
			segments = [Segment(segment,self.config,self.clinical_data) for segment in amplification['segments']]
			self.amplifications.append((Amplification(Chromosome(amplification['chromosome']), amplification['arm'], segments,self.threshold_amplification,self.clinical_data,self.config)))
			self.amplifications[-1].set_oncogenes() #TODO information is there, shouldn't need to search again 
		self.mutation_rate = Mutation_rate(data['mutation_rate']['n_1_0'],
			data['mutation_rate']['l_1_0'],
			data['mutation_rate']['n_1_1'],
			data['mutation_rate']['l_1_1'])

	def to_dict(self):
		dic = dict()
		dic['clinical_data'] = self.clinical_data.copy()
		dic['SNVs_calls_type'] = self.snv_type
		dic['segments'] = [segment.to_dict() for segment in self.segments]
		dic['major_cn_mode'] = self.major_cn_mode
		dic['threshold_amplification'] = self.threshold_amplification
		dic['amplifications'] = [amplification.to_dict() for amplification in self.amplifications]
		dic['mutation_rate'] = self.mutation_rate.to_dict()
		dic['subclonal_structure'] = self.subclonal_structure.to_dict('records')
		return dic

	def save(self):
		if not self.config['folder_preprocessed_data'].exists():
			self.config['folder_preprocessed_data'].mkdir(parents=True, exist_ok=True)
		with open(self.config['folder_preprocessed_data']/ (self.name+'.json'), 'w') as fp:
			json.dump(self.to_dict(), fp, indent=4)

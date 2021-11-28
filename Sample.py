import json
import vcf
import numpy as np
import pandas as pd
import time

from utility_classes import Segment, SNV, class_equality_attributes,Amplification,Chromosome,Mutation_rate
from decorator_equality import class_equality_attributes


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
		file = self.config['folder_snv_green'] /(self.name+self.config['suffix_snvs'])
		if file.exists():
			return self.parse_vcf(vcf.Reader(open(file,'r'))), 'green'

		file = self.config['folder_snv_gray'] /(self.name+self.config['suffix_snvs'])
		if file.exists():	
			return self.parse_vcf(vcf.Reader(open(file,'r'))),'gray'
		else:
			return [], 'missing'

	def parse_vcf(self,vcf):
		SNVs = []
		for record in vcf: 
			SNVs.append(SNV(record.CHROM,
				record.POS,
				record.INFO.get('t_ref_count',float('nan')),
				record.INFO.get('t_alt_count',float('nan')),
				record.REF,
				record.ALT,
				self.config))
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
				potential_amplifications[str(chromosome)][arm] = Amplification(Chromosome(chromosome),arm,[])

		for segment in segments:
			if segment.major_cn >=threshold_amplification:
				if segment.start.arm != segment.end.arm:
					print('Warning, amplified segment that spans both arms', self.name, segment.start, segment.end)
				potential_amplifications[segment.chromosome.chromosome][segment.start.arm].add_segment(segment)

		amplifications = []
		for chromosome in (list(range(1,23))+['X', 'Y']):
			for arm in ['p','q']:
				if len(potential_amplifications[str(chromosome)][arm].segments) != 0:
					amplifications.append(potential_amplifications[str(chromosome)][arm])
		return amplifications

	def get_mutation_rate(self,segments):
		# TODO filter out somehow subcolonal SNVs
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

	def load_preprocessed_data(self,path_preprocessed_data):
		with open(self.config['folder_preprocessed_data']/ (self.name+'.json'), 'r') as fp:
			data = json.load(fp)
		self.clinical_data = data['clinical_data']
		self.snv_type = data['SNVs_calls_type']
		self.segments = [Segment(segment,self.config,self.clinical_data) for segment in data['segments']]
		self.major_cn_mode = data['major_cn_mode']
		self.threshold_amplification = data['threshold_amplification']

		self.amplifications =  []
		for amplification in data['amplifications']:
			segments = [Segment(segment,self.config,self.clinical_data) for segment in amplification['segments']]
			self.amplifications.append((Amplification(Chromosome(amplification['chromosome']), amplification['arm'], segments)))
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
		return dic

	def save(self):
		if not self.config['folder_preprocessed_data'].exists():
			self.config['folder_preprocessed_data'].mkdir(parents=True, exist_ok=True)
		with open(self.config['folder_preprocessed_data']/ (self.name+'.json'), 'w') as fp:
			json.dump(self.to_dict(), fp, indent=4)

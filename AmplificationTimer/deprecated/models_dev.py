import pymc3 as pm
import arviz as az
import numpy as np
from scipy.special import betainc
from scipy.special import beta as beta_function
from scipy.stats import beta,binomtest,binom
import matplotlib.pyplot as plt
import theano.tensor as tt

p_value_threshold = 0.05

# TODO: should actually take bayesian approach for mutation rate to avoid mu_ML = 0, now just trick of mu = epsilon if mu_ML = 0


class Model():
	def __init__(self,amplification,mu,subclonal_structure,subclonality_modelled,cores):
		self.cores = cores
		self.subclonality_modelled = subclonality_modelled 
		self.rho = subclonal_structure["fraction_total_cells"][0]
		if subclonality_modelled:
			self.fSNV = subclonal_structure["n_snvs"]/subclonal_structure["n_snvs"].sum()
			self.fraction_tumor_cells_subclones = subclonal_structure["fraction_cancer_cells"]
			self.mu_clonal_ML = mu.get_ML()*self.fSNV[0]
		else:
			self.mu_clonal_ML = mu.get_ML()
		self.model = pm.Model()

	def get_MCMC_samples_posterior(self,n_samples):
		with self.model:
			trace = pm.sample(n_samples, return_inferencedata=False,cores=self.cores)
		return trace

	def get_approximate_MAP(self):
		return pm.find_MAP(model=self.model)

class Model1(Model):
	#Only looks at mutations that are on all copies
	#Hard assignation of mutation state copy number (i.e. before amplification)
	#ML estinate of mu

	def __init__(self,amplification,mu,subclonal_structure,subclonality_modelled,cores,filter_APOBEC):
		Model.__init__(self, amplification, mu,subclonal_structure,subclonality_modelled,cores)
		self.N = 0
		self.L = 0
		for segment in amplification.segments:
			print("SNVs in segment",len(segment.SNVs))
			M = segment.major_cn
			T = segment.major_cn + segment.minor_cn
			ploidy_nc = segment.get_ploidy_healthy()
			p_all = self.rho*M/(T*self.rho+ploidy_nc*(1-self.rho))
			p_1 = self.rho/(T*self.rho+ploidy_nc*(1-self.rho))
			for SNV in segment.SNVs:
				if (SNV.alt_count + SNV.ref_count == 0): continue
				if filter_APOBEC:
					if binomtest(SNV.alt_count, SNV.alt_count + SNV.ref_count, p_all,alternative = 'less').pvalue > p_value_threshold:
						self.N +=1
				else:
					likelihood_1 = binom.pmf(SNV.alt_count,SNV.alt_count + SNV.ref_count,p_1)
					likelihood_all = binom.pmf(SNV.alt_count,SNV.alt_count + SNV.ref_count,p_all)
					if likelihood_all>likelihood_1:
						self.N +=1
			self.L += segment.get_length()
			print("self.N end", self.N)
		print("N",self.N,"L",self.L,"mu", self.mu_clonal_ML)
		self.model_definition()

	def model_definition(self):
		#mu_alpha, mu_beta = self.mu.get_beta_posterior_parameters()
		#mu = pm.Beta("mu", alpha= mu_alpha, beta = mu_beta)
		with self.model:
			t = pm.Uniform("t", lower = 0, upper = 1)
			N_obs = pm.Binomial("N", n=self.L,p = self.mu_clonal_ML*t,observed=self.N)

	def get_analytical_posterior(self,t):
		return(self.mu_clonal_ML*beta.pdf(t*self.mu_clonal_ML, self.N, self.L-self.N)/(betainc(self.N, self.L-self.N,self.mu_clonal_ML)))



class Model2(Model):
	def __init__(self,amplification,mu,subclonal_structure,subclonality_modelled,cores):
		Model.__init__(self, amplification, mu,subclonal_structure,subclonality_modelled,cores)
		self.n_segments = len(amplification.segments)
		n_snvs = 0
		for segment in amplification.segments:
			n_snvs+= len(segment.SNVs)

		#self.N0 = np.zeros(n_segments)
		self.M_segment = np.zeros(self.n_segments)
		self.M_snv = np.zeros((n_snvs,1))
		self.m_snv = np.zeros((n_snvs,1))
		self.d = np.zeros((n_snvs,1))
		self.D = np.zeros((n_snvs,1))
		self.matrix_segment_to_snv = np.zeros((n_snvs,self.n_segments))
		
		if subclonality_modelled:
			self.q = np.zeros((n_snvs,1+len(subclonal_structure.index)))
			self.fSNV_subclonal_matrix = np.zeros((n_snvs,len(subclonal_structure.index)-1))
			for i in range(1,len(subclonal_structure.index)):
				self.fSNV_subclonal_matrix[:,i-1] = self.fSNV[i]
		else:
			self.q = np.zeros((n_snvs,2))

		c_snvs=0
		for i,segment in enumerate(amplification.segments):
			ploidy_nc = segment.get_ploidy_healthy()
			#self.N0[i] = segment.get_length()-len(segment.SNVs)
			self.M_segment[i] = segment.major_cn
			SNVs_kept = []
			for SNV in segment.SNVs:
				q_clonal_one = self.rho/((segment.major_cn + segment.minor_cn)*self.rho+ploidy_nc*(1-self.rho)) # 1 copy clonal
				q_clonal_all = q_clonal_one * segment.major_cn #all copies clonal
				if (SNV.alt_count + SNV.ref_count == 0): continue
				if (binomtest(SNV.alt_count, SNV.alt_count + SNV.ref_count, q_clonal_all,alternative = 'less').pvalue > p_value_threshold ) or (binomtest(SNV.alt_count, SNV.alt_count + SNV.ref_count, q_clonal_one,alternative = 'greater').pvalue > p_value_threshold ):
					SNVs_kept.append(SNV)
					self.d[c_snvs] = SNV.alt_count
					self.D[c_snvs] = SNV.alt_count + SNV.ref_count
					self.M_snv[c_snvs] = segment.major_cn
					self.m_snv[c_snvs] = segment.minor_cn
					self.matrix_segment_to_snv[c_snvs,i] = 1
					self.q[c_snvs,0] = q_clonal_one # 1 copy clonal
					self.q[c_snvs,1] = q_clonal_all #all copies clonal
					if subclonality_modelled:
						for j in range(1,len(subclonal_structure.index)):
							self.q[c_snvs,j+1] = self.q[c_snvs,0] * self.fraction_tumor_cells_subclones[j]
					c_snvs+=1
				segment.SNVs = SNVs_kept
		print("total snvs:", n_snvs,"snvs kept:",c_snvs)
		#Only keep the first entries since some snv were discarded
		self.M_snv = self.M_snv[:c_snvs,:]
		self.m_snv = self.m_snv[:c_snvs,:]
		self.d = self.d[:c_snvs,:]
		self.D = self.D[:c_snvs,:]
		self.q = self.q[:c_snvs,:]
		self.matrix_segment_to_snv = self.matrix_segment_to_snv[:c_snvs,:]
		if subclonality_modelled:
			self.fSNV_subclonal_matrix = self.fSNV_subclonal_matrix[:c_snvs,:]
		self.model_definition()

	def model_definition(self):
		with self.model:
			t = pm.Uniform("t", lower = 0, upper = 1)
			p_1 = (self.M_snv*(1-t))/(self.M_snv*(1-t) + t) # If SNV called #TODO: should actually add also minor
			p_M = 1-p_1
			if self.subclonality_modelled:
				w = pm.math.concatenate((self.fSNV[0]*p_1,self.fSNV[0]*p_M,self.fSNV_subclonal_matrix), axis=1)
			else:
				w = pm.math.concatenate((p_1,p_M), axis=1)
			dist = pm.Binomial.dist(n=self.D, p=self.q)
			d = pm.Mixture('d', w=w, comp_dists=dist, observed=self.d)

class Model3(Model2):
	def model_definition(self):
		include_u = self.m_snv>1
		with self.model:
			t = pm.Uniform("t", lower = 0, upper = 1)
			u = pm.Uniform("u", lower = 0, upper = 1)
			p_1 = (self.M_snv*(1-t)+(u**include_u)*self.m_snv)/(self.M_snv*(1-t) +(u**include_u)*self.m_snv+ t) # If SNV called
			p_M = 1-p_1
			if self.subclonality_modelled:
				w = pm.math.concatenate((self.fSNV[0]*p_1,self.fSNV[0]*p_M,self.fSNV_subclonal_matrix), axis=1)
			else:
				w = pm.math.concatenate((p_1,p_M), axis=1)
			dist = pm.Binomial.dist(n=self.D, p=self.q)
			d = pm.Mixture('d', w=w, comp_dists=dist, observed=self.d)
			
			#p_0 = 1-mu*t-self.M_segment*mu*(1-t)
			#N0_obs = pm.Binomial("N0", n=self.N0,p = p_0,observed=self.N0)

class Model4(Model2):
	def model_definition(self):
		include_u = self.m_snv>1
		with self.model:
			t = pm.Uniform("t", lower = 0, upper = 1)
			u = pm.Uniform("u", lower = 0, upper = 1,shape = (self.n_segments,1))
			u_snv = pm.math.dot(self.matrix_segment_to_snv, u)
			p_1 = (self.M_snv*(1-t)+(u_snv**include_u)*self.m_snv)/(self.M_snv*(1-t) +(u_snv**include_u)*self.m_snv+ t) # If SNV called
			p_M = 1-p_1
			if self.subclonality_modelled:
				w = pm.math.concatenate((self.fSNV[0]*p_1,self.fSNV[0]*p_M,self.fSNV_subclonal_matrix), axis=1)
			else:
				w = pm.math.concatenate((p_1,p_M), axis=1)
			dist = pm.Binomial.dist(n=self.D, p=self.q)
			d = pm.Mixture('d', w=w, comp_dists=dist, observed=self.d)
			
			#p_0 = 1-mu*t-self.M_segment*mu*(1-t)
			#N0_obs = pm.Binomial("N0", n=self.N0,p = p_0,observed=self.N0)

import pymc3 as pm
import arviz as az
import numpy as np
from scipy.special import betainc
from scipy.special import beta as beta_function
from scipy.stats import beta,binomtest
import matplotlib.pyplot as plt

import theano.tensor as tt

p_value_threshold = 0.05

class Model():
	def __init__(self,amplification,mu,clinical_data):
		self.mu = mu
		self.model = pm.Model()

	def get_MCMC_samples_posterior(self,n_samples):
		with self.model:
			trace = pm.sample(n_samples, return_inferencedata=False,cores=1)
		return trace

	def get_approximate_MAP(self):
		return pm.find_MAP(model=self.model)


class Model1(Model):
	def __init__(self,amplification,mu,clinical_data):
		Model.__init__(self, amplification, mu,clinical_data)
		self.N = 0
		self.L = 0
		for segment in amplification.segments:
			M = segment.major_cn
			T = segment.major_cn + segment.minor_cn
			rho = clinical_data['purity']
			ploidy_nc = segment.get_ploidy_healthy()
			p = rho*M/(T*rho+ploidy_nc*(1-rho))
			for SNV in segment.SNVs:
				if binomtest(SNV.alt_count, SNV.alt_count + SNV.ref_count, p,alternative = 'less').pvalue > p_value_threshold:
					self.N +=1
			self.L += segment.get_length()
		print('L',self.L)
		self.model_definition()

	def model_definition(self):
		mu_alpha, mu_beta = self.mu.get_beta_posterior_parameters()
		with self.model:
			t = pm.Uniform("t", lower = 0, upper = 1)
			mu = pm.Beta("mu", alpha= mu_alpha, beta = mu_beta)
			N_obs = pm.Binomial("N", n=self.L,p = mu*t,observed=self.N)

	def get_analytical_posterior(self,t):
		# using ML estimate of mu 
		mu = self.mu.get_ML()
		return(mu*beta.pdf(t*mu, self.N, self.L-self.N)/(betainc(self.N, self.L-self.N,mu)))

class Model2(Model):
	def __init__(self,amplification,mu,clinical_data):
		Model.__init__(self, amplification, mu,clinical_data)
		n_segments = len(amplification.segments)
		n_snvs = 0
		for segment in amplification.segments:
			n_snvs+= len(segment.SNVs)

		self.N0 = np.zeros(n_segments)
		self.M_segment = np.zeros(n_segments)
		self.M_snv = np.zeros((n_snvs,1))
		self.m_snv = np.zeros((n_snvs,1))
		self.d = np.zeros((n_snvs,1))
		self.D = np.zeros((n_snvs,1))
		self.q_1 = np.zeros((n_snvs,1))
		self.q_M = np.zeros((n_snvs,1))

		rho = clinical_data['purity']
		c_snvs=0
		for i,segment in enumerate(amplification.segments):
			ploidy_nc = segment.get_ploidy_healthy()
			self.N0[i] = segment.get_length()-len(segment.SNVs)
			self.M_segment[i] = segment.major_cn
			for SNV in segment.SNVs:
				self.d[c_snvs] = SNV.alt_count
				self.D[c_snvs] = SNV.alt_count + SNV.ref_count
				self.q_1[c_snvs] = rho/((segment.major_cn + segment.minor_cn)*rho+ploidy_nc*(1-rho))
				self.q_M[c_snvs] = self.q_1[c_snvs] * segment.major_cn
				self.M_snv[c_snvs] = segment.major_cn
				self.m_snv[c_snvs] = segment.minor_cn
				c_snvs+=1
		self.q = np.concatenate(( self.q_1,self.q_M),axis=1)
		self.model_definition()
		print('N0',self.N0)
		print('M_snv',self.M_snv)
		print('d',self.d)
		print('D',self.D)
		print('d/D',self.d/self.D)
		print('q_1',self.q_1)
		print('q_M',self.q_M)
		print('q',self.q)


	def model_definition(self):
		mu_alpha, mu_beta = self.mu.get_beta_posterior_parameters()
		with self.model:
			t = pm.Uniform("t", lower = 0, upper = 1)
			mu = pm.Beta("mu", alpha= mu_alpha, beta = mu_beta)
			# u = pm.Uniform("u", lower = 0, upper = 1)
			
			p_1 = (self.M_snv*(1-t))/(self.M_snv*(1-t) + t) # If SNV called
			# p_1 = (self.M_snv*(1-t)+u*self.m_snv)/(self.M_snv*(1-t) +u*self.m_snv+ t) # If SNV called
			p_M = 1-p_1
			w = pm.math.concatenate((p_1,p_M), axis=1)
			dist = pm.Binomial.dist(n=self.D, p=np.concatenate((self.q_1,self.q_M),axis=1))

			d = pm.Mixture('d', w=w, comp_dists=dist, observed=self.d)

			p_0 = 1-mu*t-self.M_segment*mu*(1-t)
			N0_obs = pm.Binomial("N0", n=self.N0,p = p_0,observed=self.N0)

	def get_analytical_posterior(self,t):
		# using ML estimate of mu 
		pass


class Model3(Model2):
	def model_definition(self):
		mu_alpha, mu_beta = self.mu.get_beta_posterior_parameters()
		with self.model:
			t = pm.Uniform("t", lower = 0, upper = 1)
			mu = pm.Beta("mu", alpha= mu_alpha, beta = mu_beta)
			u = pm.Uniform("u", lower = 0, upper = 1)
			
			p_1 = (self.M_snv*(1-t)+u*self.m_snv)/(self.M_snv*(1-t) +u*self.m_snv+ t) # If SNV called
			p_M = 1-p_1
			w = pm.math.concatenate((p_1,p_M), axis=1)
			dist = pm.Binomial.dist(n=self.D, p=np.concatenate((self.q_1,self.q_M),axis=1))

			d = pm.Mixture('d', w=w, comp_dists=dist, observed=self.d)

			p_0 = 1-mu*t-self.M_segment*mu*(1-t)
			N0_obs = pm.Binomial("N0", n=self.N0,p = p_0,observed=self.N0)



import numpy as np
import matplotlib.pyplot as plt
mu = 1; 
class NCD:

	def __init__(self,delta_t,delta_c,t,N):
		self.delta_t = delta_t
		self.delta_c = delta_c
		self.t = t
		self.N = N
	def rDelta_t(self):
		return self.delta_t
	def rDelta_c(self):
		return self.delta_c	
	def time(self):
		return self.t
	def stage(self):
		return self.N	
	def w(self, delta):
		return 1 + delta[1]*np.cos(delta[5]) + delta[2]*np.cos(delta[5])
	def s2(self,delta):
		return 1 + delta[3]**2 + delta[4]**2	
	def rho (self, delta):
		np.array([0,0,0,0,0,np.sqrt(mu*delta[0])*((w(delta)/delta[0])**2)])
	def A(delta):
		return np.sqrt(delta[0]/mu)*np.array([0, 2*delta[0]/w(delta),0],[np.sin(delta[5]),(1/w(delta))*((w(delta)+1)*np.cos(delta[5])+delta[1]), -(delta[2]/w(delta))*(delta[3]*np.sin(delta[5])-delta[4]*np.cos(delta[5]))],[-np.cos(delta[5],(1/w(delta))*((w(delta)+1)*np.sin(delta[5])+delta[2]), (delta[3]/w(delta))*(delta[3]*np.sin(delta[5])-delta[4]*np.cos(delta[5])))],[0,0, (s2(delta)/(2*w(delta))*np.cos(delta[5]))],[0,0,(s2(delta)/(2*w(delta))*np.sin(delta[5])) ],[0,0, (1/w(delta))*(delta[3]*np.sin(delta[5])-delta[4]*np.cos(delta[5]))])
	def dxi(self, delta_t, delta_c, U, Lambda_d):
		rho(delta_t)-rho(delta_c)
		
		return np.
				
	def trajectory_prop(self):
		dt = 0.1
		xi = delta_t-delta_c
	
	def control(self)
	def votage_calculation(self)
	def plot_data(self)


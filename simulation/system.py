import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from poliastro.core import elements as elem
mu = 1; 
r = [-6045, -3490, 2500] * u.km
v = [-3.457, 6.618, 2.533] * u.km / u.s
class NCD:

	def __init__(self,delta_t,delta_c,t,N,U,Lambda_d):
		self.delta_t = delta_t
		self.delta_c = delta_c
		self.t = t
		self.N = N
		self.U = U
		self.Lambda_d = Lambda_d

	def w(self, delta):
		return (1 + delta[1]*np.cos(delta[5]) + delta[2]*np.cos(delta[5]))
	def s2(self,delta):
		return 1 + delta[3]**2 + delta[4]**2	
	def rho(self, delta):
		return np.array([0,0,0,0,0,np.sqrt(mu*delta[0])*((self.w(delta)/delta[0])**2)])
	def A(self,delta):
		return np.sqrt(delta[0]/mu)*np.array([0, 2*delta[0]/self.w(delta),0],[np.sin(delta[5]),(1/self.w(delta))*((self.w(delta)+1)*np.cos(delta[5])+delta[1]), -(delta[2]/self.w(delta))*(delta[3]*np.sin(delta[5])-delta[4]*np.cos(delta[5]))],[-np.cos(delta[5],(1/self.w(delta))*((self.w(delta)+1)*np.sin(delta[5])+delta[2]), (delta[3]/self.w(delta))*(delta[3]*np.sin(delta[5])-delta[4]*np.cos(delta[5])))],[0,0, (self.s2(delta)/(2*self.w(delta))*np.cos(delta[5]))],[0,0,(self.s2(delta)/(2*self.w(delta))*np.sin(delta[5])) ],[0,0, (1/self.w(delta))*(delta[3]*np.sin(delta[5])-delta[4]*np.cos(delta[5]))])
	def dxi(self):
		rho_s = self.rho(self.delta_t)-self.rho(self.delta_c)
		return rho_s + np.dot(self.A(self.delta_t),self.U) + self.Lambda_d
	def trajectory_prop(self):
		dt = 0.1
		xi = self.delta_t-self.delta_c
	def control(self):
		pass
	def votage_calculation(self):
		pass

	def plot_data(self):
		pass


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from astropy import units as u
from poliastro.core import elements as elem
mu = 1; 
r = [-6045, -3490, 2500] * u.km
v = [-3.457, 6.618, 2.533] * u.km / u.s
Xi=[]
Time=[]
class NCD:

	def __init__(self,delta_t,delta_c,t,T,N,U,Lambda_d):
		'''initialize the class'''
		self.delta_t = delta_t
		self.delta_c = delta_c
		self.T = T
		self.t = t
		self.N = N
		self.U = U
		self.Lambda_d = Lambda_d
	def w(self, delta):
		'''w matrix'''
		return (1 + delta[1]*np.cos(delta[5]) + delta[2]*np.cos(delta[5]))
	def s2(self,delta):
		'''s2 matrix'''
		return 1 + delta[3]**2 + delta[4]**2	
	def rho(self, delta):
		'''rho matrix'''
		return np.array([0,0,0,0,0,np.sqrt(mu*delta[0])*((self.w(delta)/delta[0])**2)])
	def A(self,delta):
		'''A matrix'''
		return np.sqrt(delta[0]/mu)*np.array([0, 2*delta[0]/self.w(delta),0],[np.sin(delta[5]),(1/self.w(delta))*((self.w(delta)+1)*np.cos(delta[5])+delta[1]), -(delta[2]/self.w(delta))*(delta[3]*np.sin(delta[5])-delta[4]*np.cos(delta[5]))],[-np.cos(delta[5],(1/self.w(delta))*((self.w(delta)+1)*np.sin(delta[5])+delta[2]), (delta[3]/self.w(delta))*(delta[3]*np.sin(delta[5])-delta[4]*np.cos(delta[5])))],[0,0, (self.s2(delta)/(2*self.w(delta))*np.cos(delta[5]))],[0,0,(self.s2(delta)/(2*self.w(delta))*np.sin(delta[5])) ],[0,0, (1/self.w(delta))*(delta[3]*np.sin(delta[5])-delta[4]*np.cos(delta[5]))])
	def votage_calculation(self):
		'''calculate the votage'''
		pass
	def control(self,xi):
		'''control law'''
		self.U = xi*np.array([0,0,0])
		self.votage_calculation()
	def target(self,delta):
		'''target dynamics'''
		return np.array([0,0,0,0,0,np.sqrt(mu*delta[0])*((self.w(delta)/delta[0])**2)])  + self.Lambda_d*np.array([0,0,0,0,0,0])
	def chaser(self,delta,xi):
		'''chaser dynamics'''
		self.control(xi)
		return self.rho(delta) + np.dot(self.A(delta),self.U) + self.Lambda_d*np.array([0,0,0,0,0,0])	
	def trajectory_prop(self):
		'''propagate the trajectory of the chaser and target'''
		dt = 0.1
		xi = self.delta_t-self.delta_c
		Xi.append(xi)
		while (xi>0):
			self.delta_t = self.delta_t + self.target(self.delta_t)*dt	
			self.delta_c = self.delta_c + self.chaser(self.delta_c,xi)*dt
			self.t=self.t+dt
			xi = self.delta_t-self.delta_c
			Xi.append(xi)
			Time.append(self.t)
		
		def plot_data(self):
			'''plot XI data vs time'''
			plt.plot(Time,Xi)

'''global solution_t,solution_c
		solution_t =solve_ivp(self.target, [0,self.T], self.delta_t, t_eval = np.arange(0,self.T,dt))
		solution_c =solve_ivp(self.chaser, [0,self.T], self.delta_c, t_eval = np.arange(0,self.T,dt))'''			


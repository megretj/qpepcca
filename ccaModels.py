import numpy as np
import matplotlib.pyplot as plt
import math
import scipy

class CCA_CUBIC_MarkovChain_OG:
    def __init__(self, N=400, C = 1000, trans_err = 0.01, alpha = 1, beta = 0.5, RTT = 0.025, bit_err = 0.05):
        self.N = N # number of states
        self.W = C*RTT # in Mbyte, W = C*RTT where C is the maximum bandwidth
        self.RTT = RTT # in ms 
        self.trans_err = trans_err # a.k.a. lambda
        self.bit_err = bit_err # probability that a packet drops 
        self.alpha = alpha # window growth rate
        self.beta = beta # window reduction ratio
        self.a = (np.linspace(1,N,N)-0.5)*self.W/N # Discretisation of the window-size
        self.pi = np.ones(N)/N # Stationnary Distribution
        self.P = None # Transition Probability Matrix
        self.S = None
        self.tau = None 
        self.x = 0 # Steady State average throughput

    def D(self,x,y): # growth time
        return np.cbrt((y-x)/self.alpha)+np.cbrt((1-self.beta)*x/self.alpha)
    
    def avg_w(self,t,x):
        return 1/(t)*(x*t+1/4*(t-np.cbrt((1-self.beta)*x/self.alpha))) #not finished

    def transition_proba_CUBIC(self,i,j): # Probability to transition from i to j
        l = self.trans_err
        if j < self.beta*(i-0.5):
            return 0
        if j == self.N:
            sum = 0
            for jp in np.linspace(1,self.N-1,self.N-1):
                if jp >= self.beta*(i-0.5):
                    sum += self.transition_proba_CUBIC(i,jp)
            return 1-sum
        return np.exp(-l* np.maximum(self.D(self.a[i-1],(j-1)*self.W/self.N),0)) - np.exp(-l* self.D(self.a[i-1],j*self.W/self.N))
    
    def compute_stationnary_distribution(self):
        # 1. Compute the transition probability Matrix P
        self.P = np.zeros([self.N,self.N])
        for i in range(self.N):
            for j in range(self.N):
                if j == self.N:
                    self.P[i,j] = 1-np.sum(self.P[i,:-1])
                    continue
                self.P[i,j] = self.transition_proba_CUBIC(i+1,j+1)

        # 2. Solve the system of equation (16)&(17) 
        # piP = pi <=> pi(P-I)=0 <=> (P-I)^T pi = 0 
        # so pi is a left eigenvector of P, with eigenvalue 1
        # Furthermore, pis L1 norm needs to be equal to 1 
        w,v = np.linalg.eig(np.transpose(self.P)) # Compute eigenvalues/eigenvectors
        self.pi = np.real(v[:,0]/v[:,0].sum()) # Scale such that the values sum to 1
        return 
    
    def compute_tau_and_S(self):
        self.tau = np.zeros([self.N,self.N])
        self.S = np.zeros([self.N,self.N])
        for i in range(self.N):
            for j in range(self.N):
                self.tau[i,j] = np.maximum(self.D(self.a[i],(j-0.5)*self.W/self.N),0) # average time duration for state transistion from state i to j
                L = np.cbrt((1-self.beta)*self.a[i-1]/self.alpha)
                self.S[i,j] = self.a[i]*self.tau[i,j] + self.alpha/4*((self.tau[i,j]-L)**4-L**4)
        return
    
    def avg_throughput(self):
        if self.P == None:
            self.compute_stationnary_distribution()
        self.compute_tau_and_S()
        numerator = np.dot(self.pi,np.dot(self.P,np.transpose(self.S)).diagonal())
        denominator = np.dot(self.pi,np.dot(self.P,np.transpose(self.tau)).diagonal())
        self.x = 1/self.W*numerator/denominator
        return self.x
    
class CCA_CUBIC_MarkovChain_new:
    def __init__(self, N=400, C = 1000, alpha = 1, beta = 0.5, RTT_est = 0.025, RTT_real = 0.025, bit_err = 0.05):
        self.N = N # number of states
        self.W = C*RTT_real # in Mbyte, W = C*RTT where C is the maximum bandwidth
        self.C = C
        self.RTT = RTT_est # in ms 
        self.epsilon = bit_err # probability that a packet drops (epsilon)
        self.alpha = alpha # window growth rate
        self.beta = beta # window reduction ratio
        self.a = (np.linspace(1,N,N)-0.5)*self.W/N # Discretisation of the window-size
        self.pi = np.ones(N)/N # Stationnary Distribution
        self.P = None # Transition Probability Matrix
        self.S = None
        self.tau = None 
        self.x = 0 # Steady State average throughput

    def w(self,x,t):
        return self.alpha*(t-np.cbrt(x*(1-self.beta)/self.alpha))+x

    def T(self,x,b): # growth time
        return np.cbrt((b-x)/self.alpha)+np.cbrt(x*(1-self.beta)/self.alpha)

    def D(self,x,y): # upper bound on the number of segments sent within the growth time
        T = self.T(x,y)/self.RTT
        l = np.floor(T)
        delta = T - l
        sum = 0
        for k in range (int(l)):
            sum += self.w(x,k*self.RTT)
        return sum+delta*self.w(x,l)

    def transition_proba_CUBIC(self,i,j): # Probability to transition from i to j
        if j < self.beta*(i-0.5):
            return 0
        if j == self.N:
            sum = 0
            for jp in np.linspace(1,self.N-1,self.N-1):
                if jp >= self.beta*(i-0.5):
                    sum += self.transition_proba_CUBIC(i,jp)
            return 1-sum
        # TODO: could slightly optimise by not computing D(i,j) twice (when j' = j+1)
        nij_min = int(self.D(self.a[i-1],(j-1)*self.W/self.N))
        nij_max = int(self.D(self.a[i-1],j*self.W/self.N))
        return (1-self.epsilon)**(nij_min-1)*(1-self.epsilon)**(nij_max)
    
    def compute_stationnary_distribution(self):
        # 1. Compute the transition probability Matrix P
        # TODO: This implementation seems to introduce some mistakes, sums up to something higher than zero
        self.P = np.zeros([self.N,self.N])
        for i in range(self.N):
            for j in range(self.N):
                if j == self.N:
                    self.P[i,j] = 1-np.sum(self.P[i,:-1])
                    continue
                self.P[i,j] = self.transition_proba_CUBIC(i+1,j+1)

        # 2. Solve the system of equation (16)&(17) 
        # piP = pi <=> pi(P-I)=0 <=> (P-I)^T pi = 0 
        # so pi is a left eigenvector of P, with eigenvalue 1
        # Furthermore, pis L1 norm needs to be equal to 1 
        # w,v = np.linalg.eig(np.transpose(self.P)) # Compute eigenvalues/eigenvectors

        #TODO: Pick the correct eigenvector, the that is equal to one (not the one with index zero by default)
        # self.pi = np.real(v[:,0]/v[:,0].sum()) # Scale such that the values sum to 1

        ws,vs = scipy.sparse.linalg.eigs(A=np.transpose(self.P),k=1,sigma=1)
        self.pi = np.real(vs/vs.sum())
        print(ws)
        print(self.P)
        return 
    
    def compute_tau_and_S(self):
        self.tau = np.zeros([self.N,self.N])
        self.S = np.zeros([self.N,self.N])
        for i in range(self.N):
            for j in range(self.N):
                self.tau[i,j] = np.maximum(self.T(self.a[i],(j-0.5)*self.W/self.N),0) # average time duration for state transistion from state i to j
                L = np.cbrt((1-self.beta)*self.a[i-1]/self.alpha)
                self.S[i,j] = self.a[i]*self.tau[i,j] + self.alpha/4*((self.tau[i,j]-L)**4-L**4)
        return
    
    def avg_throughput(self):
        # TODO: re-do the derivation from the paper and see if you need S, tau, and how should you compute them
        if self.P == None:
            self.compute_stationnary_distribution()
        self.compute_tau_and_S()
        numerator = np.dot(self.pi,np.dot(self.P,np.transpose(self.S)).diagonal())
        denominator = np.dot(self.pi,np.dot(self.P,np.transpose(self.tau)).diagonal())
        self.x = 1/self.W*numerator/denominator
        return self.x
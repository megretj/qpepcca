

'''
This file provides functions for the models to study congestion control algorithms as presented in the Thesis.
'''


import numpy as np
import matplotlib.pyplot as plt
import math
import scipy
import logging, sys

class CCA_MarkovChain:
    def __init__(self, N:int = 100, C:float=1000., RTT_real:float=0.025, RTT_est:float=0.025, packet_err:float=0.01, err_rate:float = 1, beta:float = 0.7, alpha:float = 1.0):
        self.N = N # number of states
        self.C = C # Bottleneck bandwidth on the path (MSS/s)
        self.W = C*RTT_real # Maximum window size to avoid congestion. In MSS
        self.RTT_real = RTT_real # in s
        self.RTT_est = RTT_est # in s
        self.packet_err = packet_err # probability that a packet drops (for the packet model)
        self.err_rate = err_rate # error rate (for the time model)
        self.beta = beta # window reduction ratio
        self.a = (np.arange(N)+0.5)*self.W/N # Discretisation of the window-size
        self.pi = np.ones(N)/N # Stationnary Distribution
        self.P = np.zeros([self.N,self.N]) # Transition Probability Matrix
        self.S = np.zeros([self.N,self.N]) 
        self.tau = np.zeros([self.N,self.N]) 
        self.alpha = alpha # window growth rate
        self.ssThroughput = 0 # Steady State average throughput

    def update(self):
        self.W = self.C*self.RTT_real # Maximum window size to avoid congestion. In MSS
        self.a = (np.arange(self.N)+0.5)*self.W/self.N
        self.pi = np.ones(self.N)/self.N # Stationnary Distribution
        self.P = np.zeros([self.N,self.N]) # Transition Probability Matrix
        self.S = np.zeros([self.N,self.N]) 
        self.tau = np.zeros([self.N,self.N]) 
        self.ssThroughput = 0 # Steady State average throughput

    def compute_stationnary_distribution(self):
        # 1. Compute the transition probability Matrix P
        for i in range(self.N):
            for j in range(self.N):
                if j == self.N-1:
                    self.P[i,j] = 1-np.sum(self.P[i,:-1])
                else:
                    self.P[i,j] = self.transition_proba(i,j)
        # 2. Solve the system of equation (16)&(17) 
        # piP = pi <=> pi(P-I)=0 <=> (P-I)^T pi = 0 
        # so pi is a left eigenvector of P, with eigenvalue 1
        # Furthermore, pis L1 norm needs to be equal to 1 
        w,v = np.linalg.eig(np.transpose(self.P)) # Compute eigenvalues/eigenvectors
        self.pi = np.real(v[:,0]/v[:,0].sum()) # Scale such that the values sum to 1
        #ws,vs = scipy.sparse.linalg.eigs(A=np.transpose(self.P),k=1,sigma=1)
        #self.pi = np.real(vs/vs.sum())[:,0]
        return self.pi

###----------------------CUBIC START-----------------------------###

class CCA_MarkovChain_CUBIC(CCA_MarkovChain):
    def __init__(self,*args,**kwargs):
        super(CCA_MarkovChain_CUBIC,self).__init__(*args,**kwargs)

    def __str__(self):
        return "CUBIC"
    
    def T(self,x,y):
        """Time to grow from x to y

        Args:
            x (): initial cwnd (in Mb), before window size reduction of beta
            y (): final cwnd (in Mb)

        Returns:
            _type_: Time, in seconds, it takes to grow from beta*x to y.
        """
        return np.cbrt((y-x)/self.alpha)+np.cbrt((1-self.beta)*x/self.alpha)
    
    def w(self,x,t):
        """ Window size dynamics of CUBIC
        """
        return self.alpha*(t-np.cbrt(x*(1-self.beta)/self.alpha))**3+x
    
    def transition_proba(self,i:int,j:int):
        """Pij, probability to transition from state a_i to a_j

        Args:
            i (int): initial state
            j (int): final state

        Returns:
            float: NEED TO REDEFINE THAT FUNCTION FOR SPECIFIC MODELS
        """
        return 0
    
class CCA_MarkovChain_CUBIC_time(CCA_MarkovChain_CUBIC):
    
    def transition_proba(self,i,j): # Probability to transition from i to j
        l = self.err_rate
        if j+1 < self.beta*(i+0.5):
            return 0
        if j+1 == self.N:
            sum = 0
            for jp in range(self.N-1):
                if jp >= self.beta*(i+0.5):
                    sum += self.transition_proba(i,jp)
            return 1-sum
        return np.exp(-l* np.maximum(self.T(self.a[i],(j)*self.W/self.N),0)) - np.exp(-l* self.T(self.a[i],(j+1)*self.W/self.N))
    
    def compute_tau_and_S(self):
        for i in range(self.N):
            for j in range(self.N):
                self.tau[i,j] = np.maximum(self.T(self.a[i],self.a[j]),0) # average time duration for state transistion from state i to j
                L = np.cbrt((1-self.beta)*self.a[i]/self.alpha)
                self.S[i,j] = self.a[i]*self.tau[i,j] + self.alpha/4*((self.tau[i,j]-L)**4-L**4)
        return
    
    def avg_throughput(self):
        self.compute_stationnary_distribution()
        self.compute_tau_and_S()
        numerator = np.dot(self.pi,np.dot(self.P,np.transpose(self.S)).diagonal())
        denominator = np.dot(self.pi,np.dot(self.P,np.transpose(self.tau)).diagonal())
        self.ssThroughput = numerator/denominator/self.W
        return self.ssThroughput
    
class CCA_MarkovChain_CUBIC_packet(CCA_MarkovChain_CUBIC):
    def __init__(self, *args, **kwargs):
        super(CCA_MarkovChain_CUBIC_packet,self).__init__(*args,**kwargs)
        # We need additional variables
        self.Dmin = np.zeros([self.N,self.N]) # Minimum distance Matrix
        self.Dmax = np.zeros([self.N,self.N]) # Maximum distance Matrix
        self.Ptilde = np.zeros((self.N,self.N)) # Shifted transition probabilities
        self.compute_distances()
        
    def update(self):
        super(CCA_MarkovChain_CUBIC_packet,self).update()
        self.compute_distances()
    
    def compute_distances(self):
        for i in range(int(np.ceil(self.N/self.beta))):
            for j in range(i,self.N):
                if j == i:
                    self.Dmin[i,j] = 1
                else:
                    self.Dmin[i,j] = self.Dmax[i,j-1]+1
                self.Dmax[i,j] = self.D(self.a[i],(j+1)*self.W/self.N)
        return

    def D(self,x,y) -> int:
        """Number of segments sent within the transition from x (after reduction) to y (upper bounded since we assume cwnd packets are sent every RTT )

        Args:
            x (_type_): start congestion window size
            y (_type_): target congestion window size

        Returns:
            int: number of segments 
        """
        T = self.T(x*self.beta,y)/self.RTT_real
        l = int(np.floor(T))
        delta = T - l
        sum = 0
        for k in range(l):
            sum += self.w(x,k*self.RTT_real)
        return int(np.floor(sum+delta*self.w(x,l*self.RTT_real)))

    def transition_proba_tilde_CUBIC(self,i, j):
        if j<i:
            return 0
        nmin = self.Dmin[i,j]
        nmax = self.Dmax[i,j]
        #print(f"From {self.a[int(self.beta*i)]} to [{(j-1)*self.W/self.N},{j*self.W/self.N}]")
        #print(nmin,nmax)
        if j == self.N -1:
            return (1-self.packet_err)**(nmin-1)

        return (1-self.packet_err)**(nmin-1)-(1-self.packet_err)**(nmax-1)
    
    def compute_stationnary_distribution(self):
        # 1. Compute the transition probability Matrix P
        # First the shifted version
        for i in range(self.N): # Actually would only need to compute up to (N-1)beta
            for j in range(self.N):
                self.Ptilde[i,j] = self.transition_proba_tilde_CUBIC(i,j)
        # Then recover P from Ptilde
        for i in range(self.N):
            self.P[i,:] = self.Ptilde[int(i*self.beta),:]
        # 2. Solve the system of equation (16)&(17)
        ws,vs = scipy.sparse.linalg.eigs(A=np.transpose(self.P),k=1,sigma=1)
        self.pi = np.real(vs/vs.sum())[:,0]
        return 

    def compute_tau_and_S(self):
        for i in range(self.N):
            for j in range(self.N):
                self.tau[i,j] = np.maximum(self.T(self.a[i],self.a[j]),0) # average time duration for state transistion from state i to j
                L = np.cbrt((1-self.beta)*self.a[i]/self.alpha)
                self.S[i,j] = self.a[i]*self.tau[i,j] + self.alpha/4*((self.tau[i,j]-L)**4-L**4)
        return
    
    def avg_throughput(self):
        self.compute_stationnary_distribution()
        self.compute_tau_and_S()
        numerator = np.dot(self.pi,np.dot(self.P,np.transpose(self.S)).diagonal())
        denominator = np.dot(self.pi,np.dot(self.P,np.transpose(self.tau)).diagonal())
        self.ssThroughput = numerator/denominator*1/self.RTT_real
        #self.ssThroughput= numerator/denominator
        return self.ssThroughput

######------------------------- CUBIC END -----------------------------######
######------------------------- HYBLA START ---------------------------######

class CCA_MarkovChain_Hybla(CCA_MarkovChain):
    def __init__(self, RTT0 = 0.025, *args,**kwargs):
        super(CCA_MarkovChain_Hybla,self).__init__(*args,**kwargs)
        self.RTT0 = RTT0 # target RTT
        self.rho = max(self.RTT_real/self.RTT0,1) # RTT ratio

    def update(self):
        super(CCA_MarkovChain_Hybla,self).update()
        self.rho = max(self.RTT_real/self.RTT0,1)

    def __str__(self):
        return "HYBLA"
    
    def T(self,x,y):
        """Growth time

        Args:
            x (): initial cwnd (in Mb), before window size reduction of beta
            y (): final cwnd (in Mb)

        Returns:
            _type_: Time, in seconds, it takes to grow from beta*x to y.
        """
        if y <= x:
            return 0
        return self.RTT0/(self.alpha*self.rho)*(y-self.beta*x)
    
    def w(self,x,t):
        """ Window size dynamics of Hybla
        """
        return self.alpha*self.rho*t+x

    def transition_proba_Hybla(self,i,j):
        return 0
    
    def compute_stationnary_distribution(self):
        # 1. Compute the transition probability Matrix P
        for i in range(self.N):
            for j in range(self.N):
                if j == self.N-1:
                    self.P[i,j] = 1-np.sum(self.P[i,:-1])
                    continue
                self.P[i,j] = self.transition_proba_Hybla(i,j)

        # 2. Solve the system of equation (16)&(17)
        w,v = np.linalg.eig(np.transpose(self.P)) # Compute eigenvalues/eigenvectors
        self.pi = np.real(v[:,0]/v[:,0].sum()) # Scale such that the values sum to 1
        #ws,vs = scipy.sparse.linalg.eigs(A=np.transpose(self.P),k=1,sigma=1)
        #self.pi = np.real(vs/vs.sum())[:,0]
        return 

class CCA_MarkovChain_Hybla_time(CCA_MarkovChain_Hybla):
    
    def transition_proba_Hybla(self,i,j):
        l = self.err_rate
        if j < self.beta*(i-0.5):
            return 0
        if j == self.N-1:
            sum = 0
            for jp in np.linspace(1,self.N-1,self.N-1):
                if jp >= self.beta*(i-0.5):
                    sum += self.transition_proba_Hybla(i,jp)
            return 1-sum
        return np.exp(-l* np.maximum(self.T(self.a[i-1],(j-1)*self.W/self.N),0)) - np.exp(-l* self.T(self.a[i-1],j*self.W/self.N))
        
    def compute_tau_and_S(self):
        for i in range(self.N):
            for j in range(self.N):
                self.tau[i,j] = np.maximum(self.T(self.a[i],(j-0.5)*self.W/self.N),0) # average time duration for state transistion from state i to j
                self.S[i,j] = (self.a[i]*self.tau[i,j]*self.beta + self.rho*self.tau[i,j]**2/(2*self.RTT0))/self.RTT_real
        return
    
    def avg_throughput(self):
        self.compute_stationnary_distribution()
        self.compute_tau_and_S()
        numerator = np.dot(self.pi,np.dot(self.P,np.transpose(self.S)).diagonal())
        denominator = np.dot(self.pi,np.dot(self.P,np.transpose(self.tau)).diagonal())
        self.ssThroughput = 1/self.W*numerator/denominator
        return self.ssThroughput

class CCA_MarkovChain_Hybla_ssd(CCA_MarkovChain_Hybla_time):
    def avg_w(self,t,x): # packets sent every RTT (for state-depedent lambda)
        return self.beta*x+self.rho*t/(2*self.RTT0)
    
    def transition_proba_Hybla(self,i,j):
        if j < self.beta*(i-0.5):
            return 0
        if j == self.N:
            sum = 0
            for jp in np.linspace(1,self.N-1,self.N-1):
                if jp >= self.beta*(i-0.5):
                    sum += self.transition_proba_Hybla(i,jp)
            return 1-sum
        t_min = np.maximum(self.T(self.a[i-1],(j-1)*self.W/self.N),0)
        t_max = self.T(self.a[i-1],j*self.W/self.N)

        return np.exp(-t_min*self.RTT/(self.avg_w(t=t_min,x=self.a[i-1])*self.packet_err))-np.exp(-t_max*self.RTT_real/(self.avg_w(t=t_max,x=self.a[i-1])*self.packet_err))

class CCA_MarkovChain_Hybla_discrete(CCA_MarkovChain_Hybla):
    def __init__(self, *args,**kwargs):
        super(CCA_MarkovChain_Hybla_discrete,self).__init__(*args,**kwargs)
        # self.W = self.C*self.RTT0 # in Mbyte, where C is the maximum bandwidth
        # self.a = (np.linspace(1,self.N,self.N)-0.5)*self.W/self.N # Discretisation of the window-size
        self.Dmin = np.zeros([self.N,self.N])
        self.Dmax = np.zeros([self.N,self.N])
        self.Ptilde = np.zeros((self.N,self.N))
        self.N_avg = np.zeros((self.N,self.N)) #shifted version of D_avg so that i corresponds to the starting value before the window reduction
        self.compute_distances()
    
    def compute_distances(self):
        """Computes the number of packets required to transition from any state a_i to any state a_j.

        """
        for i in range(self.N):
            for j in range(i,self.N):
                if j == i:
                    self.Dmin[i,j] = 1
                else:
                    self.Dmin[i,j] = self.Dmax[i,j-1]+1
                self.Dmax[i,j] = max(self.D(self.a[i],(j+1)*self.W/self.N),1)

        D_avg = (self.Dmin+self.Dmax)/2 
        for i in range(self.N):
            self.N_avg[i,:] = D_avg[max(int((i+1)*self.beta)-1,0),:]
        return

    def transition_proba_tilde_Hybla(self,i, j):
        if j<i or self.Dmin[i,j]>self.Dmax[i,j]:
            return 0
        nmin = self.Dmin[i,j]
        nmax = self.Dmax[i,j]

        if j == self.N -1:
            return (1-self.packet_err)**(nmin-1)
        
        #print(f"From {self.a[int(self.beta*i)]} to [{(j-1)*self.W/self.N},{j*self.W/self.N}]")
        #print(nmin,nmax)
        return (1-self.packet_err)**(nmin-1)-(1-self.packet_err)**(nmax-1)
    
    def compute_stationnary_distribution(self):
        # 1. Compute the transition probability Matrix P
        # First the shifted version
        for i in range(self.N): # Actually would only need to compute up to (N-1)beta
            for j in range(self.N):
                #if j == self.N-1:
                #    self.Ptilde[i,j] = 1-np.sum(self.Ptilde[i,:-1])
                #    continue
                self.Ptilde[i,j] = self.transition_proba_tilde_Hybla(i,j)
        # Then recover P from Ptilde
        for i in range(self.N):
            self.P[i,:] = self.Ptilde[int(i*self.beta),:]
        # 2. Solve the system of equation (16)&(17)
        ws,vs = scipy.sparse.linalg.eigs(A=np.transpose(self.P),k=1,sigma=1)
        self.pi = np.real(vs/vs.sum())[:,0]
        return self.pi

    def D(self,a,b) -> int:
        """ Number of packets sent between a and b

        Args:
            a (_type_): initial window size (after the window size reduction)
            b (_type_): target window size

        Returns:
            int: number of packets sent between a and b
        """
        if b<=a:
            return 0
        U = (b-a)/self.rho**2
        Uf = np.floor(U)
        return int(U*(a+Uf*self.rho**2)-self.rho**2*(Uf+Uf**2)/2)

    def compute_times(self):
        for i in range(self.N):
            for j in range(self.N):
                self.tau[i,j] = np.maximum(self.T(self.a[i],self.a[j]),0) # average time duration for state transistion from state i to j
        return
    
    def avg_throughput(self):
        self.compute_stationnary_distribution()
        self.compute_times()
        num = 0
        denom = 0
        self.ssThroughput = 0
        for i in range(self.N):
            for j in range(self.N):
                if self.tau[i,j]>0:
                    num+= (self.beta * self.a[i]+self.a[j])/2*self.pi[i]*self.P[i,j]*self.tau[i,j]
                    denom += self.pi[i]*self.P[i,j]*self.tau[i,j]
        self.ssThroughput=num/denom/self.RTT_real
        return self.ssThroughput

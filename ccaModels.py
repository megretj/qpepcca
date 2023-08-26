

'''
This file provides functions for the models to study congestion control algorithms as presented in the Thesis.
When changing the parameters, it is safer to create a new instance of the class as many variables depend on each other and are only updated in the intitialisation. 
'''


import numpy as np
import matplotlib.pyplot as plt
import math
import scipy
import logging, sys
import warnings
import ccaModelHelper as helper
warnings.filterwarnings("error")

class CCA_MarkovChain:
    def __init__(self, N:int = 100, C:float=1000., RTT_real:float=0.025, packet_err:float=0.01, err_rate:float = 1, beta:float = 0.7, alpha:float = 1.0):
        self.N = N # number of states
        self.C = C # Bottleneck bandwidth on the path (MTU/s)
        self.W = C*RTT_real # Maximum window size to avoid congestion. In MTU, so again upper bounding the real maximum window size.
        self.RTT_real = RTT_real # in seconds
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
        self.W = self.C*self.RTT_real # Maximum window size to avoid congestion. In MTU
        self.a = (np.arange(self.N)+0.5)*self.W/self.N
        self.pi = np.ones(self.N)/self.N # Stationnary Distribution
        self.P = np.zeros([self.N,self.N]) # Transition Probability Matrix
        self.S = np.zeros([self.N,self.N]) 
        self.tau = np.zeros([self.N,self.N]) 
        self.ssThroughput = 0 # Steady State average throughput

    def compute_transition_matrix(self):
        for i in range(self.N):
            for j in range(self.N):
                if j == self.N-1:
                    self.P[i,j] = 1-np.sum(self.P[i,:-1])
                else:
                    self.P[i,j] = self.transition_proba(i,j)
        return self.P
    
    def compute_stationnary_distribution(self):
        # 1. Compute the transition probability Matrix P
        self.compute_transition_matrix()
        # 2. Find the stationnary distribution
        
        # TODO: Could potentially check for all absorbing states by looking on the diagonal and set the corresponding values in pi to 1. 
        if abs(self.P[0,0]-1)<1e-5: # First state is absorbing
            self.pi = np.array([1]+[0]*(self.N-1))
            return self.pi
        cols = np.arange(self.N)
        zero_cols = []
        reduP = self.P.copy()
        while np.where(~reduP.any(axis=0))[0].size > 0:
            # To solve the issue with unreachable states, we check if a column is zero
            # if so, we remove the corresponding row and column and solve the system on the reduced matrix
            # We need to do this recursively until we have non-zero columns
            temp_zero_cols = np.where(~reduP.any(axis=0))[0] # find zero columns
            zero_cols.extend([cols[i] for i in temp_zero_cols]) # add original indices to the list of zero columns
            reduP = np.delete(reduP,temp_zero_cols,axis=0) # remove zero columns
            reduP = np.delete(reduP,temp_zero_cols,axis=1) # remove zero rows
            cols = np.delete(cols,temp_zero_cols,axis=0) # remove indices that correspond to zero columns
        try:
            # w,v = np.linalg.eig(np.transpose(self.P)) # Compute eigenvalues/eigenvectors
            # self.pi = np.real(v[:,0]/v[:,0].sum()) # Scale such that the values sum to 1
            ws,vs = scipy.sparse.linalg.eigs(A=np.transpose(reduP),k=1,sigma=1)
            redupi = np.real(vs/vs.sum())[:,0]
        except:
            # This makes it possible to continue the simulation and see on the graph that there was an error
            print(f"Error in computing the reduced stationnary distribution.")
            print(f"P had rank {np.linalg.matrix_rank(self.P)}. Reduced P has rank {np.linalg.matrix_rank(reduP)} Zero columns where at indices: {zero_cols}")
            print(f"W= {self.W}, C= {self.C}, RTT={self.RTT_real}, packet_err= {self.packet_err}")
            print(reduP)
            return self.pi
        self.pi = helper.replace_indices_with_zeros(redupi,zero_cols)
        #helper.well_defined_pi(self.pi)

        ### Simpler version that does not take into account the zero columns ###
        # try:
        #     ws,vs = scipy.sparse.linalg.eigs(A=np.transpose(self.P),k=1,sigma=1)
        #     self.pi = np.real(vs/vs.sum())[:,0]
        # except: 
        #     print(f"Could not compute eigenvalues. Setting pi to uniform. Parameters were: W= {self.W}, C= {self.C}, RTT={self.RTT_real}, packet_err= {self.err_rate}")
        #     self.pi = np.ones(self.N)/self.N
        #######################################################################
        return self.pi
    
    def avg_throughput(self):
        """Computes the average throughput of a segment instance.

        Returns:
            _type_: _description_
        """
        self.compute_stationnary_distribution()
        self.compute_tau_and_S()
        numerator = np.dot(self.pi,np.dot(self.P,np.transpose(self.S)).diagonal())
        denominator = np.dot(self.pi,np.dot(self.P,np.transpose(self.tau)).diagonal())
        self.ssThroughput = numerator/denominator/self.W
        return self.ssThroughput
    
    def sscwnd(self):
        if self.ssThroughput == 0:
            self.avg_throughput()
        return self.ssThroughput*self.W

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
    
    def D(self,a,b) -> int:
        """Number of segments sent within the transition from a (after reduction) to b (upper bounded since we assume cwnd packets are sent every RTT )

        Args:
            x (_type_): start congestion window size
            y (_type_): target congestion window size

        Returns:
            int: number of segments 
        """
        if b<=a:
            return 0
        T = self.T(a/self.beta,b)/self.RTT_real
        l = int(np.floor(T))
        delta = T - l
        sum = 0
        for k in range(l):
            sum += self.w(a/self.beta,k*self.RTT_real)
        return sum+delta*self.w(a/self.beta,l*self.RTT_real)
    
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
        self.Dmin = np.zeros([self.N,self.N]) # Minimum distance Matrix
        self.Dmax = np.zeros([self.N,self.N]) # Maximum distance Matrix
        self.Ptilde = np.zeros((self.N,self.N)) # Shifted transition probabilities
        self.compute_distances()
    
    def compute_distances(self):
        for i in range(self.N):
            for j in range(i,self.N):
                if j == i:
                    self.Dmin[i,j] = 0
                    continue
                self.Dmin[i,j] = max(self.D(self.a[i],j*self.W/self.N),self.Dmin[i,j-1])
                self.Dmax[i,j-1] = self.Dmin[i,j]
        return

    def transition_proba_tilde_CUBIC(self,i, j):
        if j<i:
            return 0
        nmin = self.Dmin[i,j]
        nmax = self.Dmax[i,j]
        #print(f"From {self.a[int(self.beta*i)]} to [{(j-1)*self.W/self.N},{j*self.W/self.N}]")
        #print(nmin,nmax)
        # if j == self.N-1:
        #     return (1-self.packet_err)**(nmin-1)

        # return (1-self.packet_err)**(nmin-1)-(1-self.packet_err)**(nmax-1)
        if j == self.N-1:
            return np.exp(-self.packet_err*nmin)
        return np.exp(-self.packet_err*nmin) - np.exp(-self.packet_err*nmax)
    
    def compute_transition_matrix(self):
        # First the shifted version
        for i in range(self.N): # Actually would only need to compute up to (N-1)beta
            for j in range(self.N):
                self.Ptilde[i,j] = self.transition_proba_tilde_CUBIC(i,j)
        # Then recover P from Ptilde
        for i in range(self.N):
            self.P[i,:] = self.Ptilde[int(i*self.beta),:]
        #self.P = self.P.round(3)
        return self.P

    def compute_tau_and_S(self):
        for i in range(self.N):
            for j in range(self.N):
                self.tau[i,j] = np.maximum(self.T(self.a[i],self.a[j]),0) # average time duration for state transistion from state i to j
                L = np.cbrt((1-self.beta)*self.a[i]/self.alpha)
                self.S[i,j] = self.a[i]*self.tau[i,j] + self.alpha/4*((self.tau[i,j]-L)**4-L**4)
        return

class CCA_MarkovChain_CUBIC_bit(CCA_MarkovChain_CUBIC_packet):
    def __init__(self, MTU = 1500, *args, **kwargs):
        self.packet_size = MTU*8 # in bits
        super(CCA_MarkovChain_CUBIC_bit,self).__init__(*args,**kwargs)
        self.bit_err = self.packet_err/(MTU*8) # probability that a bit is corrupted

    def compute_distances(self):
        for i in range(self.N):
            for j in range(i,self.N):
                if j == i:
                    self.Dmin[i,j] = 0
                    continue
                try:
                    self.Dmin[i,j] = max(self.D(self.a[i],j*self.W/self.N),self.Dmin[i,j-1])
                except RuntimeWarning:
                    print(f"Error in computing Dmin = {self.D(self.a[i],j*self.W/self.N)} times {self.packet_size} at i={i},j={j}")
                self.Dmax[i,j-1] = self.Dmin[i,j]
        self.Dmax[:,-1] = self.Dmin[:,-1] # Doesn't matter what we put here, since we never use it
        self.Dmin *= self.packet_size
        self.Dmax *= self.packet_size
        return
    
    def transition_proba_tilde_CUBIC(self,i, j):
        if j<i:
            return 0
        nmin = int(self.Dmin[i,j])
        nmax = int(self.Dmax[i,j])
        if j == self.N-1:
            return (1-self.bit_err)**(nmin-1)
        return (1-self.bit_err)**(nmin-1)-(1-self.bit_err)**(nmax-1)

class CCA_MarkovChain_CUBIC_Poojary():
    def __init__(self, N:int = 100, C:float=1000., RTT_real:float=0.025, packet_err:float=0.01, err_rate:float = 1, beta:float = 0.7, alpha:float = 1.0):
        self.N = N # number of states
        self.C = C # Bottleneck bandwidth on the path (MTU/s)
        self.W = C*RTT_real # Maximum window size to avoid congestion. In MTU
        self.RTT_real = RTT_real # in seconds
        self.packet_err = packet_err # probability that a packet drops (for the packet model)
        self.beta = beta # window reduction ratio
        self.alpha = alpha # window growth rate
        self.ssThroughput = 0 # Steady State average throughput
    
    def __str__(self):
        return "CUBIC Poojary"
    
    def wn_cub(self,x,t,r):
        return self.alpha*(r*t-np.cbrt(x*(1-self.beta)/self.alpha))**3+x

    def simulate(self, T:int = 1000):
        """Simulates the Markov Chain for T steps

        Args:
            T (int, optional): Number of steps to simulate. Defaults to 1000.
        """
        states = np.zeros(1000)
        T = 1000
        initial_state = 2
        epoch_state = initial_state
        epoch_time = 0

        for t in range(T):
            if t == 0:
                states[t] = initial_state
                continue
            prob_no_loss = (1-self.packet_err)**(states[t-1])
            if np.random.choice([0,1],p=[prob_no_loss,1-prob_no_loss])==1:
                # Packet was lost
                states[t] = states[t-1]*self.beta
                epoch_state = states[t-1]
                epoch_time = 0
            elif self.wn_cub(epoch_state,epoch_time,self.RTT_real) > self.W:
                # Congestion
                states[t] = states[t-1]*self.beta
                epoch_state = states[t-1]
                epoch_time = 0
            else:
                # Packet was not lost
                states[t] = self.wn_cub(epoch_state,epoch_time,self.RTT_real)
                epoch_time += 1
        return states
######------------------------- CUBIC END -----------------------------######
######------------------------- HYBLA START ---------------------------######

class CCA_MarkovChain_Hybla(CCA_MarkovChain):
    def __init__(self, RTT0 = 0.05, *args,**kwargs):
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
        return self.RTT_real*(y-self.beta*x)/(self.rho**2)

    def D(self,a,b) -> int:
        """ Number of packets sent between a and b

        Args:
            a (_type_): initial window size (after the window size reduction)
            b (_type_): target window size

        Returns:
            int: number of packets sent between a and b
        """
        return (b*b-a*a)/(2*self.rho*self.rho)

class CCA_MarkovChain_Hybla_time(CCA_MarkovChain_Hybla):
    
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
                self.S[i,j] = (self.a[i]*self.tau[i,j]*self.beta + self.rho*self.tau[i,j]**2/(2*self.RTT0))
        return

class CCA_MarkovChain_Hybla_packet(CCA_MarkovChain_Hybla):
    def __init__(self, *args,**kwargs):
        super(CCA_MarkovChain_Hybla_packet,self).__init__(*args,**kwargs)
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
                    continue
                self.Dmin[i,j] = self.D(self.a[i],j*self.W/self.N)
                self.Dmax[i,j-1] = self.Dmin[i,j]-1

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
        # OLD VERSION
        # if b<=a:
        #     return 0
        # U = (b-a)/self.rho**2
        # Uf = np.floor(U)
        # return int(U*(a+Uf*self.rho**2)-self.rho**2*(Uf+Uf**2)/2)
        return int((b*b-a*a)/(2*self.rho*self.rho))

    def compute_tau_and_S(self):
        for i in range(self.N):
            for j in range(self.N):
                self.tau[i,j] = np.maximum(self.T(self.a[i],self.a[j]),0) # average time duration for state transition from state i to j
                self.S[i,j] = self.tau[i,j]*(self.beta*self.a[i]+self.a[j])/2
        return
    
    def avg_throughput(self):
        self.compute_stationnary_distribution()
        self.compute_tau_and_S()
        numerator = np.dot(self.pi,np.dot(self.P,np.transpose(self.S)).diagonal())
        denominator = np.dot(self.pi,np.dot(self.P,np.transpose(self.tau)).diagonal())
        self.ssThroughput = numerator/denominator*1/self.RTT_real
        return self.ssThroughput

class CCA_MarkovChain_Hybla_packet_new(CCA_MarkovChain_Hybla):
    def __init__(self, *args,**kwargs):
        super(CCA_MarkovChain_Hybla_packet_new,self).__init__(*args,**kwargs)
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
                    self.Dmin[i,j] = 0
                    continue
                self.Dmin[i,j] = max(self.D(self.a[i],j*self.W/self.N),1)
                self.Dmax[i,j-1] = self.Dmin[i,j]
        return

    def transition_proba_tilde_Hybla(self,i, j):
        if j<i:
            return 0
        nmin = self.Dmin[i,j]
        nmax = self.Dmax[i,j]

        # if j == self.N -1:
        #     return (1-self.packet_err)**(nmin-1)

        # return (1-self.packet_err)**(nmin-1)-(1-self.packet_err)**(nmax-1)
        if j == self.N-1:
            return np.exp(-self.packet_err*nmin)
        return np.exp(-self.packet_err*nmin) - np.exp(-self.packet_err*nmax)
    
    def compute_transition_matrix(self):
        for i in range(self.N): # Actually would only need to compute up to (N-1)beta
            for j in range(self.N):
                self.Ptilde[i,j] = self.transition_proba_tilde_Hybla(i,j)
        # Then recover P from Ptilde
        for i in range(self.N):
            self.P[i,:] = self.Ptilde[int(i*self.beta),:]
        return self.P

    def compute_tau_and_S(self):
        for i in range(self.N):
            for j in range(self.N):
                self.tau[i,j] = np.maximum(self.T(self.a[i],self.a[j]),0) # average time duration for state transition from state i to j
                self.S[i,j] = self.tau[i,j]*(self.beta*self.a[i]+self.a[j])/2
        return
    
class CCA_MarkovChain_Hybla_bit(CCA_MarkovChain_Hybla_packet_new):
    def __init__(self, MTU = 1500, *args, **kwargs):
        self.packet_size = MTU*8 # in bits
        super(CCA_MarkovChain_Hybla_bit,self).__init__(*args,**kwargs)
        self.bit_err = self.packet_err/(MTU*8) # probability that a bit is corrupted

    def compute_distances(self):
        """Computes the number of bits required to transition from any state a_i to any state a_j.
        """
        for i in range(self.N):
            for j in range(i,self.N):
                if j == i:
                    self.Dmin[i,j] = 0
                    continue
                try:
                    self.Dmin[i,j] = max(self.D(self.a[i],j*self.W/self.N),self.Dmin[i,j-1])
                except RuntimeWarning:
                    print(f"Error in computing Dmin = {self.D(self.a[i],j*self.W/self.N)} times {self.packet_size} at i={i},j={j}")
                self.Dmax[i,j-1] = self.Dmin[i,j]
        self.Dmax[:,-1] = self.Dmin[:,-1] # Doesn't matter what we put here, since we never use it
        self.Dmin *= self.packet_size
        self.Dmax *= self.packet_size
        return
    
    def transition_proba_tilde_Hybla(self,i, j):
        if j<i:
            return 0
        nmin = int(self.Dmin[i,j])
        nmax = int(self.Dmax[i,j])
        if j == self.N-1:
            return (1-self.bit_err)**(nmin-1)
        return (1-self.bit_err)**(nmin-1)-(1-self.bit_err)**(nmax-1)
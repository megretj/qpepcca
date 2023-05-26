import numpy as np
import matplotlib.pyplot as plt
import math
import scipy

class CCA_MarkovChain:
    def __init__(self, N:int = 400, C:float=1000., RTT_real:float=0.025, RTT_est:float=0.025, packet_err:float=0.01, err_rate:float = 1):
        self.N = N # number of states
        self.W = C*RTT_real # in Mbyte, W = C*RTT where C is the maximum bandwidth
        self.RTT_real = RTT_real
        self.RTT_est = RTT_est # in ms
        self.packet_err = packet_err # probability that a packet drops 
        self.err_rate = err_rate # error rate
        self.a = (np.linspace(1,N,N)-0.5)*self.W/N # Discretisation of the window-size
        self.pi = np.ones(N)/N # Stationnary Distribution
        self.P = None # Transition Probability Matrix
        self.S = None
        self.tau = None 
        self.ssThroughput = 0 # Steady State average throughput

class CCA_MarkovChain_CUBIC(CCA_MarkovChain):
    def __init__(self,alpha = 1, beta = 0.5,*args,**kwargs):
        super(CCA_MarkovChain_CUBIC,self).__init__(*args,**kwargs)
        self.alpha = alpha # window growth rate
        self.beta = beta # window reduction ratio
    
    def T(self,x,y):
        """Growth time

        Args:
            x (): initial cwnd (in Mb), before window size reduction of beta
            y (): final cwnd (in Mb)

        Returns:
            _type_: Time, in seconds, it takes to grow from beta*x to y.
        """
        return np.cbrt((y-x)/self.alpha)+np.cbrt((1-self.beta)*x/self.alpha)
    
    def w(self,x,t):
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
    
    def compute_stationnary_distribution(self):
        # 1. Compute the transition probability Matrix P
        self.P = np.zeros([self.N,self.N])
        for i in range(self.N):
            for j in range(self.N):
                if j == self.N:
                    self.P[i,j] = 1-np.sum(self.P[i,:-1])
                    continue
                self.P[i,j] = self.transition_proba(i+1,j+1)

        # 2. Solve the system of equation (16)&(17) 
        # piP = pi <=> pi(P-I)=0 <=> (P-I)^T pi = 0 
        # so pi is a left eigenvector of P, with eigenvalue 1
        # Furthermore, pis L1 norm needs to be equal to 1 
        # w,v = np.linalg.eig(np.transpose(self.P)) # Compute eigenvalues/eigenvectors
        # self.pi = np.real(v[:,0]/v[:,0].sum()) # Scale such that the values sum to 1
        ws,vs = scipy.sparse.linalg.eigs(A=np.transpose(self.P),k=1,sigma=1)
        self.pi = np.real(vs/vs.sum())[:,0]
        return self.pi

    
class CCA_MarkovChain_CUBIC_OG(CCA_MarkovChain_CUBIC):
    # Don't need to redefine the init as it is inherited from the CUBIC MarkovChain
    def transition_proba(self,i,j): # Probability to transition from i to j
        l = self.err_rate
        if j < self.beta*(i-0.5):
            return 0
        if j == self.N:
            sum = 0
            for jp in np.linspace(1,self.N-1,self.N-1):
                if jp >= self.beta*(i-0.5):
                    sum += self.transition_proba(i,jp)
            return 1-sum
        return np.exp(-l* np.maximum(self.T(self.a[i-1],(j-1)*self.W/self.N),0)) - np.exp(-l* self.T(self.a[i-1],j*self.W/self.N))
    
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
        if self.P == None:
            self.compute_stationnary_distribution()
        self.compute_tau_and_S()
        numerator = np.dot(self.pi,np.dot(self.P,np.transpose(self.S)).diagonal())
        denominator = np.dot(self.pi,np.dot(self.P,np.transpose(self.tau)).diagonal())
        self.x = 1/self.W*numerator/denominator
        return self.x
    
class CCA_CUBIC_MarkovChain_new(CCA_MarkovChain_CUBIC):
    def __init__(self, distribution = "bernoulli", *args, **kwargs):
        super(CCA_CUBIC_MarkovChain_new,self).__init__(*args,**kwargs)
        self.distribution = distribution

    def D(self,x,y):
        """Number of segments sent within the transition from x to y (upper bounded since we assume cwnd packets are sent every RTT )

        Args:
            x (_type_): _description_
            y (_type_): _description_

        Returns:
            _type_: _description_
        """
        T = self.T(x,y)/self.RTT_real
        #print(f"Time from {x} to {y} is {self.T(x,y)}")
        l = int(np.floor(T))
        delta = T - l
        # print(f"So it has time to send {T} = {l} + {delta} cwnds")
        sum = 0
        for k in range(l):
            sum += self.w(x,k*self.RTT_real)
            # if k%10 == 0:
            #     print(sum)
        return np.floor(sum+delta*self.w(x,l*self.RTT_real))

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
        nij_min = max(int(self.D(self.a[i-1],(j-1)*self.W/self.N)),0)
        nij_max = int(self.D(self.a[i-1],j*self.W/self.N))
        #print(f"a_i = {self.a[i-1]}, a_j = {self.a[j-1]} in ({(j-1)*self.W/self.N},{j*self.W/self.N}), n_min ={nij_min}, n_max={nij_max}")
        if self.distr == "bernoulli":
            return (1-self.epsilon)**(nij_min-1)*(1-self.epsilon)**(nij_max)
        else:
            s = 0
            for k in range(nij_min,nij_max+1):
                s+= (1/self.epsilon)**k/math.factorial(k)

            return np.exp(-1/self.epsilon)*s
    
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
    

class CCA_Hybla_MarkovChain:
    def __init__(self, N=400, C = 100, trans_err = 0.01, alpha = 1, beta = 0.5, RTT = 600,RTT0 = 25, packet_err = 0.05, new = True):
        self.N = N # number of states
        self.C = C
        self.W = C * RTT # in Mbyte, W = C*RTT where C is the maximum bandwidth
        self.RTT = RTT # in ms 
        self.RTT0 = RTT0
        self.rho = max([RTT/RTT0,1])
        self.trans_err = trans_err # a.k.a. lambda 
        self.packet_err = packet_err # percentage chance that a packet is lost
        self.alpha = alpha # window growth rate
        self.beta = beta # window reduction ratio
        self.a = (np.linspace(1,N,N)-0.5)*self.W/N # Discretisation of the window-size
        self.pi = np.ones(N)/N # Stationnary Distribution
        self.P = None # Transition Probability Matrix
        self.S = None
        self.tau = None 
        self.x = 0 # Steady State average throughput
        self.new = new

    def D(self,x,y): # growth time in seconds
        return self.RTT0**2*(y-self.beta*x/self.RTT)
    
    def avg_w(self,t,x): # packets sent every RTT
        return self.beta*x+self.rho*t*self.alpha/(2*self.RTT0)

    def transition_proba_Hybla(self,i,j): # Probability to transition from i to j
        l = self.trans_err
        if j < self.beta*(i-0.5):
            return 0
        if j == self.N:
            sum = 0
            for jp in np.linspace(1,self.N-1,self.N-1):
                if jp >= self.beta*(i-0.5):
                    sum += self.transition_proba_Hybla(i,jp)
            return 1-sum
        return np.exp(-l* np.maximum(self.D(self.a[i-1],(j-1)*self.W/self.N),0)) - np.exp(-l* self.D(self.a[i-1],j*self.W/self.N))
    
    def transition_proba_Hybla_new(self,i,j):
        if j < self.beta*(i-0.5):
            return 0
        if j == self.N:
            sum = 0
            for jp in np.linspace(1,self.N-1,self.N-1):
                if jp >= self.beta*(i-0.5):
                    sum += self.transition_proba_Hybla_new(i,jp)
            return 1-sum
        t_min = np.maximum(self.D(self.a[i-1],(j-1)*self.W/self.N),0)
        t_max = self.D(self.a[i-1],j*self.W/self.N)

        return np.exp(-t_min*self.RTT/(self.avg_w(t=t_min,x=self.a[i-1])*self.packet_err))-np.exp(-t_max*self.RTT/(self.avg_w(t=t_max,x=self.a[i-1])*self.packet_err))
    
    def compute_stationnary_distribution(self):
        # 1. Compute the transition probability Matrix P
        self.P = np.zeros([self.N,self.N])
        for i in range(self.N):
            for j in range(self.N):
                if j == self.N:
                    self.P[i,j] = 1-np.sum(self.P[i,:-1])
                    continue
                if not self.new:
                    self.P[i,j] = self.transition_proba_Hybla(i+1,j+1)
                else:
                    self.P[i,j] = self.transition_proba_Hybla_new(i+1,j+1)

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
                self.S[i,j] = (self.a[i]*self.tau[i,j]*self.beta + self.rho*self.tau[i,j]**2/(2*self.RTT0))/self.RTT
        return
    
    def avg_throughput(self):
        if self.P == None:
            self.compute_stationnary_distribution()
        self.compute_tau_and_S()
        numerator = np.dot(self.pi,np.dot(self.P,np.transpose(self.S)).diagonal())
        denominator = np.dot(self.pi,np.dot(self.P,np.transpose(self.tau)).diagonal())
        self.x = 1/self.C*numerator/denominator
        return self.x
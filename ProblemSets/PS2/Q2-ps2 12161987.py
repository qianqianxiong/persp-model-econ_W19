#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#Problem 2.1


# In[1]:


def integ_g_x(g, a, b, N, methods ={'midpoint','trapezoid','Simpsons'} ):
        inter_mid = [a + (b-a) *(2*i+1)/(2 * N) for i in range(N)]
        appro_mid=sum(g(i) for i in inter_mid)*(b-a)/N
        inter_trap=[a+i*(b-a)/N for i in range (N+1)]
        appro_trap=(2*sum(g(i) for i in inter_trap[1:-1])+ g(inter_trap[0])+g(inter_trap[-1]))*(b-a)/(N*2)
        inter_sims=[a+i*(b-a)/(N*2) for i in range(2*N+1)]
        appro_sims=(g(inter_sims[0])+g(inter_sims[-1])+4*sum(g(i) for i in inter_sims[1:-1:2])+2*sum(g(i) for i in inter_sims[2:-2:2]))*(b-a)/(6*N)
        return appro_mid, appro_trap, appro_sims


# In[2]:


g =lambda x: 0.1 * (x ** 4) - 1.5 * (x ** 3) + 0.53 * (x ** 2) + 2 * x + 1
value = integ_g_x(g,-10, 10, 2000)
print("Midpoint method value:",value[0])
print("Trapezoid method value:",value[1])
print("Simpsons method value:",value[2])
print("Difference between Midpoint value and true value:", abs(value[0] - 4373 - 1/3))
print("Difference between Trapezoid value and true value:", abs(value[1] - 4373 - 1/3))
print("Difference between Simpsons value and true value:", abs(value[2] - 4373 - 1/3))


# In[ ]:


#Problem 2.2


# In[3]:


import scipy.stats as sts
import pandas as ps
import numpy as np

def normal_dis(mean, std, N, k):
    nodes = np.linspace(mean - std * k, mean + std * k, N)
    mids = (nodes[1:] + nodes[:-1]) / 2
    w0 = sts.norm.cdf(mids[0], mean, std)
    w1 = sts.norm.cdf(mids, mean, std)
    w1 = np.diff(w1).tolist()
    w2 = 1 - sts.norm.cdf(mids[-1], mean, std)
    weights = [w0] + w1 + [w2]
    nodes = nodes.tolist()
    table = {'Nodes' : nodes, 'Weights' : weights}
    psf = ps.DataFrame(table)
    return psf


# In[4]:


normal_dis(0,1,11,3)


# In[ ]:


#Problem 2.3


# In[5]:


def normal_log(mean, std, N, k):
    nodes = np.linspace(mean-std*k, mean+std*k, N)
    mids = (nodes[1:] + nodes[:-1]) / 2
    w1 = sts.norm.cdf(mids[0], mean, std)
    w2 = sts.norm.cdf(mids, mean, std)
    w2 = np.diff(w2).tolist()
    w3 = 1 - sts.norm.cdf(mids[-1], mean, std)
    weights = [w1] + w2 + [w3]
    nodes = np.exp(nodes).tolist()
    table = {'Nodes' : nodes, 'Weights' : weights}
    psf = ps.DataFrame(table)    
    return psf


# In[6]:


normal_log(0,1,11,3)


# In[ ]:


#Problem 2.4


# In[7]:


appro_log=normal_log(10.5,0.8,19,3)
appro = sum(appro_log.Nodes * appro_log.Weights)
exact=np.exp(10.5+0.8*0.8/2)
print("The exact value is:",exact)
print("Approximation is:",appro)
print("Difference between approximation and the exact value is :",abs(appro-exact))


# In[ ]:


#Problem 3.1


# In[8]:


g=lambda x: 0.1 * x ** 4 - 1.5 * x ** 3 + 0.53 * x ** 2 + 2 * x + 1 
from scipy import optimize
def h_x(x, a=-1, b=1):

    return [x[0] + x[1] + x[2] - b + a,
           x[0] * x[3] + x[1] * x[4] + x[2] * x[5] - 0.5*b**2 + 0.5*a**2,
           x[0] * x[3]**2 + x[1] * x[4]**2 + x[2] * x[5]**2 - (1/3)*b**3 + (1/3)*a**3,
           x[0] * x[3]**3 + x[1] * x[4]**3 + x[2] * x[5]**3 - (1/4)*b**4 + (1/4)*a**4,
           x[0] * x[3]**4 + x[1] * x[4]**4 + x[2] * x[5]**4 - (1/5)*b**5 + (1/5)*a**5,
           x[0] * x[3]**5 + x[1] * x[4]**5 + x[2] * x[5]**5 - (1/6)*b**6 + (1/6)*a**6]


# In[9]:


array = optimize.root(h_x, [0, 0, 0, 0, 0, 0])
array.x


# In[10]:


w1=array.x[0]
w2=array.x[1]
w3=array.x[2]
x1=array.x[3]
x2=array.x[4]
x3=array.x[5]
print(w1,w2,w3,x1,x2,x3) 


# In[11]:


a=-10
b=10
p = (b - a)/2
q = (b + a)/2
appro=(b - a) * 0.5 * (w1 * g(p*x1+q)+w2*g(p*x2+q)+w3*g(p*x3+q))
exact=4373+1/3
print("Approximation:",appro)
print("the exact value:",exact)
print("Absolute Error:",abs(exact-appro))


# In[ ]:


#Problem 3.2


# In[12]:


from scipy.integrate import quad
g=lambda x: 0.1 * x ** 4 - 1.5 * x ** 3 + 0.53 * x ** 2 + 2 * x + 1
appro=quad(g,-10,10)[0]
exact=4373+1/3
print("Approximated integral using Scipy Gaussian quadrature:",appro)
print("the exact value:",exact)
print("Absolute Error:",abs(exact-appro))


# In[ ]:


#Problem 4.1


# In[13]:


def cir_inside(x, y):
    if x ** 2 + y ** 2 <= 1:
        return 1
    else:
        return 0


# In[14]:


import numpy as np
np.random.seed=25
N = 1
while N > 0:
    mc_x = np.random.uniform(-1, 1, N)
    mc_y = np.random.uniform(-1, 1, N)
    frac = list(map(cir_inside, mc_x, mc_y))
    area = 4 * np.mean(frac)
    if round(area, 4) == 3.1415:
        break
    else:
        N += 1


# In[15]:


print('The smallest number of random draws:', N)
print('Circle Area:',area)


# In[ ]:


#Problem 4.2


# In[16]:


def isPrime(n):
    for i in range(2, int(np.sqrt(n) + 1)):
        if n % i == 0:
            return False

    return True


# In[17]:


def primes_ascend(N, min_val=2):
    
    primes_vec = np.zeros(N, dtype=int)
    MinIsEven = 1 - min_val % 2
    MinIsGrtrThn2 = min_val > 2
    curr_prime_ind = 0
    if not MinIsGrtrThn2:
        i = 2
        curr_prime_ind += 1
        primes_vec[0] = i
    i = min(3, min_val + (MinIsEven * 1))
    while curr_prime_ind < N:
        if isPrime(i):
            curr_prime_ind += 1
            primes_vec[curr_prime_ind - 1] = i
        i += 2

    return primes_vec


# In[18]:


def equi_dis_nth(n, d, seq_name):
    prime_vec = primes_ascend(d)
    if seq_name == 'weyl':
        seq = np.sqrt(prime_vec) * n
        
    elif seq_name == 'haber':
        seq = np.sqrt(prime_vec) * (n * (n + 1) / 2)
    
    elif seq_name == 'niederreiter':
        seq = [n * 2 ** (i / (n + 1)) for i in range(1, d+1)]
        seq = np.array(seq)
    
    elif seq_name == 'baker':
        seq = n * np.exp(prime_vec)
    
    sequence_name = seq - np.floor(seq)
    return sequence_name       


# In[19]:


equi_dis_nth(5, 2, 'weyl')


# In[20]:


equi_dis_nth(5, 2, 'haber')


# In[21]:


equi_dis_nth(5, 2, 'niederreiter')


# In[22]:


equi_dis_nth(5, 2, 'baker')


# In[ ]:


#Problem 4.3


# In[23]:


import numpy as np
np.random.seed = 25
N = 1
while N >= 0:
    x= [2 * equi_dis_nth(i, 2, 'weyl')[0] - 1 for i in range(N)]
    y= [2 * equi_dis_nth(i, 2, 'weyl')[1] - 1 for i in range(N)]
    frac = list(map(cir_inside,x,y ))
    area = 4 * np.mean(frac)
    if round(area, 4) == 3.1415:
        print('Smallest number of draws for Weyl:', N)
        break
    else:
        N += 1


# In[24]:


N = 1
while N >= 0:
    x = [2 * equi_dis_nth(i, 2, 'haber')[0] - 1 for i in range(N)]
    y = [2 * equi_dis_nth(i, 2, 'haber')[1] - 1 for i in range(N)]
    frac = list(map(cir_inside, x, y))
    area = 4 * np.mean(frac)
    if round(area, 4) == 3.1415:
        print('Smallest number of draws for Haber:', N)
        break
    else:
        N += 1


# In[ ]:


N = 1
while N >= 0:
    x = [2 * equi_dis_nth(i, 2, 'niederreiter')[0] - 1 for i in range(N)]
    y = [2 * equi_dis_nth(i, 2, 'niederreiter')[1] - 1 for i in range(N)]
    frac = list(map(cir_inside, x, y))
    area = 4 * np.mean(frac)
    if round(area, 4) == 3.1415:
        print('Smallest number of draws for Niederreiter:', N)
        break
    elif N > 10000:
        print('Niederreiter not converging')
    else:
        N += 1


# In[ ]:


N = 1
while N >= 0:
    x = [2 * equi_dis_nth(i, 2, 'baker')[0] - 1 for i in range(N)]
    y = [2 * equi_dis_nth(i, 2, 'baker')[1] - 1 for i in range(N)]
    frac = list(map(cir_inside, x, y))
    area = 4 * np.mean(frac)
    if round(area, 4) == 3.1415:
        print('Smallest number of draws for Baker:', N)
        break
    else:
        N += 1


# In[ ]:





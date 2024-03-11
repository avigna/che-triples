#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import math
import random
import h5py

from os import path
import multiprocessing
import scipy as sp
import psutil
from mpi4py import MPI

import myconstants as const


# In[2]:


#initial and final values of distributions

filename = "inipop_triple.hdf5"

CONST_R_SUN = 0.004649130343817401

N_sample = 10000 #population strength 
N_CPU = 4 
nbin_p = nbin_mass = nbin_q = 1000 # number of sampling for rejection sampling

p_i=0.2 #initial value of log p(days) for sampling
p_f=7.9 #final value of log p(days) for sampling

m_i=8 #initial value of m(in solar masses) for sampling
m_f=100 #final value of m(in solar masses) for sampling

q_i=0.1 #initial value of q and q_out for sampling
f_q_in = 1 #final value of q_in for sampling
q_f=10 #final value of q_out for sampling (from Tokovinin MSC data Sept 2021)
lambda_par =1.05 # exponent for declining exponentional (estimated by fitting Tokovinin MSC data Sept 2021)

m_lower_limit = 0.08
m_upper_limit = 100

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


# In[3]:


#mass-distribution kroupa et al 2001
#Rejection sampling for mass
#Initial Mass Function

def mass_sampling(nbin_mass,m_i,m_f):
    def IMF(m,m_i,m_f):
        alpha1 = -1.3
        alpha2 = -2.3
        alpha3 = -2.3
        m1 = 0.08
        m2 = 0.5
        m3 = 1.0
        m4 = 100.0
        C1 = 0.25293414013708077
        C2 = 0.1264670700685404
        C3 = 0.1264670700685404
        if 0.08 < m <= 0.5:
            return C1*(m**(alpha1))
        elif 0.5 < m <= 1:
            return C2*(m**(alpha2))
        else:
            return C3*(m**(alpha3)) 

    m = np.linspace(m_i,m_f,nbin_mass) # random-sample m to input in the function
    pdf_m=[]
    for i in range(len(m)):
        pdf_m.append(IMF(m[i],m_i,m_f))
    k= (max(pdf_m))  #maximum of the function

    x = np.random.uniform(m_i,m_f) # random sample (one per iteration) variable
    u = np.random.uniform(0,k)  #random sample pdf
    #print(u)
    #print(IMF(n_star,x))
    while u >= IMF(x,m_i,m_f):
        x = np.random.uniform(m_i,m_f) # random sample (one per iteration) variable
        u = np.random.uniform(0,k)  #random sample pdf
        if u <= IMF(x,m_i,m_f):
            #print(u,x)
            break
        else:
            continue
    return x


# In[4]:


# period distribution moe di stefano et al 2017

# moe di stefano et al 2017 function to sample period
def Ipf(m,period):
    def p_lessthan1(m):  
        log10_M1 = np.log10(m)
        return 0.02 + 0.04*log10_M1 + 0.07*log10_M1**2
    def p_equals2pt7(m):
        log10_M1 = np.log10(m)
        return 0.039 + 0.07*log10_M1 + 0.01*log10_M1**2
    def p_equals5pt5(m):
        log10_M1 = np.log10(m)
        return 0.078 - 0.05*log10_M1 + 0.04*log10_M1**2
        
    delP=0.7 #moe di stefano et al 2017
    alpha = 0.018 #moe di stefano et al 2017
   
    IPF=[]
    P=[]
    for p in period:
        if 0.2<=p<1.0:
            IPF.append(p_lessthan1(m))
            P.append(p)
        elif 1.0<=p<(2.7-delP):
            IPF.append(p_lessthan1(m)+(((p-1)/(1.7-delP))*(p_equals2pt7(m)-p_lessthan1(m)-alpha*delP)))
            P.append(p)
        elif (2.7-delP)<=p<(2.7+delP):
            IPF.append(p_equals2pt7(m)+(alpha*(p-2.7)))
            P.append(p)
        elif (2.7+delP)<=p<5.5:
            IPF.append(p_equals2pt7(m)+(alpha*delP)+(((p-2.7-delP)/(2.8-delP))*(p_equals5pt5(m)-p_equals2pt7(m)-(alpha*delP))))
            P.append(p)
        elif 5.5<=p<8.0:
            IPF.append(p_equals5pt5(m)*np.exp(-0.3*(p-5.5)))
            P.append(p)

        return IPF,P 


# In[5]:


#rejection sampling for period - Ref : moe di stefano et 2017

def period_sampling(nbin_period,p_i,p_f,mass):
    period = np.linspace(p_i,p_f,nbin_period)
    
    # moe di stefano et al 2017 function to sample period
    def Ipf(m,period):
        def p_lessthan1(m):
            log10_M1 = np.log10(m)
            return 0.02 + 0.04*log10_M1 + 0.07*log10_M1**2
        def p_equals2pt7(m):
            log10_M1 = np.log10(m)
            return 0.039 + 0.07*log10_M1 + 0.01*log10_M1**2
        def p_equals5pt5(m):
            log10_M1 = np.log10(m)
            return 0.078 - 0.05*log10_M1 + 0.04*log10_M1**2

        delP=0.7 #moe di stefano et al 2017
        alpha = 0.018 #moe di stefano et al 2017
        
        IPF=[]
        P=[]
        for p in period:
            if 0.2<=p<1.0:
                IPF.append(p_lessthan1(m))
                P.append(p)
            elif 1.0<=p<(2.7-delP):
                IPF.append(p_lessthan1(m)+(((p-1)/(1.7-delP))*(p_equals2pt7(m)-p_lessthan1(m)-alpha*delP)))
                P.append(p)
            elif (2.7-delP)<=p<(2.7+delP):
                IPF.append(p_equals2pt7(m)+(alpha*(p-2.7)))
                P.append(p)
            elif (2.7+delP)<=p<5.5:
                IPF.append(p_equals2pt7(m)+(alpha*delP)+(((p-2.7-delP)/(2.8-delP))*(p_equals5pt5(m)-p_equals2pt7(m)-(alpha*delP))))
                P.append(p)
            elif 5.5<=p<8.0:
                IPF.append(p_equals5pt5(m)*np.exp(-0.3*(p-5.5)))
                P.append(p)

            return IPF,P 
    k=max(Ipf(mass,period)[0]) #maximum of the functiom
    #print(k)

    x = np.random.uniform(p_i,p_f) # random sample (one per iteration) variable
    #print(x)
    u = np.random.uniform(0,k)  # random sample pdf (one per iteration) 

    while u >= Ipf(mass,[x])[0][0]:
        x = np.random.uniform(p_i,p_f) # random sample (one per iteration) variable
        u = np.random.uniform(0,k)  #random sample pdf
        if u <= Ipf(mass,[x])[0][0]:
            #print(u,x)
            break
        else:
            continue
    return x  


# In[6]:


# mass_ratio distribution - Ref: moe di stefano et al 2017

#exponent as a function of period

# gamma small - exponent in the q_range (0.1<q<0.3)
def gamma_small_mass_range1(period):
    return 0.3
def gamma_small_mass_range2(period):
    if 0.2<=period<2.5:
        return 0.2
    elif 2.5<=period<5.5:
        return 0.2-0.3*(period-2.5)
    elif 5.5<=period<8:
        return -0.7-0.2*(period-5.5)
def gamma_small_mass_range3(period):
    if 0.2<=period<1:
        return 0.1
    elif 1<=period<3.0:
        return (0.1-0.15*(period-1))
        Period_1=(period)
    elif 3<=period<5.6:
        return (-0.2-0.5*(period-3))
    elif 5.6 <= period<8:
        return -1.5

# gamma large - exponent in the q_range (0.3<q<1)
def gamma_large_mass_range1(period):
    if 0.2 <= period < 5.0:
        return -0.5
    elif  5.0<= period < 8:
        return -0.5 - 0.3 * (period-5.0)
def gamma_large_mass_range2(period):
    if 0.2<=period<1:
        return -0.5
    elif 1<=period<4.5:
        return -0.5 - 0.2 * (period-1)
    elif 4.5<=period<6.5:
        return -1.2 - 0.4 * (period-4.5)
    elif 6.5 <=period<8:
        return -2.0   
def gamma_large_mass_range3(period):
    if 0.2<=period<1:
        return -0.5
    elif 1<=period<2.0:
        return -0.5 - 0.9 * (period-1)
    elif 2.0<=period<4:
        return -1.4 - 0.3 * (period-2)
    elif 4<=period<8:
        return -2.0


# In[7]:


#exponent as a function of period and mass(interpolated linearly in the unknown mass ranges) - ref : Moe di stefano 2017 

# gamma small - exponent in the q_range (0.1<q<0.3)
def moe_di_stefano_gamma_small(M1,period):
    
    if M1 >= 0.08 and M1 <= 1.2:
        return gamma_small_mass_range1(period)
    if M1 > 1.2 and M1 < 3.5:
        y1 = gamma_small_mass_range1(period)
        y2 = gamma_small_mass_range2(period)
        x1 = 1.2
        x2 = 3.5
        return y1 + (M1 - x1) * (y2-y1)/(x2-x1)
    if M1 >= 3.5 and M1 < 6.0:
        y1 = gamma_small_mass_range2(period)
        y2 = gamma_small_mass_range3(period)
        x1 = 3.5
        x2 = 6.0
        return y1 + (M1 - x1) * (y2-y1)/(x2-x1)
    
    if M1 > 6.0:
        return gamma_small_mass_range3(period)
    
# gamma large - exponent in the q_range (0.3<q<1)    
def moe_di_stefano_gamma_large(M1,period):
    
    if M1 >= 0.08 and M1 <= 1.2:
        return gamma_large_mass_range1(period)
    if M1 > 1.2 and M1 < 3.5:
        y1 = gamma_large_mass_range1(period)
        y2 = gamma_large_mass_range2(period)
        x1 = 1.2
        x2 = 3.5
        return y1 + (M1 - x1) * (y2-y1)/(x2-x1)
    if M1 >= 3.5 and M1 < 6.0:
        y1 = gamma_large_mass_range2(period)
        y2 = gamma_large_mass_range3(period)
        x1 = 3.5
        x2 = 6.0
        return y1 + (M1 - x1) * (y2-y1)/(x2-x1)
    
    if M1 > 6.0:
        return gamma_large_mass_range3(period)


# In[8]:


# plot to check the interpolation

'''
m=np.linspace(1,6.5,1000)
period=np.linspace(0.2,7.9,1000)

gamma_small=[]
gamma_large=[]
for i in range(1000):
    gamma_small.append(moe_di_stefano_gamma_small(m[i],period[i]))
    gamma_large.append(moe_di_stefano_gamma_large(m[i],period[i]))

plt.plot(m,gamma_small)
plt.plot(m,gamma_large)

'''


# In[9]:


# moe di stefano et al 2017 function to determine fraction of twin samples(q = 0.95 - 1 )
def q_twin(m,p):   
    def P_twin(m):
        if m <= (6.5):
            return (8 - m)
        else:
            return  1.5
#twin fraction for masses less than one
    def F_twin_Plessthan1(m): 
        return (0.30 - (0.15*(np.log10(m))))
#to define twin fraction as a function of mass and period
    def F_twin(m,p):
        if p < 1:
            return F_twin_Plessthan1(m)
        elif 1<=p<P_twin(m):
            return (F_twin_Plessthan1(m)*(1-((p-1)/(P_twin(m)-1))))
        elif p >= P_twin(m):
            return 0.0
    return F_twin(m,p)


# In[10]:


# Estimate the correction factor to normalize p_distribution over the complete q range. Ref : moe di stefano et al 2017

# constrcuting a pdf function for the period
def second_period(mass,period): # pdf of p in the complete q_range
    g1= moe_di_stefano_gamma_small(mass,period)
    g2 = moe_di_stefano_gamma_large(mass,period)
    F_twin=q_twin(mass,period)
    correction_factor=(((0.3**(g2-g1))*(((0.3**(g1+1))-(0.1**(g1+1)))/(g1+1)))+(((1-(0.3**(g2+1)))/(g2+1))*(1/(1-F_twin))))/(((1-0.3**(g2+1))/(g2+1))*(1/(1-F_twin)))
    return Ipf(mass,[period])[0][0]*correction_factor


# In[11]:


#check plot

'''

a=[]
b=[]
c=[]
d=[]
e=[]
period=np.linspace(0.2,7.9,1000)
for i in range(1000):
    mass=1
    a.append(second_period(mass,period[i]))
    mass=3
    b.append(second_period(mass,period[i]))
    mass=5
    c.append(second_period(mass,period[i]))
    mass=10
    d.append(second_period(mass,period[i]))
    mass=27vi: error while loading shared libraries: libffi.so.6: cannot open shared object file: No such file or directory
    e.append(second_period(mass,period[i]))
    
plt.plot(period,a,'--',label="mass = 1")
plt.plot(period,b,'--',label="mass = 3")
plt.plot(period,c,'--',label="mass = 5")
plt.plot(period,d,'--',label="mass = 10")
plt.plot(period,e,'--',label="mass = 27")
plt.legend()
plt.xlabel('log(P) in days')
plt.ylabel('P(p)')
plt.savefig('Period_distribution')

'''


# In[12]:


#p sampling using the pdf after correction factor

def p_qgreaterthanpt1(mass,p_i,p_f,nbin_p):
    period=np.linspace(p_i,p_f,nbin_p)
    P_qpt1=[]
    for i in range(len(period)):
        P_qpt1.append(second_period(mass,period[i]))
    k=max(P_qpt1) #maximum of the functiom
    #print(k)

    x = np.random.uniform(p_i,p_f) # random sample (one per iteration) variable
    #print(x)
    u = np.random.uniform(0,k)  # random sample pdf (one per iteration) 

    while u >= second_period(mass,x):
        x = np.random.uniform(p_i,p_f) # random sample (one per iteration) variable
        u = np.random.uniform(0,k)  #random sample pdf
        if u <= second_period(mass,x):
            #print(u,x)
            break
        else:
            continue
    return x


# In[13]:


#check plot

'''
nbin_p=1000
P_solar=[]
P_Atype=[]
P_MidB=[]
for i in period:
    mass=1
    P_solar.append(p_qgreaterthanpt1(mass,p_i,p_f,nbin_p))
    mass=4
    P_Atype.append(p_qgreaterthanpt1(mass,p_i,p_f,nbin_p))
    mass=8
    P_MidB.append(p_qgreaterthanpt1(mass,p_i,p_f,nbin_p))
    
plt.plot(period,a,'--',color='red',label='Solar')
plt.plot(period,b,'--',color='orange',label='Atype')
plt.plot(period,c,'--',color='green',label='MidB')
plt.hist(P_solar,histtype='step',density='True',color='red')
plt.hist(P_Atype,histtype='step',density='True',color='orange')
plt.hist(P_MidB,histtype='step',density='True',color='green')
plt.legend()
plt.title('p_greaterthanpt1')
plt.xlabel('LogP in days')
plt.ylabel('P(Period)')
plt.savefig('p_greaterthanpt1_distribution')

'''


# In[14]:


#mass_ratio sampling Ref : Moe di stefano et al 2017

def mass_ratio_sampling(mass,period,q_i,q_f,nbin_q):
    
    def mass_ratio(q,mass,period):
        g1=moe_di_stefano_gamma_small(mass,period)
        g2=moe_di_stefano_gamma_large(mass,period)
        F_twin=q_twin(mass,period)
        c2=1/(((0.3)**(g2-g1))*(((0.3**(g1+1))-(0.1**(g1+1)))/(g1+1))+((1-(0.3**(g2+1)))/(g2+1))+((F_twin/(1-F_twin))*(1-(0.3**(g2+1)))/(g2+1)))
        c1=c2*(0.3**(g2-g1))
        c3=c2*((F_twin/(1-F_twin))*((1-(0.3**(g2+1)))/(1-(0.95**(g2+1)))))
        if 0.1<=q<=0.3:
            return c1*(q**g1)
        elif 0.3<q<0.95:
            return c2*(q**g2)
        elif 0.95<=q<=1:
            return c3*(q**g2)

    
    q=np.linspace(q_i,q_f,nbin_q)
    q_arr =[]
    for i in q:
        q_arr.append(mass_ratio(i,mass,period))
    k=max(q_arr)
    
    x = np.random.uniform(q_i,q_f) # random sample (one per iteration) variable
    #print(x)
    u = np.random.uniform(0,k)  # random sample pdf (one per iteration) 

    while u >= mass_ratio(x,mass,period):
        x = np.random.uniform(q_i,q_f) # random sample (one per iteration) variable
        u = np.random.uniform(0,k)  #random sample pdf
        if u <= mass_ratio(x,mass,period):
            #print(u,x)
            break
        else:
            continue
    return x    


# In[15]:


#mass_ratio sampling for q_out by extrapolating Moe di stefano et al 2017 for q greater than 1 values 

def mass_ratio_sampling_MoediStefano_extrapolation(mass,period,q_i,q_f,nbin_q):
    def mass_ratio(q,mass,period):
        g1=moe_di_stefano_gamma_small(mass,period)
        #print(g1)
        g2=moe_di_stefano_gamma_large(mass,period)
        #print(g2)
        c1 = 1/(((0.3**(g1+1) - 0.1**(g1+1))/(g1+1))+(0.3**(g1-g2) * ((q_f**(g2+1))-(0.3**(g2+1)))/(g2+1)))
        #print(c1)
        c2 = c1*(0.3**(g1-g2))
        #print(c2)
        if 0.1<=q<=0.3:
            return c1*(q**g1)
        elif 0.3<q<=q_f:
            return c2*(q**g2)

    
    q=np.linspace(q_i,q_f,nbin_q)
    q_arr =[]
    for i in q:
        #print(mass_ratio(i,mass,period))
        q_arr.append(mass_ratio(i,mass,period))
    k=max(q_arr)
    
    x = np.random.uniform(q_i,q_f) # random sample (one per iteration) variable
    #print(x)
    u = np.random.uniform(0,k)  # random sample pdf (one per iteration) 
    o=1
    while u >= mass_ratio(x,mass,period):
        x = np.random.uniform(q_i,q_f) # random sample (one per iteration) variable
        u = np.random.uniform(0,k)  #random sample pdf
        if u <= mass_ratio(x,mass,period):
            #print(u,x)
            break
        else:
            o = o+1
            #print(o)
            continue
    return x,mass_ratio(x,mass,period)    


# In[16]:

#check plot 

'''

a=[]
period=np.linspace(0.2,7.9,1000)
for i in range(1000):
    mass=1
    a.append(q_greater_than1_sampling(mass,period[i],q_i,q_f,nbin_q,lambda_par=1.05)[0])
    
plt.hist(a,histtype='step',density='True')
plt.xlabel('Mass_ratio')
plt.ylabel('P(q)')

'''


# In[18]:


#eccentricity sampling Ref : Moe di stefano et al 2017
def moe_di_stefano_sample_e(log10_M1,log_P):
    
    e_max = moe_di_stefano_sample_e_e_max(log_P)
    
    eta = moe_di_stefano_sample_e_eta(log10_M1,log_P)
    if (eta < -5.0): # this is done to avoid numerical error
        eta = -5.0
    
    etap1 = eta + 1.0
    
    y = random.random()

    e1 = 1.0e-10
    e2 = e_max

    e = pow( y*(pow(e2,etap1) - pow(e1,etap1)) + pow(e1,etap1), 1.0/etap1)
    #e = e_max * pow(y, 1.0/etap1)

    if (e>=1.0):
        print("ERROR in moe_di_stefano_sample_e")
        print("emax",e_max,"eta",eta,e)
        exit(0)
    return e,eta,e_max
    
def moe_di_stefano_sample_e_e_max(log_P):
    P_day = pow(10.0,log_P)
    if P_day < 2.0:
        return 0.0
    else:
        return 1.0 - pow(P_day/2.0,-2.0/3.0)
    
def moe_di_stefano_sample_e_eta(log10_M1,log_P):
    M1 = pow(10.0,log10_M1)

    if M1 >= 0.08 and M1 <= 3.0:
        return moe_di_stefano_sample_e_aux1(log_P)
    if M1 > 3.0 and M1 <= 7.0:
        y1 = moe_di_stefano_sample_e_aux1(log_P)
        y2 = moe_di_stefano_sample_e_aux2(log_P)
        x1 = 3.0
        x2 = 7.0
        return y1 + (M1 - x1) * (y2-y1)/(x2-x1)
    if M1 > 7.0:
        return moe_di_stefano_sample_e_aux2(log_P)
    
def moe_di_stefano_sample_e_aux1(log_P):
    return 0.6 - 0.7/(log_P - 0.5)

def moe_di_stefano_sample_e_aux2(log_P):
    return 0.9 - 0.2/(log_P - 0.5)


# In[19]:


# check plot 
'''
ecc =[]
for i in range(1000):
    mass = mass_sampling(nbin_mass,m_i,m_f)
    period = p_qgreaterthanpt1(mass,p_i,p_f,nbin_p)
    #mass_ratio = mass_ratio_sampling(m1,p_in,q_i,1,nbin_q)
    log10_M1 = np.log10(mass)
    log_P = period
    #print(log10_M1,log_P)
    ecc.append(moe_di_stefano_sample_e(log10_M1,log_P)[0])
    #ecc.append(e)
plt.hist(ecc)

'''


# In[20]:


#AP and LAN sampling 
AP_lower = LAN_lower = 0.01*np.pi/180.0
AP_upper = LAN_upper = 359.99*np.pi/180.0


# In[21]:


#convert period in seconds to separations in cms
def separation(m,T):
    return ((const.CGS_G*m*const.CGS_MSUN)*((T)**2)/(4*(np.pi**2)))**(1./3)


# In[22]:


#inclination distribution
def inclination():
    return math.degrees(np.arccos(np.random.uniform(np.cos(0),np.cos(math.radians(180)))))


# In[23]:


#stability criteria

def stability_criterion_MA01(m1,m2,m3,a_in,a_out,e_out,i_rel):
    q_out = m3/(m1+m2)
    rp_out = a_out*(1.0 - e_out)
    rp_out_crit = a_in * 2.8*pow( (1.0 + q_out)*(1.0 + e_out)/np.sqrt(1.0-e_out),2.0/5.0 ) * (1.0 - 0.3*i_rel/np.pi)
    if (rp_out < rp_out_crit):
        stable = False
    else:
        stable = True
    
    return stable


# In[24]:


def estimate_MS_radius_RSun(M_MSun):
    return pow(M_MSun,0.7)


# In[25]:


def compute_Eggleton_fq(q):
    ### q is defined as m_primary/m_secondary ###
    q_pow_one_third = pow(q,1.0/3.0)
    q_pow_two_third = q_pow_one_third*q_pow_one_third
    return 0.49*q_pow_two_third/(0.6*q_pow_two_third + np.log(1.0 + q_pow_one_third))

def compute_mutual_inclination(INCL_k,INCL_l,LAN_k,LAN_l):
    cos_INCL_rel = np.cos(INCL_k)*np.cos(INCL_l) + np.sin(INCL_k)*np.sin(INCL_l)*np.cos(LAN_k-LAN_l)
    return np.arccos(cos_INCL_rel)


import sklearn

from triple_stability.classify_trip import mlp_classifier
# from quadruple_stability.classify_quad_3p1 import mlp_classifier_3p1

mlp_triple = "./triple_stability/mlp_model_trip_ghost_v1.0.2.pkl"
# mlp_3p1_pfile = "./quadruple_stability/classify_quad_3p1.py/mlp_model_3p1_ghost.pkl"


# In[26]:


pop=[]
print("len_pop",len(pop))

sample_size = range(N_sample)

def pop_sample(i):
    np.random.seed(4*i)
    print(i)
    #print(i,multiprocessing.current_process()._identity[0])
    sampled = False
    while sampled == False:
        m1 = mass_sampling(nbin_mass,m_i,m_f)
        p_in = p_qgreaterthanpt1(m1,p_i,p_f,nbin_p)
        a_in = separation(m1,pow(10.0,p_in)*86400)/const.CGS_AU
        LAN_in = np.random.uniform(LAN_lower,LAN_upper)
        LAN_out = np.random.uniform(LAN_lower,LAN_upper)
        AP_in = np.random.uniform(AP_lower,AP_upper)
        AP_out = np.random.uniform(AP_lower,AP_upper)
        q_in = mass_ratio_sampling(m1,p_in,q_i,f_q_in,nbin_q)
        log10_M1 = np.log10(m1)
        log_P = p_in
        e_in = moe_di_stefano_sample_e(log10_M1,log_P)[0]
        m2 = q_in * m1
        rp1 = a_in*(1.0-e_in)
        RL1 = rp1*compute_Eggleton_fq(q_in)
        RL2 = rp1*compute_Eggleton_fq(1.0/q_in)
        R1 = estimate_MS_radius_RSun(m1)*CONST_R_SUN
        R2 = estimate_MS_radius_RSun(m2)*CONST_R_SUN
        #print(R1,RL1,R2,RL2)
        p_out = p_qgreaterthanpt1(m1+m2,p_i,p_f,nbin_p)
        log_P = p_out
        #print('p_out',p_out)
        a_out = separation(m1+m2,pow(10.0,p_out)*86400)/const.CGS_AU
        q_out = mass_ratio_sampling_MoediStefano_extrapolation(m1+m2,p_out,q_i,q_f,nbin_q)[0]
        e_out =moe_di_stefano_sample_e(log10_M1,log_P)[0]
        #print('e',e_out)
        m3=q_out*(m1+m2)
        i_1=inclination()
        i_2=inclination()
        i_rel = compute_mutual_inclination(i_1,i_2,LAN_in,LAN_out)
        
        al = a_in/a_out
        
        if mlp_classifier(mlp_triple, q_in, q_out, al, e_in, e_out, i_rel) == True:
            if R1 < RL1 or R2 < RL2:
                if m_upper_limit >= m2 > m_lower_limit and m_upper_limit >= m3 > m_lower_limit:
                    sampled = True
        else:
            sampled = False
            print("false :(",m2,m3)
        if sampled == True:
            return([m1,m2,m3,a_in,a_out,e_in,e_out,i_1,i_2,LAN_in,LAN_out,AP_in,AP_out,q_in,q_out])
            break
        else:
            continue  
            
#print(pop_sample(1))

pop=[]
print("len_pop",len(pop))

if __name__ == '__main__':

    pool = multiprocessing.Pool(processes=N_CPU)

    result_list = pool.map(pop_sample, sample_size)

    print(len(result_list))

initial_pop = np.array(result_list)

m1s = initial_pop[:,0]
m2s = initial_pop[:,1]
m3s = initial_pop[:,2]

a1s = initial_pop[:,3]
a2s = initial_pop[:,4]

e1s = initial_pop[:,5]
e2s = initial_pop[:,6]

INCL1s = initial_pop[:,7]
INCL2s = initial_pop[:,8]

AP1s = initial_pop[:,11]
AP2s = initial_pop[:,12]

LAN1s = initial_pop[:,9]
LAN2s = initial_pop[:,10]


f = h5py.File(filename, 'w')

f.create_dataset("m1s", data=m1s)
f.create_dataset("m2s", data=m2s)
f.create_dataset("m3s", data=m3s)

f.create_dataset("a1s", data=a1s)
f.create_dataset("e1s", data=e1s)
f.create_dataset("INCL1s", data=INCL1s)
f.create_dataset("AP1s", data=AP1s)
f.create_dataset("LAN1s", data=LAN1s)

f.create_dataset("a2s", data=a2s)
f.create_dataset("e2s", data=e2s)
f.create_dataset("INCL2s", data=INCL2s)
f.create_dataset("AP2s", data=AP2s)
f.create_dataset("LAN2s", data=LAN2s)
f.close()

print("Congrats! sample successfully saved in file",filename)


from __future__ import division
#Python 2.7.5 Using Brian
#Model from Ullah et al. 2009
from brian import *
#from numpy import *
#from scipy import *

#Params
N			= 100
dim			= 15
N1			= 10000000
beta		= 7.0        
n_k0		= 14.0		* mmole
diffusion	= 250.0		* umeter2/second
epsilon		= 0.07/0.3	* hertz
G_glia		= 5.0/0.3	* mmole/second
deltax		= 10.0		* umeter
v_b			= -50.0		* mvolt
gamma		= 0.0
gamma_tilde	= 0.4	                       
nu			= 5
C			= 1.0		* ufarad#/cmeter2
phi			= 3.0		* hertz
tau_e		= 4.0		* msecond
tau_i		= 8.0		* msecond
A			= 20.0		# maximum transmitter concentration
g_ca		= 0.1		* msiemens#/meter2 
v_ca		= 120.0		* mvolt
g_l			= 0.05		* msiemens#/meter2
El			= -65.0		* mvolt
g_na		= 100.0		* msiemens#/meter2
ENa		= 55.0		* mvolt
g_kd			= 40.0		* msiemens#/meter2
g_ahp		= 0.01		* msiemens#/meter2
EK			= -80.0		* mvolt
v_ee		= 0.0		* mvolt
v_ie		= -80.0		* mvolt
v_ei		= 0.0		* mvolt
v_ii		= -80.0		* mvolt
v_sp		= 40		* mvolt
alpha_ee	= 0.12		# for synaptic weights
alpha_ie	= 0.06		# for synaptic weights
alpha_ei	= 0.2 		# for synaptic weights
alpha_ii	= 0.02		# for synaptic weights
alpha_g	= 0.0

#Equations

eqs_e = Equations('''
dv/dt = (g_l*(El-v) + g_na*(m_inf*m_inf*m_inf)*h*(v-ENa)-g_kd*(n*n*n*n)*(v-EK))/C : volt
#dm/dt = alpha_m*(1-m)-beta_m*m : 1
dn/dt = alpha_n*(1-n)-beta_n*n : 1
dh/dt = alpha_h*(1-h)-beta_h*h : 1
# dge/dt = -ge*(1./taue) : siemens
# dgi/dt = -gi*(1./taui) : siemens
alpha_m = 0.1*(mV**-1)*(v+30*mV)/ \
    (1.0 - exp((v+30*mV)/(10*mV)))/ms : Hz
beta_m = 4.0*exp((v+55*mV)/(-18*mV))/ms : Hz
alpha_n = 0.01*(mV**-1)*(34*mV+v)/ \
    (1.0 - exp((v+34*mV)/(-10*mV)))/ms : Hz
beta_n = .125*exp((v+44.0*mV)/(-80*mV))/ms : Hz
alpha_h = 0.07*exp((v+44*mV)/(-20*mV))/ms : Hz
beta_h = 1./(1.0+exp((14*mV+v)/(-10*mV)))/ms : Hz
m_inf = alpha_m/(alpha_m + beta_m) : 1

#I_k	= ((g_kd/msiemens) * (n * n * n * n) + g_ahp/msiemens * Cai/(1.0+Cai) )*(v/mV-EK/mV)*uamp	: amp
''')

eqs_i = Equations('''
dv/dt = (g_l*(El-v) + g_na*(m_inf*m_inf*m_inf)*h*(v-ENa)-g_kd*(n*n*n*n)*(v-EK))/C : volt
#dm/dt = alpha_m*(1-m)-beta_m*m : 1
dn/dt = alpha_n*(1-n)-beta_n*n : 1
dh/dt = alpha_h*(1-h)-beta_h*h : 1
# dge/dt = -ge*(1./taue) : siemens
# dgi/dt = -gi*(1./taui) : siemens
alpha_m = 0.1*(mV**-1)*(v+30*mV)/ \
    (1.0 - exp((v+30*mV)/(10*mV)))/ms : Hz
# betam = 0.28*(mV**-1)*(v-VT-40*mV)/ \
# (exp((v-VT-40*mV)/(5*mV))-1)/ms : Hz
beta_m = 4.0*exp((v+55*mV)/(-18*mV))/ms : Hz
alpha_n = 0.01*(mV**-1)*(34*mV+v)/ \
    (1.0 - exp((v+34*mV)/(-10*mV)))/ms : Hz
beta_n = .125*exp((v+44.0*mV)/(-80*mV))/ms : Hz
alpha_h = 0.07*exp((v+44*mV)/(-20*mV))/ms : Hz
beta_h = 1./(1.0+exp((14*mV+v)/(-10*mV)))/ms : Hz
m_inf = alpha_m/(alpha_m + beta_m) : 1
''')

#eqs_s = '''
#I_syn = (-1.0*w*s_pre*chi_pre*(v_post/mV - v_ee/mV)/N)*uamp : amp
#w : 1.0
#'''

print 'Creating cells and synapses...'
#Create cells
myclock=Clock(dt=0.02*ms)
E = NeuronGroup(N=10, model=eqs_e, threshold=EmpiricalThreshold(state='v',threshold=-20*mV,refractory=2*ms),method='Euler',clock=myclock)
I = NeuronGroup(N=10, model=eqs_i, threshold=EmpiricalThreshold(state='v',threshold=-20*mV,refractory=2*ms),method='Euler',clock=myclock)

#Define synapses
Cee = Connection(E, E, 'v')
Cei = Connection(E, I, 'v')
Cii = Connection(I, I, 'v')
Cie = Connection(I, E, 'v')

#Create synapses
Cee.connect_random(sparseness=0.1, weight=0.5*mvolt)
Cei.connect_random(sparseness=0.1, weight=0.5*mvolt)
Cie.connect_random(sparseness=0.1, weight=0.5*mvolt)
Cii.connect_random(sparseness=0.1, weight=0.5*mvolt)
print 'done.'

#Set up external inputs
print 'setting up external input...'
spiketimes = [(1,20*ms),(1,23*ms),(1,25*ms)]
G = SpikeGeneratorGroup(1, spiketimes)
CGe = Connection(G,E)
CGi = Connection(G,I)
CGe.connect_full(weight=50*volt)
CGe.connect_full(weight=50*volt)
print 'done.'

#Set up synaptic weights
# for j in range(0,N): # From cell j in each group
	# for k in range(0,N): # To cell k in each group
		# if j == k: # No connection to itself
			# S_ee.w[j,k] = 0.0 #*nS
			# S_ei.w[j,k] = alpha_ei*sq30#*nS
			# S_ie.w[j,k] = alpha_ie*sq30#*nS
			# S_ii.w[j,k] = 0.0 #*nS
		# else:
			# distance = float((j-k)/N)**2.0
			# S_ee.w[j,k] = alpha_ee*sq100*exp(-100.0*distance)#*nS
			# S_ei.w[j,k] = alpha_ei*sq30*exp(-30.0*distance)#*nS
			# S_ie.w[j,k] = alpha_ie*sq30*exp(-30.0*distance)#*nS
			# S_ii.w[j,k] = alpha_ii*sq30*exp(-30.0*distance)#*nS
# print 'done.'

# S_ee.delay = 0 * ms
# S_ei.delay = 0 * ms
# S_ie.delay = 0 * ms
# S_ii.delay = 0 * ms

# TODO : Implement I_rand, I_ext, diffusion
# thresh	= 0.0

# Initialization
E.v = -65.1034452264476897 * mvolt
E.n = 0.645021617064742286*0.1         #potassium channel activation variable n for e network
E.h = 0.981306641820766101             #sodium channel inactivation variable h for e network
#E.s = 0.217441399671948571d-05         !variable s for temporal evolution of synaptic efficacy emanating from e network 
#E.Cai = 0.978772044450795857*0.0000001 *mole         #calcium concentration for e network

I.v = -65.1034452264476897 *mvolt             #membrane potential of inhibitory network
I.n = 0.645021617064742286 * 0.1         #potassium channel activation variable n for e network
I.h = 0.981306641820766101             #sodium channel inactivation variable h for inhibitory network
          # x((i-1)*dim+9) = 0.217441399671948571d-05         !variable s for temporal evolution of synaptic efficacy emanating from inhibitory network   
          # x((i-1)*dim+10) = 7.22                            !extracellular potassium concentration for e network
          # x((i-1)*dim+11) = 18.54                           !intracellular sodium concentration for e network
          # x((i-1)*dim+12) = 7.22                            !extracellular potassium concentration for inhibitory network
          # x((i-1)*dim+13) = 18.54                           !intracellular sodium concentration for inhibitory network
          # x((i-1)*dim+14) = 0.0								!variable eta, modeling the synaptic block due to the depolarization for e network
          # x((i-1)*dim+15) = 0.0								!variable eta, modeling the synaptic block due to the depolarization for inhibitory network 

M_e = SpikeMonitor(E)
M_i = SpikeMonitor(I)
Mv_e = StateMonitor(E, 'v', record=True)
Mv_i = StateMonitor(I, 'v', record=True)

try:
	print 'running simulation...'
	run(50 * ms)
except Exception as e:
	print str(type(e))
	print str(e.args)
	print e.message
	print str(e)

print M_e.nspikes
print 'done.'
se0 = [ M_e.spikes[i][0] for i in xrange(len(M_e.spikes)) ]
se1 = [ M_e.spikes[i][1] for i in xrange(len(M_e.spikes)) ]
scatter(se1,se0,color='r')
si0 = [ M_i.spikes[i][0] for i in xrange(len(M_i.spikes)) ]
si1 = [ M_i.spikes[i][1] for i in xrange(len(M_i.spikes)) ]
scatter(si1,si0,color='b')
plot(Mv_e[1])
plot(Mv_i[2],color='k')
show()
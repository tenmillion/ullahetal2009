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
A			= 20.0
g_ca		= 0.1		* msiemens#/meter2 
v_ca		= 120.0		* mvolt
g_l			= 0.05		* msiemens#/meter2
v_l			= -65.0		* mvolt
g_na		= 100.0		* msiemens#/meter2
v_na		= 55.0		* mvolt
g_k			= 40.0		* msiemens#/meter2
g_ahp		= 0.01		* msiemens#/meter2
v_k			= -80.0		* mvolt
v_ee		= 0.0		* mvolt
v_ie		= -80.0		* mvolt
v_ei		= 0.0		* mvolt
v_ii		= -80.0		* mvolt
v_sp		= 40		* mvolt
alpha_ee	= 0.12
alpha_ie	= 0.06
alpha_ei	= 0.2 
alpha_ii	= 0.02
alpha_g	= 0.0

#Equations
eqs_e = Equations('''
alpha_n_v = 0.01 * (v/mV + 34.0)/( 1.0 - exp(-0.1 * (v/mV + 34.0)) ) /second : Hz
alpha_h_v = 0.07 * exp(-1.0*(v/mV + 44.0)/20.0) /second : Hz
beta_n_v = 0.125 * exp(-1.0*(v/mV + 44.0)/80.0) /second : Hz
alpha_m_v = 0.1 * (v/mV+30.0)/( 1.0 - exp(-0.1 * (v/mV + 30.0)) ) /second : Hz
beta_m_v = 4.0 * exp(-1.0*(v/mV + 55.0)/18.0) /second : Hz
beta_h_v = 1.0 * 1.0 /( 1.0 + exp(-0.1 * (v/mV + 14.0)) ) /second : Hz
m_inf_v = alpha_m_v/(alpha_m_v + beta_m_v) : 1

I_k	= ((g_k/msiemens) * (n * n * n * n) + g_ahp/msiemens * Cai/(1.0+Cai) )*(v/mV-v_k/mV)*uamp	: amp
I_na	= (g_na/msiemens) * (m_inf_v * m_inf_v * m_inf_v) * h * (v/mV-v_na/mV)*uamp		: amp
I_mem	= -(g_l/msiemens) * (v/mV - v_l/mV)*uamp - I_na - I_k						: amp

dv/dt	= (1/C) * (I_mem)		: volt
# + I_ext + I_rand + I_syn
dn/dt	= phi * ( alpha_n_v*second * (1.0-n) - beta_n_v*second * n )	: 1.0
dh/dt	= phi * ( alpha_h_v*second * (1.0-h) - beta_h_v*second * h )	: 1.0
dCai/dt = (-0.002 * (g_ca/msiemens) * (v/mV - v_ca/mV)/( 1.0 + exp( -1.0*(v/mV+25.0)/2.5 ) ) - Cai/80.0 )/second : 1.0
''')

eqs_i = Equations('''
alpha_n_v = 0.01 * (v/mV + 34.0)/( 1.0 - exp(-0.1 * (v/mV + 34.0)) ) /second : Hz
alpha_h_v = 0.07 * exp(-1.0*(v/mV + 44.0)/20.0) /second : Hz
beta_n_v = 0.125 * exp(-1.0*(v/mV + 44.0)/80.0) /second : Hz
alpha_m_v = 0.1 * (v/mV+30.0)/( 1.0 - exp(-0.1 * (v/mV + 30.0)) ) /second : Hz
beta_m_v = 4.0 * exp(-1.0*(v/mV + 55.0)/18.0) /second : Hz
beta_h_v = 1.0 * 1.0 /( 1.0 + exp(-0.1 * (v/mV + 14.0)) ) /second : Hz
m_inf_v = alpha_m_v/(alpha_m_v + beta_m_v) : 1

I_k		= (g_k/msiemens) * (n * n * n * n)*(v/mV-v_k/mV)*uamp	: amp
I_na	= (g_na/msiemens) * (m_inf_v * m_inf_v * m_inf_v) * h * (v/mV-v_na/mV)*uamp	: amp
I_mem	= -(g_l/msiemens) * (v/mV - v_l/mV)*uamp - I_na - I_k						: amp

dv/dt	= (1/C) * I_mem		: volt
# + I_ext + I_rand
dn/dt	= phi * ( alpha_n_v*second * (1.0-n) - beta_n_v*second * n )	: 1.0
dh/dt	= phi * ( alpha_h_v*second * (1.0-h) - beta_h_v*second * h )	: 1.0
''')

#eqs_s = '''
#I_syn = (-1.0*w*s_pre*chi_pre*(v_post/mV - v_ee/mV)/N)*uamp : amp
#w : 1.0
#'''

print 'Creating cells and synapses...'
#Create cells
Excitatory = NeuronGroup(N=100, model=eqs_e, threshold=EmpiricalThreshold(state='v',threshold=-20*mV,refractory=1*ms))
Inhibitory = NeuronGroup(N=100, model=eqs_i, threshold=EmpiricalThreshold(state='v',threshold=-20*mV,refractory=1*ms))

#Define synapses
#S_ee = Synapses(Excitatory, Excitatory, model=eqs_s, pre='I_mem_post+=I_syn')
#S_ei = Synapses(Excitatory, Inhibitory, model=eqs_s, pre='I_mem_post+=I_syn')
#S_ie = Synapses(Inhibitory, Excitatory, model=eqs_s, pre='I_mem_post+=I_syn')
#S_ii = Synapses(Inhibitory, Inhibitory, model=eqs_s, pre='I_mem_post+=I_syn')

Cee = Connection(Excitatory, Excitatory, 'v', sparseness=0.6, weight=(60 * 0.27 / 10) * mV)
Cei = Connection(Excitatory, Inhibitory, 'v', sparseness=0.6, weight=(60 * 0.27 / 10) * mV)
Cii = Connection(Inhibitory, Inhibitory, 'v', sparseness=0.6, weight=(30 * 0.27 / 10) * mV)
Cie = Connection(Inhibitory, Excitatory, 'v', sparseness=0.6, weight=(30 * 0.27 / 10) * mV)

#Create synapses
# S_ee[:,:] = True
# S_ei[:,:] = True
# S_ie[:,:] = True
# S_ii[:,:] = True
Cee[:,:] = True
Cei[:,:] = True
Cie[:,:] = True
Cii[:,:] = True
print 'done.'

#Set up synaptic weights
# print 'Setting up synaptic weights...'
# sq100 = sqrt(100.0/pi)
# sq30 = sqrt(30.0/pi)

# print 'done.'

# spiketimes = [(0, 1 * ms), (0, 4 * ms),
              # (1, 2 * ms), (1, 3 * ms)]
# G1 = SpikeGeneratorGroup(2, spiketimes)

# G = Connection(G1,Excitatory,sparse

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

M_e = SpikeMonitor(Excitatory)
M_i = SpikeMonitor(Inhibitory)
# Ma = StateMonitor(Excitatory, 'v', record=True)
# Mb = StateMonitor(Inhibitory, 'v', record=True)

try:
	print 'running simulation...'
	run(500 * ms)
except Exception as e:
	print str(type(e))
	print str(e.args)
	print e.message
	print str(e)

print 'done.'
se0 = [ M_e.spikes[i][0] for i in xrange(len(M_e.spikes)) ]
se1 = [ M_e.spikes[i][1] for i in xrange(len(M_e.spikes)) ]
scatter(se1,se0,color='r')
si0 = [ M_i.spikes[i][0] for i in xrange(len(M_i.spikes)) ]
si1 = [ M_i.spikes[i][1] for i in xrange(len(M_i.spikes)) ]
scatter(si1,si0,color='b')
show()
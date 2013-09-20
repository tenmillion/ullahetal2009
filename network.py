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
v_l_init	= -65.0		* mvolt
g_na		= 100.0		* msiemens#/meter2
v_na_init	= 55.0		* mvolt
g_k			= 40.0		* msiemens#/meter2
g_ahp		= 0.01		* msiemens#/meter2
v_k_init	= -80.0		* mvolt
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
sigma_v = 1.0/( 1.0 + exp(-1.0*(v/mV + 20.0)/(4.0)) ) /second : Hz
m_inf_v = alpha_m_v/(alpha_m_v + beta_m_v) : 1
I_pump = (1.25/(1.0 + exp((25.0 - Nai)/3.0)))*(1.0/(1.0 + exp(8.0 - Ko))) * uamp : amp
I_diff = ( epsilon*second*(Ko-n_k0/mmole) + (G_glia*second/mmole)/(1.0 + exp((18.0-Ko)/2.5)) ) * uamp : amp
Kin	= (140.0 +(18.0-Nai) ) : 1.0
Nao	= (144.0 - beta * (Nai - 18.0) ) : 1.0

v_k	= 26.64 * log(Ko/Kin) * mV 		: volt
v_na	= 26.64 * log(Nao/Nai) * mV		: volt
v_l	= 26.64 * log ((Ko + 0.065*Nao + 0.6*6.0) / (Kin + 0.065*Nai + 0.6*130.0)) * mV	: volt

I_k	= ((g_k/msiemens) * (n * n * n * n) + g_ahp/msiemens * Cai/(1.0+Cai) )*(v/mV-v_k/mV)*uamp	: amp
I_na	= (g_na/msiemens) * (m_inf_v * m_inf_v * m_inf_v) * h * (v/mV-v_na/mV)*uamp		: amp
I_mem	= -(g_l/msiemens) * (v/mV - v_l/mV)*uamp - I_na - I_k						: amp

dv/dt	= (1/C) * (I_mem)		: volt
# + I_ext + I_rand + I_syn
dn/dt	= phi * ( alpha_n_v*second * (1.0-n) - beta_n_v*second * n )	: 1.0
dh/dt	= phi * ( alpha_h_v*second * (1.0-h) - beta_h_v*second * h )	: 1.0
ds/dt	= (1/tau_e) * (A * sigma_v*second * (1.0 - s) - s)	: 1.0
dCai/dt = (-0.002 * (g_ca/msiemens) * (v/mV - v_ca/mV)/( 1.0 + exp( -1.0*(v/mV+25.0)/2.5 ) ) - Cai/80.0 )/second : 1.0

dKo/dt	= 0.001*((I_k/uamp)/3.0 - 7.0* 2.0*(I_pump/uamp) - I_diff/uamp)/second	: 1.0
dNai/dt	= 0.001*(-1.0*(I_na/uamp)/(7.0*3.0) - 3.0*I_pump/uamp)/second			: 1.0

deta/dt	= (gamma * (v/mV - v_b/mV) - gamma_tilde * eta*second)/(second**2)	: Hz

chi = exp(-1.0*eta*second/nu) : 1.0
''')

eqs_i = Equations('''
alpha_n_v = 0.01 * (v/mV + 34.0)/( 1.0 - exp(-0.1 * (v/mV + 34.0)) ) /second : Hz
alpha_h_v = 0.07 * exp(-1.0*(v/mV + 44.0)/20.0) /second : Hz
beta_n_v = 0.125 * exp(-1.0*(v/mV + 44.0)/80.0) /second : Hz
alpha_m_v = 0.1 * (v/mV+30.0)/( 1.0 - exp(-0.1 * (v/mV + 30.0)) ) /second : Hz
beta_m_v = 4.0 * exp(-1.0*(v/mV + 55.0)/18.0) /second : Hz
beta_h_v = 1.0 * 1.0 /( 1.0 + exp(-0.1 * (v/mV + 14.0)) ) /second : Hz
sigma_v = 1.0/( 1.0 + exp(-1.0*(v/mV + 20.0)/(4.0)) ) /second : Hz
m_inf_v = alpha_m_v/(alpha_m_v + beta_m_v) : 1
I_pump = (1.25/(1.0 + exp((25.0 - Nai)/3.0)))*(1.0/(1.0 + exp(8.0 - Ko))) * uamp : amp
I_diff = ( epsilon*second*(Ko-n_k0/mmole) + (G_glia*second/mmole)/(1.0 + exp((18.0-Ko)/2.5)) ) * uamp : amp
Kin	= (140.0 +(18.0-Nai) ) : 1.0
Nao	= (144.0 - beta * (Nai - 18.0) ) : 1.0

v_k	= 26.64 * log(Ko/Kin) * mV 		: volt
v_na	= 26.64 * log(Nao/Nai) * mV		: volt
v_l	= 26.64 * log ((Ko + 0.065*Nao + 0.6*6.0) / (Kin + 0.065*Nai + 0.6*130.0)) * mV	: volt

I_k		= (g_k/msiemens) * (n * n * n * n)*(v/mV-v_k/mV)*uamp	: amp
I_na	= (g_na/msiemens) * (m_inf_v * m_inf_v * m_inf_v) * h * (v/mV-v_na/mV)*uamp	: amp
I_mem	= -(g_l/msiemens) * (v/mV - v_l/mV)*uamp - I_na - I_k						: amp

dv/dt	= (1/C) * I_mem		: volt
# + I_ext + I_rand
dn/dt	= phi * ( alpha_n_v*second * (1.0-n) - beta_n_v*second * n )	: 1.0
dh/dt	= phi * ( alpha_h_v*second * (1.0-h) - beta_h_v*second * h )	: 1.0
ds/dt	= (1/tau_i) * (A * sigma_v*second * (1.0 - s) - s)	: 1.0

dKo/dt	= 0.001*((I_k/uamp)/3.0 - 7.0* 2.0*(I_pump/uamp) - I_diff/uamp)/second	: 1.0
dNai/dt	= 0.001*(-1.0*(I_na/uamp)/(7.0*3.0) - 3.0*I_pump/uamp)/second			: 1.0

deta/dt	= (gamma * (v/mV - v_b/mV) - gamma_tilde * eta*second)/(second**2)	: Hz
chi = exp(-1.0*eta*second/nu) : 1.0
''')

eqs_s = '''
I_syn = (-1.0*w*s_pre*chi_pre*(v_post/mV - v_ee/mV)/N)*uamp : amp
w : 1.0
'''

print 'Creating cells and synapses...'
#Create cells
Excitatory = NeuronGroup(N=100, model=eqs_e, threshold=EmpiricalThreshold(state='v',threshold=-20*mV,refractory=1*ms))
Inhibitory = NeuronGroup(N=100, model=eqs_i, threshold=EmpiricalThreshold(state='v',threshold=-20*mV,refractory=1*ms))

#Define synapses
S_ee = Synapses(Excitatory, Excitatory, model=eqs_s, pre='I_mem_post+=I_syn')
S_ei = Synapses(Excitatory, Inhibitory, model=eqs_s, pre='I_mem_post+=I_syn')
S_ie = Synapses(Inhibitory, Excitatory, model=eqs_s, pre='I_mem_post+=I_syn')
S_ii = Synapses(Inhibitory, Inhibitory, model=eqs_s, pre='I_mem_post+=I_syn')

#Create synapses
S_ee[:,:] = True
S_ei[:,:] = True
S_ie[:,:] = True
S_ii[:,:] = True
print 'done.'

#Set up synaptic weights
print 'Setting up synaptic weights...'
sq100 = sqrt(100.0/pi)
sq30 = sqrt(30.0/pi)

for j in range(0,N): # From cell j in each group
	for k in range(0,N): # To cell k in each group
		if j == k: # No connection to itself
			S_ee.w[j,k] = 0.0 #*nS
			S_ei.w[j,k] = alpha_ei*sq30#*nS
			S_ie.w[j,k] = alpha_ie*sq30#*nS
			S_ii.w[j,k] = 0.0 #*nS
		else:
			distance = float((j-k)/N)**2.0
			S_ee.w[j,k] = alpha_ee*sq100*exp(-100.0*distance)#*nS
			S_ei.w[j,k] = alpha_ei*sq30*exp(-30.0*distance)#*nS
			S_ie.w[j,k] = alpha_ie*sq30*exp(-30.0*distance)#*nS
			S_ii.w[j,k] = alpha_ii*sq30*exp(-30.0*distance)#*nS
print 'done.'

S_ee.delay = 0 * ms
S_ei.delay = 0 * ms
S_ie.delay = 0 * ms
S_ii.delay = 0 * ms

# TODO : Implement I_rand, I_ext, diffusion
# thresh	= 0.0

Monitor_e = SpikeMonitor(Excitatory)
Monitor_i = SpikeMonitor(Inhibitory)

try:
	print 'running simulation...'
	run(1 * second)
except Exception as e:
	print str(type(e))
	print str(e.args)
	print e.message
	print str(e)

print 'done.'
print Monitor_e.nspikes

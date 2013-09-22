from brian import *

# Parameters
area = 20000 * umetre ** 2
Cm = (1 * ufarad * cm ** -2) * area
gl = (5e-5 * siemens * cm ** -2) * area
El = -65 * mV #-65
EK = -80 * mV #-80
ENa = 55 * mV #55
g_na = (100 * msiemens * cm ** -2) * area
g_kd = (30 * msiemens * cm ** -2) * area
VTm = 55 * mV #55
VTh = 44 * mV #44
VTn = 44 * mV #44
# Time constants
taue = 5 * ms # Ah this is ok because dse/dt is in s**-1!
taui = 10 * ms
A = 20 * ms # dunno, but using ms to get s unitless
# Reversal potentials
Ee = 0 * mV
Ei = -80 * mV
# we = 6 * nS # excitatory synaptic weight (voltage)
# wi = 67 * nS # inhibitory synaptic weight
alpha_ee	= 0.12
alpha_ie	= 0.06
alpha_ei	= 0.2 
alpha_ii	= 0.02

# The model

eqs = Equations('''
dv/dt = (gl*(El-v)+Iesyn+Iisyn-\
    g_na*(m*m*m)*h*(v-ENa)-\
    g_kd*(n*n*n*n)*(v-EK))/Cm : volt
dm/dt = alpham*(1-m)-betam*m : 1
dn/dt = alphan*(1-n)-betan*n : 1
dh/dt = alphah*(1-h)-betah*h : 1
dse/dt = (A*sigma*(1-se)-se)/taue : 1
dsi/dt = (A*sigma*(1-si)-si)/taui : 1

alpham = 0.1*(mV**-1)*(25*mV-(v+VTm))/ \
    (exp((25*mV-(v+VTm))/(10*mV))-1.)/ms : Hz
betam = 4.0*exp((-(v+VTm))/(18*mV))/ms : Hz
alphah = 0.07*exp((-(v+VTh))/(20*mV))/ms : Hz
betah = 1./(1+exp((30*mV-(v+VTh))/(10*mV)))/ms : Hz
alphan = 0.01*(mV**-1)*(10*mV-(v+VTn))/ \
    (exp((10*mV-(v+VTn))/(10*mV))-1.)/ms : Hz
betan = .125*exp((-(v+VTn))/(80*mV))/ms : Hz
sigma = 1./(1+exp((-(v+20*mV))/(4*mV)))/ms : Hz

Iesyn : amp
Iisyn : amp
''')

eqs_esyn = '''
g_jk : nS	# synaptic weight
Iesyn = (g_jk * se_pre)*(Ee-v_post) : amp
'''

eqs_isyn = '''
g_jk : nS	# synaptic weight
Iisyn = (g_jk * si_pre)*(Ei-v_post) : amp
'''
N=400
P = NeuronGroup(N, model=eqs,
    threshold=EmpiricalThreshold(threshold= -20 * mV,
                                 refractory=3 * ms),
    implicit=True, freeze=True)
Pe = P.subgroup(N/5*4)
Pi = P.subgroup(N/5)

Se = Synapses(Pe, P, model=eqs_esyn)
Si = Synapses(Pi, P, model=eqs_isyn)
Se[:,:] = True
Si[:,:] = True
P.Iesyn = Se.Iesyn
P.Iisyn = Si.Iisyn

# Initialization
P.v = El + (randn(len(P)) * 3 - 3) * mV
# Se.g_jk = abs(randn(len(Se))+6) * nS   # Were these weights too big?
# Si.g_jk = abs(randn(len(Si))+67) * nS

#Set up synaptic weights
print "Setting up synaptic weights"
sq100 = sqrt(100.0/pi)
sq30 = sqrt(30.0/pi)
for j in range(0,N): # From cell j in each group
	for k in range(0,N): # To cell k in each group
		if j == k: # No connection to itself
			Se.g_jk[j,k] = 0.0 *nS
			Se.g_jk[j,k] = alpha_ei*sq30*nS
			Si.g_jk[j,k] = alpha_ie*sq30*nS
			Si.g_jk[j,k] = 0.0 *nS
		else:
			distance = float((j-k)/N)**2.0
			Se.g_jk[j,k] = alpha_ee*sq100*exp(-100.0*distance)*nS
			Si.g_jk[j,k] = alpha_ii*sq30*exp(-30.0*distance)*nS
print 'done.'

# External input
spiketimes = [(0,100*ms),(0,100*ms)]
G = SpikeGeneratorGroup(1, spiketimes)
Pe_in = Pe.subgroup(N/5)
Input = Connection(G,Pe_in,weight=5*mV,sparseness=0.3)

# Record the number of spikes and a few traces
trace = StateMonitor(Pe, 'v', record=arange(0,40))
trace2 = StateMonitor(Pi, 'v', record=arange(0,40))
Me = SpikeMonitor(Pe)
Mi = SpikeMonitor(Pi)

run(500 * msecond)
print Me.nspikes
print Mi.nspikes

subplot(411)
raster_plot(Me)
xlim(0,500)
subplot(412)
raster_plot(Mi)
xlim(0,500)
subplot(413)
for i in arange(0,40):
	plot(trace[i])
#xlim(len(trace[1])*99.5/500,len(trace[1])*115/500)
subplot(414)
for i in arange(0,40):
	plot(trace2[i])
#xlim(len(trace[1])*99.5/500,len(trace[1])*115/500)
savefig('foo.png')

# f = open('spikes.txt', 'w')
# for i in arange(0,400):
	# f.write(M[i])
# f.close()
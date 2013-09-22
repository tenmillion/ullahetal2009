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
taue = 5 * ms
taui = 10 * ms
A = 20 * ms # dunno, but using ms to get s unitless
# Reversal potentials
Ee = 0 * mV
Ei = -80 * mV
we = 6 * nS # excitatory synaptic weight (voltage)
wi = 67 * nS # inhibitory synaptic weight

# The model

eqs = Equations('''
dv/dt = (gl*(El-v)+Iesyn+Iisyn-\
    g_na*(m*m*m)*h*(v-ENa)-\
    g_kd*(n*n*n*n)*(v-EK))/Cm : volt
dm/dt = alpham*(1-m)-betam*m : 1
dn/dt = alphan*(1-n)-betan*n : 1
dh/dt = alphah*(1-h)-betah*h : 1
dse/dt = ((A*ms)*sigma*(1-se)-se)/taue : 1
dsi/dt = ((A*ms)*sigma*(1-si)-si)/taui : 1

alpham = 0.1*(mV**-1)*(25*mV-(v+VTm))/ \
    (exp((25*mV-(v+VTm))/(10*mV))-1.)/ms : Hz
betam = 4.0*exp((-(v+VTm))/(18*mV))/ms : Hz
alphah = 0.07*exp((-(v+VTh))/(20*mV))/ms : Hz
betah = 1./(1+exp((30*mV-(v+VTh))/(10*mV)))/ms : Hz
alphan = 0.01*(mV**-1)*(10*mV-(v+VTn))/ \
    (exp((10*mV-(v+VTn))/(10*mV))-1.)/ms : Hz
betan = .125*exp((-(v+VTn))/(80*mV))/ms : Hz
sigma = 1./(1+exp((-(v+20*mV))/(4*mV)))/ms : Hz

Iesyn : A
Iisyn : A
''')

eqs_esyn = '''
g_jk : 1	# synaptic weight
Iesyn = (g_jk * se_pre)*(Ee-v_post) : A
'''

eqs_isyn = '''
g_jk : 1	# synaptic weight
Iisyn = (g_jk * si_pre)*(Ei-v_post) : A
'''

P = NeuronGroup(4000, model=eqs,
    threshold=EmpiricalThreshold(threshold= -20 * mV,
                                 refractory=3 * ms),
    implicit=True, freeze=True)
Pe = P.subgroup(3200)
Pi = P.subgroup(800)

Se = Synapses(Pe, P, model=eqs_esyn)
Si = Synapses(Pi, P, model=eqs_isyn)
Se[:,:] = True
P.Iesyn = Se.Iesyn
P.Iisyn = Si.Iisyn

# Initialization
P.v = El + (randn(len(P)) * 5 - 5) * mV
Se.g_jk = (randn(len(P)) * 1.5 + 4) * 10. * nS
Si.g_jk = (randn(len(P)) * 12 + 20) * 10. * nS

# External input
spiketimes = [(0,100*ms)]
G = SpikeGeneratorGroup(1, spiketimes)
Input = Connection(G,Pe,weight=30*mV,sparseness=0.5)

# Record the number of spikes and a few traces
trace = StateMonitor(P, 'v', record=arange(0,40))
trace2 = StateMonitor(P, 'gi', record=arange(0,40))

M = SpikeMonitor(P)
run(500 * msecond)
print M.nspikes
subplot(311)
raster_plot(M)
subplot(312)
for i in arange(0,40):
	plot(trace[i])
subplot(313)
for i in arange(0,40):
	plot(trace2[i])
show()
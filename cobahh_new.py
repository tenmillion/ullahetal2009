from brian import *

# Parameters
area = 20000 * umetre ** 2
Cm = (1 * ufarad * cm ** -2) * area
gl = (5e-5 * siemens * cm ** -2) * area
El = -60 * mV
EK = -90 * mV
ENa = 50 * mV
g_na = (100 * msiemens * cm ** -2) * area
g_kd = (30 * msiemens * cm ** -2) * area
VT = -50 * mV
# Time constants
taue = 5 * ms
taui = 10 * ms
# Reversal potentials
Ee = 0 * mV
Ei = -80 * mV
we = 6 * nS # excitatory synaptic weight (voltage)
wi = 67 * nS # inhibitory synaptic weight

# The model
# +ge*(Ee-v)+gi*(Ei-v)
eqs = Equations('''
# dv/dt = (gl*(El-v)-\
    # g_na*(m*m*m)*h*(v-ENa)-\
    # g_kd*(n*n*n*n)*(v-EK))/Cm : volt
# dm/dt = alpham*(1-m)-betam*m : 1
# dn/dt = alphan*(1-n)-betan*n : 1
# dh/dt = alphah*(1-h)-betah*h : 1
# dge/dt = -ge*(1./taue) : siemens
# dgi/dt = -gi*(1./taui) : siemens
# alpham = 0.1*(mV**-1)*(25*mV-v+VT)/ \
    # (exp((25*mV-v+VT)/(10*mV))-1.)/ms : Hz
# betam = 4.0*exp((-v+VT)/(18*mV))/ms : Hz
# alphah = 0.07*exp((-v+VT)/(20*mV))/ms : Hz
# betah = 1./(1+exp((30*mV-v+VT)/(10*mV)))/ms : Hz
# alphan = 0.01*(mV**-1)*(10*mV-v+VT)/ \
    # (exp((10*mV-v+VT)/(10*mV))-1.)/ms : Hz
# betan = .125*exp((-v+VT)/(80*mV))/ms : Hz
# ''')



P = NeuronGroup(4000, model=eqs,
    threshold=EmpiricalThreshold(threshold= -20 * mV,
                                 refractory=3 * ms),
    implicit=True, freeze=True)
Pe = P.subgroup(3200)
Pi = P.subgroup(800)
# Ce = Connection(Pe, P, 'ge', weight=we, sparseness=0.2)
# Ci = Connection(Pi, P, 'gi', weight=wi, sparseness=0.02)
# Initialization
P.v = El + (randn(len(P)) * 5 - 5) * mV
# P.ge = (randn(len(P)) * 1.5 + 4) * 10. * nS
# P.gi = (randn(len(P)) * 12 + 20) * 10. * nS

# External input
spiketimes = [(0,100*ms)]
G = SpikeGeneratorGroup(1, spiketimes)
Input = Connection(G,Pe,weight=30*mV,sparseness=0.5)

# Record the number of spikes and a few traces
trace = StateMonitor(P, 'v', record=arange(37,50))

M = SpikeMonitor(P)
run(500 * msecond)
print M.nspikes
subplot(211)
raster_plot(M)
subplot(212)
for i in arange(37,50):
	plot(trace[i])
show()
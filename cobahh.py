from brian import *

# Parameters- Gutkin
# area = 20000 * umetre ** 2
area = 1 * cm ** 2
Cm = (1 * ufarad * cm ** -2) * area
gl = (5e-5 * siemens * cm ** -2) * area
El = -65 * mV
EK = -80 * mV
ENa = 100 * mV
g_na = (100 * msiemens * cm ** -2) * area
g_kd = (40 * msiemens * cm ** -2) * area
VT = 0 * mV
# Time constants
taue = 5 * ms
taui = 10 * ms
# Reversal potentials
Ee = 0 * mV
Ei = -80 * mV
we = 6 * nS # excitatory synaptic weight (voltage)
wi = 67 * nS # inhibitory synaptic weight
phi = 3.0
g_ca = (0.1 * msiemens * cm ** -2) * area
ECa = 120 * mV

# Parameters - Original
# area = 20000 * umetre ** 2
# Cm = (1 * ufarad * cm ** -2) * area
# gl = (5e-5 * siemens * cm ** -2) * area
# El = -60 * mV
# EK = -90 * mV
# ENa = 50 * mV
# g_na = (100 * msiemens * cm ** -2) * area
# g_kd = (30 * msiemens * cm ** -2) * area
# VT = -63 * mV
# Time constants
# taue = 5 * ms
# taui = 10 * ms
# Reversal potentials
# Ee = 0 * mV
# Ei = -80 * mV
# we = 6 * nS # excitatory synaptic weight (voltage)
# wi = 67 * nS # inhibitory synaptic weight


# The model - Gutkin
# Not using m_inf
# dv/dt = (gl*(El-v)-\
    # g_na*(m*m*m)*h*(v-ENa)-\
    # g_kd*(n*n*n*n)*(v-EK))/Cm : volt
# dm/dt = phi*(alpham*(1-m)-betam*m) : 1
# dn/dt = phi*(alphan*(1-n)-betan*n) : 1
# dh/dt = phi*(alphah*(1-h)-betah*h) : 1
# dge/dt = -ge*(1./taue) : siemens
# dgi/dt = -gi*(1./taui) : siemens
eqs = Equations('''
dv/dt = (gl*(El-v)-\
    g_na*(m_inf*m_inf*m_inf)*h*(v-ENa)-\
    (g_kd*n*n*n*n+(g_ahp*Cai/(1+Cai)))*(v-EK))/Cm : volt
#dm/dt = phi*(alpham*(1-m)-betam*m) : 1
dn/dt = phi*(alphan*(1-n)-betan*n) : 1
dh/dt = phi*(alphah*(1-h)-betah*h) : 1
dCai/dt = (-0.002 * (g_ca/msiemens) * (v/mV - ECa/mV)/( 1.0 + exp( -1.0*(v/mV+25.0)/2.5 ) ) - Cai/80.0 )/ms : 1.0

#dge/dt = -ge*(1./taue) : siemens
#dgi/dt = -gi*(1./taui) : siemens

alpham = 0.1*(mV**-1)*(-30*mV-v+VT)/ \
    (exp((-30*mV-v+VT)/(10*mV))-1.)/ms : Hz
betam = 4.0*exp((-55*mV-v+VT)/(18*mV))/ms : Hz
alphah = 0.07*exp((44*mV-v+VT)/(20*mV))/ms : Hz
betah = 1./(1+exp((-14*mV-v+VT)/(10*mV)))/ms : Hz
alphan = 0.01*(mV**-1)*(-34*mV-v+VT)/ \
    (exp((-34*mV-v+VT)/(10*mV))-1.)/ms : Hz
betan = .125*exp((-44*mV-v+VT)/(80*mV))/ms : Hz

m_inf = alpham /(alpham+betam) : 1
''')

# The model - Original
# eqs = Equations('''
# dv/dt = (gl*(El-v)+ge*(Ee-v)+gi*(Ei-v)-\
    # g_na*(m*m*m)*h*(v-ENa)-\
    # g_kd*(n*n*n*n)*(v-EK))/Cm : volt
# dm/dt = alpham*(1-m)-betam*m : 1
# dn/dt = alphan*(1-n)-betan*n : 1
# dh/dt = alphah*(1-h)-betah*h : 1
# dge/dt = -ge*(1./taue) : siemens
# dgi/dt = -gi*(1./taui) : siemens
# alpham = 0.32*(mV**-1)*(13*mV-v+VT)/ \
    # (exp((13*mV-v+VT)/(4*mV))-1.)/ms : Hz
# betam = 0.28*(mV**-1)*(v-VT-40*mV)/ \
    # (exp((v-VT-40*mV)/(5*mV))-1)/ms : Hz
# alphah = 0.128*exp((17*mV-v+VT)/(18*mV))/ms : Hz
# betah = 4./(1+exp((40*mV-v+VT)/(5*mV)))/ms : Hz
# alphan = 0.032*(mV**-1)*(15*mV-v+VT)/ \
    # (exp((15*mV-v+VT)/(5*mV))-1.)/ms : Hz
# betan = .5*exp((10*mV-v+VT)/(40*mV))/ms : Hz
# ''')

myclock=Clock(dt=0.02*ms)
P = NeuronGroup(40, model=eqs,
    threshold=EmpiricalThreshold(threshold= -20 * mV,
                                 refractory= 3 * ms),
    implicit=False,order=2, freeze=True, clock=myclock)
    # implicit=True, freeze=True, clock=myclock)
Pe = P.subgroup(32)
Pi = P.subgroup(8)
#Ce = Connection(Pe, P, 'ge', weight=we, sparseness=0.2,seed=20548)
#Ci = Connection(Pi, P, 'gi', weight=wi, sparseness=0.2,seed=23234)

# print Ce[:,1]

# Initialization
P.v = El + (randn(len(P)) * 5 - 5) * mV
#P.ge = (randn(len(P)) * 1.5 + 4) * 10. * nS
#P.gi = (randn(len(P)) * 12 + 20) * 10. * nS

# Record the number of spikes and a few traces
trace = StateMonitor(P, 'v', record=arange(0,40))
M = SpikeMonitor(P)
run(.5 * second)
print M.nspikes
subplot(211)
raster_plot(M)
subplot(212)
for i in arange(0,40):
	plot(trace[i])
show()
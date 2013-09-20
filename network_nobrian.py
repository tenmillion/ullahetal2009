from __future__ import division
#Python 2.7.5 Using Brian
#Model from Ullah et al. 2009
from brian import *
from math import *
from numpy import *
#from scipy import *

#Params
N			= 100
dim			= 15
N1			= 10000000
beta		= 7.0        
n_k0		= 14.0		# mmole
diffusion	= 250.0		# umeter2/second
epsilon		= 0.07/0.3	# hertz
G_glia		= 5.0/0.3	# mmole/second
deltax		= 10.0		# umeter
v_b			= -50.0		# mvolt
gamma_tilde	= 0.4	                       
nu			= 5
C			= 1.0		# ufarad/cmeter2
phi			= 3.0		# hertz
tau_e		= 4.0		# msecond
A			= 20.0
g_ca		= 0.1		# msiemens/meter2 
V_ca		= 120.0		# mvolt
g_l			= 0.05		# msiemens/meter2
V_l			= -65.0		# mvolt
g_na		= 100.0		# msiemens/meter2
V_na		= 55.0		# mvolt
g_k			= 40.0		# msiemens/meter2
g_ahp		= 0.01		# msiemens/meter2
V_k			= -80.0		# mvolt
tau_i		= 8.0		# msecond
V_ee		= 0.0		# mvolt
V_ie		= -80.0		# mvolt
V_ei		= 0.0		# mvolt
V_ii		= -80.0		# mvolt
V_sp		= 40		# mvolt

def fullmodel( t, y ):
	"The Full model from Cressman et al. 2009"
	V=y[0] #V, the membrane voltage
	n=y[1] #n, gating variable
	h=y[2] #h, gating variable
	Ko=y[3] #[K]_o, the extracellular potassium concentration
	Nai=y[4] #[Na]_i, the intracellular sodium concentration
	Cai=y[5] #[Ca]_i, the calcium concentration

# Gating variables
	if V != -34:
		alpha_n = 0.01 * (V+34.0)/( 1.0 - np.exp(-0.1 * (V+34.0)) )
	else:
		alpha_n = 0.1

	beta_n = 0.125 * np.exp(-(V+44.0)/80.0)

	if V != -30:
		alpha_m = 0.1 * (V+30.0)/( 1.0 - np.exp(-0.1 * (V+30.0)) )
	else:
		alpha_m = 1

	beta_m = 4.0 * np.exp(-(V+55.0)/18.0)
	alpha_h = 0.07 * np.exp(-(V+44.0)/20.0)
	beta_h = 1.0/( 1.0 + np.exp(-0.1 * (V+14.0)) )

	m_inf = alpha_m/(alpha_m + beta_m)

# Ion concentrations
	Ki = (158.0-Nai)
	Nao = (144.0-beta*(Nai-18.0))

# Reversal potentials
	E_k = 26.64 * np.log((Ko/Ki))
	E_na = 26.64 * np.log((Nao/Nai))

# Membrane currents
	Ina = g_na*(m_inf*m_inf*m_inf)*h*(V-E_na) + g_naL*(V-E_na)
	Ik = (g_k*n*n*n*n + g_ahp*Cai/(1+Cai))*(V-E_k) + g_kL*(V-E_k)
	Icl = g_clL*(V-E_cl)

# Extracellular currents
	Ipump = (rho/(1.0+np.exp((25.0-Nai)/3.0)))*(1/(1+np.exp(5.5-Ko)))
	Iglia = (glia/(1.0+np.exp((18.0-Ko)/2.5)))
	Idiff = epsilon*(Ko-kbath)

# Derivatives
	dVdt = (1.0/Cm)*(-Ina -Ik -Icl)
	dndt = phi*(alpha_n*(1-n)-beta_n*n)
	dhdt = phi*(alpha_h*(1-h)-beta_h*h)
	dKodt = (1/tau)*(gam*beta*Ik -2.0*beta*Ipump -Iglia -Idiff)
	dNaidt = (1/tau)*(-gam*Ina -3.0*Ipump)
	dCaidt = -0.002*g_ca*(V-E_ca)/(1+np.exp(-(V+25.0)/2.5)) -Cai/80.0

	return np.array([dVdt, dndt, dhdt, dKodt, dNaidt, dCaidt], dtype=np.float64)

#Equations
#excitatory
alpha_n_Ve = 0.01 * (Ve+34.0)/( 1.0 - exp(-0.1 * (Ve + 34.0)) ) # volt
alpha_h_Ve = 0.07 * exp(-1.0*(Ve + 44.0)/20.0) # volt
beta_n_Ve = 0.125 * exp(-1.0*(Ve + 44.0)/80.0) # volt
alpha_m_Ve = 0.1 * (Ve+30.0)/( 1.0 - exp(-0.1 * (Ve + 30.0)) ) # volt
beta_m_Ve = 4.0 * exp(-1.0*(Ve + 55.0)/18.0) # volt
beta_h_Ve = 1.0/( 1.0 + exp(-0.1 * (Ve + 14.0)) ) # volt
sigma_Ve = 1.0/( 1.0 + exp(-1.0*(Ve + 20.0)/4.0) ) # volt
m_inf_Ve = alpha_m_Ve/(alpha_m_Ve + beta_m_Ve) # volt
I_pump_e = (1.25/(1.0 + exp((25.0 - Nai_e)/3.0)))*(1.0/(1.0 + exp(8.0 - Ko_e))) # amp
I_diff_e = epsilon*(Ko_e-n_k0) + G_glia/(1.0 + exp((18.0-Ko_e)/2.5)) # amp
Kin_e	= 140.0+(18.0-Nai_e) # mole
Nao_e	= 144.0 - beta * (Nai_e - 18.0) # mole

V_k_e	= 26.64 * log(Ko_e/Kin_e) 							# volt
V_na_e	= 26.64 * log(Nao_e/Nai_e) 							# volt
V_l_e	= 26.64 * log ( (Ko_e + 0.065*Nao_e + 0.6*6.0) / (Kina_e + 0.065*Nai_e + 0.6*130.0))	# volt

Ie_k	= (g_k * (n_e * n_e * n_e * n_e) + g_ahp * Cai_e/(1.0+Cai_e) )*(Ve-V_k_e)	# amp
Ie_na	= g_na * (m_inf_Ve * m_inf_Ve * m_inf_Ve) * h_e * (Ve-V_na_e)			# amp
Ie_mem	= -g_l * (Ve - V_l_e) - Ie_na - Ie_k						# amp

dVedt	= (1/C) * (Ie_mem)		# volt # + Ie_ext + Ie_rand + Ie_syn
dn_edt	= phi * ( alpha_n_Ve * (1.0-n_e) - beta_n_Ve * n_e )	# volt
dh_edt	= phi * ( alpha_ra_Ve * (1.0-h_e) - beta_h_Ve * h_e )	# volt
ds_edt	= (1/tau_e) * (A * sigma_Ve * (1.0 - s_e) - s_e)	# volt
dCai_edt	= -0.002 * g_ca * (Ve - V_ca)/( 1.0 + exp( -1.0*(Ve+25.0)/2.5 ) ) - Cai_e/80.0	# mole

dKo_edt	= 0.001*(Ie_k/3.0 - 7.0* 2.0*I_pump_e - I_diff_e + diffusion_e)	# mole
dNai_edt	= 0.001*(-1.0*Ie_na/(7.0*3.0) - 3.0*I_pump_e)			# mole

deta_edt	= gamma_e * (Ve - v_b) - gamma_tilde * eta_e	# volt

#inhibitory
alpha_n_Vi	= 0.01 * (Vi+34.0)/( 1.0 - exp(-0.1 * (Vi+34.0)) ) # volt
beta_n_Vi	= 0.125 * exp(-1.0*(Vi+44.0)/80.0) # volt
alpha_m_Vi	= 0.1 * (Vi+30.0)/( 1.0 - exp(-0.1 * (Vi+30.0)) ) # volt
beta_m_Vi	= 4.0 * exp(-1.0*(Vi+55.0)/18.0) # volt
alpha_h_Vi	= 0.07 * exp(-1.0*(Vi+44.0)/20.0) # volt
beta_h_Vi	= 1.0/( 1.0 + exp(-0.1 * (Vi+14.0)) ) # volt
sigma_Vi	= 1.0/( 1.0 + exp(-1.0*(Vi+20.0)/4.0) ) # volt
m_inf_Vi	= alpha_m_Vi/(alpha_m_Vi + beta_m_Vi) # volt
I_pump_i	= (1.25/(1.0+exp((25.0-Nai_i)/3.0)))*(1.0/(1.0+exp(8.0-Ko_i)))	# amp
I_diff_i	= epsilon*(Ko_i-n_k0)+G_glia/(1.0 + exp((18.0-Ko_i)/2.5))		# amp
Kin_i		= 140.0+(18.0-Nai_i) 	# mole
Nao_i		= 144.0 - beta * (Nai_i - 18.0)	# mole

V_k_i	= 26.64 * log(Ko_i/Kin_i)	# volt
V_na_i	= 26.64 * log(Nao_i/Nai_i)	# volt
V_l_i	= 26.64 * log ( (Ko_i + 0.065*Nao_i + 0.6*6.0) / (Kin_i + 0.065*Nai_i + 0.6*130.0))	# volt

Ii_k	= (g_k * (n_i * n_i * n_i * n_i) + g_ahp * Cai_i/(1.0+Cai_i) )*(Vi-V_k_i)	# amp
Ii_na	= g_na * (m_inf_Vi * m_inf_Vi * m_inf_Vi) * h_i * (Vi-V_na_i) 	# amp
Ii_mem	= -g_l * (Vi - V_l_i) - Ii_na - Ii_k							# amp

dVidt	= (1/C) * (Ii_mem + Ii_ext + Ii_rand + Ii_syn)			# volt
dn_idt	= phi * ( alpha_n_Vi * (1.0-n_i) - beta_n_Vi * n_i )	# volt
dh_idt	= phi * ( alpha_h_Vi * (1.0-h_i) - beta_h_Vi * h_i )	# volt
ds_idt	= (1/tau_i) * (A * sigma_Vi * (1.0 - s_i) - s_i)		# volt

dKo_idt	= 0.001*(Ii_k/3.0 - 7.0* 2.0*I_pump_i - I_diff_i + diffusion_i) # mole
dNai_idt	= 0.001*(-1.0*Ii_na/(7.0*3.0) - 3.0*I_pump_i)				# mole

deta_idt	= gamma_i * (Vi - v_b) - gamma_tilde * eta_i	# volt


# TODO : Implement Ie_rand, Ie_syn, Ie_ext
# TODO : Implement Ii_rand, Ii_syn, Ii_ext
#specnum1	= sqrt(100.0/pi)
#specnum2	= sqrt(30.0/pi)		put in directly
# thresh		= 0.0
# alpha_ee	= 0.12
# alpha_ie	= 0.06
# alpha_ei	= 0.2 
# alpha_ii	= 0.02
# alpha_g	= 0.0
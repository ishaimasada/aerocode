import math

"""
NOTES

- Station Numbering convention: SAE ARP 755A
- 3 Spool Engine
- if a cooling flow does work --> enters at the inlet of a turbine
- if a cooling flow does not do work --> enters at the exit of the turbine
"""

def corrected_mdot(mdot, Tt, Pt):
    m_dot_corrected_a = mdot*math.sqrt(Tt/Ta)/(Pt/Pa)

# Inputs
Ta = 293 # K
Pa = 101.325 # kPa
Tt4 = 1610
hpc_ipc_split = 0.4
OPR = 35
BPR = 2.5
FPR = 2.1
CPR = OPR / FPR

pi_ipc = CPR**hpc_ipc_split
pi_hpc = CPR**(1 - hpc_ipc_split)
pi_b = 0.95
OGV_loss = 0.03
bypass_duct_loss = 0.01
fan_intake_tip_recovery = 0.98 # polyropic effciency of the fan tip intake
fan_intake_root_recovery =	0.997 # polyropic effciency of the fan root intake
e_tt_fan_tip = 0.9
e_tt_fan_root = 0.9
ess_Pt_loss = 0.04 # ESS total pressure loss
swan_loss =	0.02 # swan neck total pressure loss
e_tt_core_comp = 0.91 # core compressor total to total polytropic efficiency
eta_b =	0.99 # isentropic combustion efficiency
eta_tt_turbine = 0.91
eta_m =	0.99

# Unmixed Nozzle
'''
CD_core = 0.9
CV_core = 0.99
CD_fan = 0.9
CV_fan = 0.99
'''

# Mixed Nozzle
CV = 0.99
CD = 0.9

# Secondary Air System
hpt_cooling = 6.818181818 #does work (enters at inlets of stages)
hpt_packing = 0.681818182 # does no work
ipt_cooling = 4.090909091 # does work
ipt_packing = 0.681818182 # does no work
lpt_cooling = 2.045454545 # does work
lpt_packing = 0.681818182 # does no work
customer_offtake = 1
# Total SAS [%]	16


# Assumed Machs to get static properties
Ma = 0.95
M1 = 0.65
M2 = 0.5
M12 = 0.5
M3 = 0.3
M13 = 0.2

Cp = 10**3 # replace with polynomial correlation
gamma = 1.4 # replace with polynomial correlation
R = Cp * (gamma - 1) / gamma # replace with polynomial correlation
QR = 45*(10**6)


# Station 0: Ambient
mdota = 30 # kg/sec #Provide either mass flow or inlet area
Tta = 1 + ((gamma - 1) / 2) * Ma**2
Pta =  Pa * (Tta/Ta) ** (gamma / (gamma - 1))
mdot_corrected_a = corrected_mdot(mdota, Tta, Pta)

# Station 1: Inlet face
mdot1 = mdota
mdot_corrected_1 = mdot_corrected_a
Tt1 = Tta
Pt1 = Pta
ht1 = Tt1 * Cp

# Station 12: Fan Tip Inlet
mdot12 = mdota * (BPR / (BPR + 1))
Pt12 = Pt1 * fan_intake_tip_recovery
Tt12 = Tt1
ht12 = Tt12 * Cp
mdot_corrected_12 = corrected_mdot(mdot12, Tt12, Pt12)

# Station 13: Fan Tip Exit (Bypass Duct Inlet)
mdot13 = mdot12
Pt13 = Pt12 * FPR
Tt13 = Tt12 * (Pt13/Pt12)**((gamma - 1)/(gamma*e_tt_fan_tip))
ht13 = Tt13 * Cp
mdot_corrected_13 = corrected_mdot(mdot13, Tt13, Pt13)
fan_tip_power = mdot13*ht13 - mdot12*ht12


# Station 2: Fan Root 
mdot2 = mdot1
Tt2 = Tt1
ht2 = Tt2 * Cp
Pt2 =  Pt1 * fan_intake_root_recovery
mdot_corrected_2 = corrected_mdot(mdot2, Tt2, Pt2)

# Station 2.05: Fan Root Exit 
mdot205 = mdot2
Pt205 =  Pt1 * FPR
Tt205 = Tt2 * (Pt205/Pt2)**((gamma - 1) / (gamma * e_tt_fan_root))
ht205 = Tt205 * Cp
mdot_corrected_205 = corrected_mdot(mdot205, Tt205, Pt205)
fan_root_power = mdot205*ht205 - mdot2*ht2

# Station 2.1: IPC Inlet
mdot21 = mdot205
Pt21 =  Pt205*(1 - ess_Pt_loss)
Tt21 = Tt205 * (Pt21/Pt205)**((gamma - 1) / (gamma * e_tt_fan_root))
ht21 = Tt21 * Cp
mdot_corrected_21 = corrected_mdot(mdot21, Tt21, Pt21)

# Station 2.4: IPC Exit
m_offtake24 = -((ipt_cooling + ipt_packing + lpt_cooling + lpt_packing) / 100) * mdot2
mdot24 = mdot21 + m_offtake24
Pt24 =  Pt21*pi_ipc
Tt24 = Tt21 * (Pt24/Pt21)**((gamma - 1) / (gamma * e_tt_core_comp))
ht24 = Tt24 * Cp
mdot_corrected_24 = corrected_mdot(mdot24, Tt24, Pt24)
ipc_power = (mdot21*ht21 - mdot24*ht24) / eta_m

# Station 2.6: HPC Inlet
mdot26 = mdot24
Pt26 =  Pt21*pi_ipc
Tt26 = Tt24
ht26 = Tt26 * Cp
mdot_corrected_26 = corrected_mdot(mdot26, Tt26, Pt26)

# Station 3: HPC Exit
m_offtake3 = -((hpt_packing + hpt_cooling) / 100) * mdot2
mdot3 = mdot24 + m_offtake3
Pt3 = Pt26 * pi_hpc
Tt3 = Tt26 * (Pt3/Pt26 ** ((gamma - 1) / (gamma * e_tt_core_comp)))
ht3 = Tt3 * Cp
mdot_corrected_3 = corrected_mdot(mdot3, Tt3, Pt3)
hpc_power = (mdot26*ht26 - mdot3*ht3) / eta_m

# Station 3.1: Combustor Inlet
mdot31 = mdot24
Tt31 = Tt3
Pt31 = Pt3
mdot_corrected_31 = corrected_mdot(mdot31, Tt31, Pt31)

# Station 4: Combustor Exit
f4 = (Cp * (Tt4 - Tt31)) / ((QR * eta_b) - (Cp * Tt4))
mdot4 = mdot31 * (1 + f4)
Pt4 = Pt3 * pi_b
mdot_corrected_4 = corrected_mdot(mdot4, Tt4, Pt4)

# Station 405: HPT NGV Throat
m_offtake405 = (hpt_cooling / 100) * mdot2
mdot405 = mdot4 + m_offtake405
Pt405 = Pt4
Tt405 = (Tt4*mdot4 + Tt3*m_offtake405) / mdot405
mdot_corrected_405 = corrected_mdot(mdot405, Tt405, Pt405)
f405 = mdot31*f4 / (mdot405 - mdot31*f4)

# Station 41: HPT Rotor Inlet
mdot41 = mdot405
Pt41 = Pt4
Tt41 = Tt405
ht41 = Tt41 * Cp
mdot_corrected_41 = corrected_mdot(mdot41, Tt41, Pt41)
f41 = mdot31*f4 / (mdot41 - mdot31*f4)

# Station 44: HPT Exit
mdot44 = mdot405
ht44 = (ht41*mdot41 - hpc_power) / mdot41
Tt44 = ht44 / Cp
hpt_expansion = 1 / (Tt44 / Tt41)**((gamma - 1) / (gamma * eta_tt_turbine))
Pt44 = Pt4
mdot_corrected_44 = corrected_mdot(mdot44, Tt44, Pt44)
f44 = mdot31*f4 / (mdot44 - mdot31*f4)

# Station 45: IPT Inlet
m_offtake45 = (ipt_cooling + ipt_packing) * mdot21 / 100
mdot45 = mdot44 + m_offtake45
Tt45 = (Tt44*mdot44 + Tt24*m_offtake45) / mdot45
Pt45 = Pt4
ht45 = Tt45 * Cp
mdot_corrected_45 = corrected_mdot(mdot45, Tt45, Pt45)

# Station 47: IPT Exit
mdot47 = mdot45
Tt47 = (Tt44*mdot44 + Tt24*m_offtake45) / mdot45
ht47 = (Tt45*mdot45 - ipc_power) / mdot45
ipt_expansion = 1 / (Tt47 / Tt45)**((gamma - 1) / (gamma * eta_tt_turbine))
Pt47 = Pt45 / ipt_expansion
mdot_corrected_47 = corrected_mdot(mdot47, Tt47, Pt47)

# Station 48: LPT Inlet
m_offtake48 = (lpt_cooling + ipt_packing) * mdot21 / 100
mdot48 = mdot44 + m_offtake45
Pt48 = Pt47
Tt48 = (Tt47*mdot47 + Tt24*m_offtake48) / mdot48
ht48 = Tt48 * Cp
mdot_corrected_48 = corrected_mdot(mdot48, Tt48, Pt48)

# Station 5: LPT Outlet (Mixer Inlet from Core)
''' LPT packing flow is included here to avoid redundant stations '''
m_offtake5 = lpt_packing / 100 * mdot2
mdot5 = mdot48 + m_offtake5
ht5_pre_cooling = (Tt48*mdot48 - (fan_tip_power + fan_root_power)) / mdot48
Tt5_pre_cooling = ht5_pre_cooling / Cp
Tt5 = (Tt5_pre_cooling*mdot48 + Tt24*m_offtake5) / mdot5
ht5 = Tt5 * Cp
lpt_expansion = 1 / (Tt48 / Tt47)**((gamma - 1) / (gamma * eta_tt_turbine))
Pt5 = Pt48 / lpt_expansion
mdot_corrected_5 = corrected_mdot(mdot5, Tt5, Pt5)

# Station 15: Mixer Inlet from Fan
mdot15 = mdot13
Pt15 = Pt12 * (1 - OGV_loss - bypass_duct_loss)
Tt15 = Tt13
ht15 = Tt15 * Cp
mdot_corrected_15 = corrected_mdot(mdot15, Tt15, Pt15)

# Station 6: Mixer Exit
# Station 7: Nozzle Inlet
# Station 8: Nozzle Exit

Pt5_Pc = 1 / (1 - ((gamma - 1) / (gamma + 1)))**(gamma / (gamma - 1))
u_i = Ma * math.sqrt(gamma * R * Ta)
if (Pt5/Pa) > Pt5_Pc:
    print("Nozzle is choked\n")
    u_e = math.sqrt(2 * Cp * Tt5 * (1 - (1 / Pt5_Pc) ** ((gamma - 1) / gamma)))
else:
    print("Nozzle is unchoked\n")
    P5 = Pa
    T5 = Tt5 * (P5/Pt5)**((gamma - 1) / gamma)
    print(Tt5, T5)
    u_e = math.sqrt(2 * Cp * (Tt5 - T5))


""" Performance Parameters """
# Specific Thrust
T_ma = u_e * (1 + f44) - u_i
# Thrust
thrust = T_ma * mdota
# Thrust Specific Fuel Consumption (TSFC)
TSFC = (f44 / T_ma) * 10**3
# Propulsive Efficiency
eta_p = 2 / (1 + (u_e / u_i))
# Thermal Efficiency
eta_th = (((1 + f44) * u_e ** 2) - u_i**2) / (f44 * QR)
# Overall Efficiency
eta_o = eta_p * eta_th

print(f"TSFC: {TSFC}")
print(f"Specific Thrust: {T_ma}")
print(f"Thrust: {thrust}")
"""
Module for cycle Analysis and preliminary component design 
"""

import os
import sys
import numpy

# Change the current working directory to the file location
filepath = os.path.abspath(__file__)
directory = os.path.dirname(filepath)
os.chdir(directory)

# Add to search locations
sys.path.append(r'..\aerodynamics')

from compressible import *
from atmosphere import Ambient


class Station:
    """ Correlations for combustion products from Walsh and Fletcher """
    A = {
         "A0" : 992.313, "A1" : 236.688, "A2" : -1852.148, "A3" : 6083.152, "A4" : -8893.933, "A5" : 7097.112, 
         "A6" : -3234.725, "A7" : 794.571, "A8" : -81.873, "A9" : 422.178, "A10" : 1.053
        }

    B = {
         "B0" : -718.874, "B1": 8747.481, "B2": -15863.157, "B3": 17254.096, "B4": -10233.795, "B5": 3081.778, 
         "B6": -361.112, "B7": -3.919, "B8": 55.593, "B9": 1.6079
        }

    C = {"C0": 1.0001, "C1": 0.9248, "C2": -2.2078}
    
    REFH0 = 422.2202178
    Tref = 288.15
    Pref = 101325

    def __init__(self, W, Tt, Pt, name:dict=None):
        self.W = W
        self.Tt = Tt
        self.Pt = Pt
        self.name = name
        self.Wc = self.W * numpy.sqrt(self.Tt / self.Tref) / (self.Pt / self.Tref)

        # Statics
        self.M = self.T = self.P = self.V = self.rho = self.area = 0

        # Gas Model
        self.FAR = 0

    def __str__(self):
        return f"Station {self.name}\nW = {self.W}\nCorrected W = {self.Wc}\nTt = {self.Tt}\nPt = {self.Pt}"

    @property
    def R(self): return self.get_R(self.FAR)

    @property
    def cp(self): return self.get_cp(self.Tt, self.FAR)

    @property
    def gamma(self): return self.get_gamma()

    @property
    def ht(self): return self.get_ht(self.Tt, self.FAR)

    @R.setter
    def R(self, value): self._R = value

    @gamma.setter
    def gamma(self, value): self._gamma = value

    @ht.setter
    def ht(self, value): self._ht = value

    @cp.setter
    def cp(self, value): self._cp = value

    def get_R(self, FAR): return 287.05 - 0.0099 * FAR + 0.0000001 * FAR**2 # Formula 3.22 for Kerosene

    def get_gamma(self): return self.cp / (self.cp - self.R)

    def get_cp(self, T, FAR):
        TZ = T/1000
        cp = (
                self.A["A0"] + self.A["A1"] * TZ + (self.A["A2"] * TZ**2) + (self.A["A3"] * TZ**3) +
                (self.A["A4"] * TZ**4) + (self.A["A5"] * TZ**5) + (self.A["A6"] * TZ**6) +
                (self.A["A7"] * TZ**7) + (self.A["A8"] * TZ**8) + (self.B["B0"] + self.B["B1"] * TZ +
                (self.B["B2"] * TZ**2) + (self.B["B3"] * TZ**3) + (self.B["B4"] * TZ**4) +
                (self.B["B5"] * TZ**5) + (self.B["B6"] * TZ**6) + (self.B["B7"] * TZ**7)) * (FAR / (1 + FAR))
                ) # Formula 3.24
        return cp

        
    def set_statics(self, M):
        [_, Tt_T, Pt_P, _, _] = isentropic(M, lookup_key="M")
        self.M = M
        self.T = (1 / Tt_T) * self.Tt
        self.P = (1 / Pt_P) * self.Pt
        self.V = self.M * numpy.sqrt(self.gamma * self.R * self.T)
        self.rho = self.P / (self.R * self.T)
        self.area = self.W / (self.rho * self.V)


    def get_ht(self, T, FAR):
        """ Formula 3.27 """
        TZ = T/1000
        ht = (
              self.A["A0"] * TZ + self.A["A1"] * TZ**2 / 2 + (self.A["A2"] * TZ**3) / 3 + 
              (self.A["A3"] * TZ**4) / 4 + (self.A["A4"] * TZ**5) / 5 + (self.A["A5"] * TZ**6) / 6 +
              (self.A["A6"] * TZ**7) / 7 + (self.A["A7"] * TZ**8) / 8 + (self.A["A8"] * TZ**9) / 9 +
              self.A["A9"] + (self.B["B0"] * TZ + self.B["B1"] * TZ**2 / 2 + (self.B["B2"] * TZ**3) / 3 +
              (self.B["B3"] * TZ**4) / 4 + (self.B["B4"] * TZ**5) / 5 + (self.B["B5"] * TZ**6) / 6 +
              (self.B["B6"] * TZ**7) / 7 + (self.B["B7"] * TZ**8) / 8 + self.B["B8"]) * (FAR / (1 + FAR))
            ) - self.REFH0
        return ht


    def T_from_H(self, h, FAR, Thi, Tlo):
        """ Finding temperature from enthalpy polynomial using the Bisection Method """
        # Initial calculation of hmid to prevent uninitialized use
        Tmid = (Thi + Tlo) / 2
        hmid = self.get_ht(Tmid, FAR)
        error = abs(hmid - h) / h
        iterations = 0
        
        # Perform bisection
        while error > 0.001:
            Tmid = (Thi + Tlo) / 2
            hmid = self.get_ht(Tmid, FAR)
            
            if hmid < h: Tlo = Tmid
            elif hmid > h: Thi = Tmid
            
            iterations += 1
            T = Tmid
            error = abs(hmid - h) / h

        return T


class Inlet:
    def __init__(self, upstream:Station, Pt_recovery):
        self.upstream = upstream
        self.exit_W = self.upstream.W
        self.exit_Pt = self.upstream.Pt * Pt_recovery
        self.exit_Tt = self.upstream.Tt
        self.downstream = Station(self.exit_W, self.exit_Tt, self.exit_Pt)

class Compressor:
    def __init__(self, upstream:Station, PR):
        self.upstream = upstream
        self.exit_W = self.upstream.W
        self.exit_Pt = self.upstream.Pt * PR
        self.exit_Tt = self.upstream.Tt
        self.downstream = Station(self.exit_W, self.exit_Tt, self.exit_Pt)

class Burner:
    def __init__(self, upstream:Station, TET, LHV, pressure_loss):
        self.upstream = upstream
        self.eta_b = 1
        self.exit_FAR = self.upstream.get_FAR(TET, self.inlet.Tt, LHV)
        self.exit_W = self.upstream.W * (1 + self.exit_FAR)
        self.exit_Pt = upstream.Pt * (1 - pressure_loss)
        self.exit_Tt = TET
        self.downstream = Station(self.exit_W, self.exit_Tt, self.exit_Pt)
    

class Turbine:
    def __init__(self, cycle_parameters, cycle_specification=None):
        # CYCLE ANALYSIS
        self.upstream = cycle_parameters["upstream"] # Station object
        TET = cycle_parameters["TET"]
        e_tt = cycle_parameters["polytropic efficiency"]
        mdot_cool = cycle_parameters["coolant"]
        power = cycle_parameters["turbine power consumption"] * 1000 # MW to kW

        self.exit_W = self.upstream.W + mdot_cool
        self.exit_ht = (self.upstream.ht * self.upstream.W - power) / self.upstream.W
        FAR = self.upstream.W * self.upstream.FAR / (self.exit_W - (self.upstream.W * self.upstream.FAR))
        self.exit_Tt = self.upstream.T_from_H(self.exit_ht, FAR, TET, 0)
        ER = (self.exit_Tt  / self.upstream.Tt)**(-self.upstream.gamma / ((self.upstream.gamma - 1)*e_tt))
        self.exit_Pt = self.upstream.Pt / ER
        self.downstream = Station(self.exit_W, self.exit_Tt, self.exit_Pt)

        # COMPONENT DESIGN
        if cycle_specification != None:
            self.phi = cycle_specification["load coefficient"]
            self.psi = cycle_specification["work coefficient"]

    def velocity_triangle(self, R: list, rpm, Vax, Vu_mid, Rmid):
        omega = rpm * (2*numpy.pi / 60)
        U = omega * numpy.array(R)
        Vu = Vu_mid * (R / Rmid) # Free vortex assumption
        Wu = Vu - U
        V = numpy.sqrt(Vax**2 + Vu**2)
        W = numpy.sqrt(Vax**2 + Wu**2)
        alpha = numpy.atan(Vu / Vax)
        beta = numpy.atan(Wu / Vax)
        reaction = (W[2]**2 - W[1]**2) / (V[2]**2 - V[1]**2 + W[2]**2 - W[1]**2)
        #M_absolute = V / numpy.sqrt(self.downstream.gamma * self.downstream.R * self.downstream.T)
        #statics = isentropic()


class Nozzle:
    def __init__(self, upstream:Station):
        self.upstream = upstream
        self.downstream = Station(self.upstream.W, self.upstream.Pt, self.upstream.Tt, self.upstream.ht)

        def statics(self,  Pinf):
            gamma = self.downstream.gamma
            cp = self.downstream.cp
            R = self.downstream.R
            critical_NPR = (1 + ((gamma-1) / 2))**(gamma / (gamma-1))
            NPR = self.upstream.Pt / Pinf
            if NPR > critical_NPR:
                self.downstream.M = 1
                self.downstream.V = numpy.sqrt(gamma * R * self.downstream.T)
                self.downstream.P = self.downstream.Pt * critical_NPR
            else:
                self.downstream.P = Pinf
                self.downstream.T = self.downstream * (self.downstream.P/self.downstream.Pt)**((gamma - 1) / gamma)
                self.downstream.V = numpy.sqrt(2 * cp * (self.downstream.Tt - self.downstream.T))
                self.downstream.M = self.downstream.V / numpy.sqrt(gamma * R * self.downstream.T)

class Diffuser:
    def __init__(self, upstream:Station):
        self.upstream = upstream
        self.downstream = Station(self.upstream.W, self.upstream.Pt, self.upstream.Tt, self.upstream.ht)

        def set_statics(self,  Pa):
            pass

class Bleed:
    def init(self):
        pass
class Duct:
    def init(self):
        pass

class Engine:
    def __init__(self):
        # Design Parameters (Defaults are in metric units)
        self.OPR = 20
        self.TET = 1600
        self.BPR = 4
        self.FPR = 2.1
        self.LHV = 43.15*10**6
        self.CPR = self.OPR / self.FPR
        self.Minf = 0.85
        self.altitude = 13106.4

        # Other Parameters
        self.ambient = Ambient(self.Minf, self.altitude, self.alpha)
        self.front_diameter = 42 # inches
        self.FAR = 0
        self.alpha = 0

        # Component Efficiency Assumptions
        self.inlet_Pt_recovery_tip = 0.98
        self.inlet_Pt_recovery_hub = 0.97
        self.OGV_Pt_loss = 0.03
        self.bypass_duct_Pt_loss = 0.01
        self.e_tt_fan_tip = 0.9
        self.e_tt_fan_hub = 0.9
        self.swan_neck_loss = 0.02
        self.e_tt_core_comp = 0.91
        self.eta_b = 0.99
        self.burner_Pt_loss = 0.05
        self.e_tt_turbine = 0.8
        self.CD = 0.97
        self.CV = 0.99
        self.e_mechanical = 0.99

        # Prepare cycle analysis and component design inputs

        # Architecture
        ambient = Station()
        self.inlet = Inlet(upstream = ambient)
        self.compressor = Compressor(upstream = self.inlet.downstream)
        self.burner = Burner(upstream = self.compressor.downstream)
        self.turbine =  Turbine(upstream = self.burner.downstream)
        self.exhaust =  Nozzle(upstream = self.turbine.downstream)



    """ Performance Parameters """
    '''
    # Specific Thrust
    T_ma = u_e * (1 + f44) - u_i
    # Thrust Specific Fuel Consumption (TSFC)
    TSFC = (f44 / T_ma) * 10**3
    # Thrust
    thrust = T_ma * mdota

    # Propulsive Efficiency
    eta_p = 2 / (1 + (u_e / u_i))
    # Thermal Efficiency
    eta_th = (((1 + f44) * u_e ** 2) - u_i**2) / (f44 * QR)
    # Overall Efficiency
    eta_o = eta_p * eta_th
    '''

    def validate_cycle(self):
        pass
    def optimize(self, performance_parameter):
        pass
    def sensitivity_study(self):
        pass
    def off_design(self):
        pass
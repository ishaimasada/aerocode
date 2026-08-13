""" Module for Cycle Analysis and preliminary component design """

#NOTE: Axial Turbomachinery component classes assumes radially-constant work distribution (dht/dr = 0) and Free Vortex Design (dVax/dr = 0)

import matplotlib.pyplot as plt
import pandas
import numpy
import json
import copy
import sys
import os

filepath = os.path.abspath(__file__)
directory = os.path.dirname(filepath)
parent_directory = "\\".join(directory.split("\\")[:-1])
sys.path.append(parent_directory + r'\aerodynamics')

from compressible import bisection, isentropic, get_Tt_Ttstar # type: ignore
from atmosphere import Ambient # type: ignore
#from curves import *
os.chdir(directory)


class Station:
    """ Correlations for combustion products from Walsh and Fletcher """
    A = {"A0" : 992.313, "A1" : 236.688, "A2" : -1852.148, "A3" : 6083.152, "A4" : -8893.933, "A5" : 7097.112, 
         "A6" : -3234.725, "A7" : 794.571, "A8" : -81.873, "A9" : 422.178, "A10" : 1.053}

    B = {"B0" : -718.874, "B1": 8747.481, "B2": -15863.157, "B3": 17254.096, "B4": -10233.795, "B5": 3081.778, 
         "B6": -361.112, "B7": -3.919, "B8": 55.593, "B9": 1.6079}

    C = {"C0": 1.0001, "C1": 0.9248, "C2": -2.2078}
    
    REFH0 = 422.2202178
    Tref = 288.15
    Pref = 101325

    column_names = [
                    "Station",
                    "W [kg/sec]", 
                    "Tt [K]",
                    "Pt [kPa]",
                    "ht [kJ/kg]",
                    "Wc [kg/sec]", 
                    "Wf [kg/sec]", 
                    "FAR",
                    "M",
                    "Ts [K]",
                    "hs [kJ/kg]",
                    "Ps [kPa]",
                    "V [m/sec]",
                    "rho [kg/m^3]",
                    "Area [m^2]"
                   ]

    def __init__(self, W, Tt, Pt, FAR=None, M=None, idx=None):
        self.W = W
        self.Tt = Tt
        self.Pt = Pt
        self.Wc = self.get_Wc(self.W, self.Tt, self.Pt)
        self.M = M
        self.idx = idx

        if FAR != None: 
            self.FAR = FAR
            self.Wa = self.W / (1 + FAR)
            self.Wf = FAR * self.Wa
        else:
            self.FAR = 0
            self.Wf = 0
        if self.M != None: 
            self.set_statics(self.M)
        self.ht = self.get_ht(self.Tt, self.FAR)

    def __str__(self):
        return f"Station {self.idx}\nW = {self.W}\nWc = {self.Wc}\nCorrected W = {self.Wc}\nTt = {self.Tt}\nPt = {self.Pt/1000}\nM = {self.M}"

    @property
    def Wc(self): return self.get_Wc(self.W, self.Tt, self.Pt)

    @Wc.setter
    def Wc(self, value): self._Wc = value

    def get_Wc(self, W, T, P): return W * numpy.sqrt(T / self.Tref) / (P / self.Pref)

    @property
    def R(self): return self.get_R(self.FAR)

    @R.setter
    def R(self, value): self._R = value

    def get_R(self, FAR): return 287.05 - 0.0099 * FAR + 0.0000001 * FAR**2 # Formula 3.22 for Kerosene

    @property
    def gamma(self): return self.get_gamma()

    @gamma.setter
    def gamma(self, value): self._gamma = value

    def get_gamma(self): return self.cp / (self.cp - self.R)

    @property
    def FAR(self): return self.get_FAR()

    @FAR.setter
    def FAR(self, value): self._FAR = value

    def get_FAR(self): return self.Wf / (self.W - self.Wf)

    @property
    def cp(self): return self.get_cp(self.Tt, self.FAR)

    @cp.setter
    def cp(self, value): self._cp = value

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

    @property
    def ht(self): return self.get_ht(self.Tt, self.FAR)

    @ht.setter
    def ht(self, value): self._ht = value

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
        return ht*1000

        
    def set_statics(self, M):
        [_, Tt_T, Pt_P, _, _] = isentropic(M, self.gamma, lookup_key="M")
        self.M = M
        self.T = (1 / Tt_T) * self.Tt
        self.h = self.T * self.cp
        self.P = (1 / Pt_P) * self.Pt
        self.V = self.M * numpy.sqrt(self.gamma * self.R * self.T)
        self.rho = self.P / (self.R * self.T)
        if self.V == 0:
            self.area = numpy.inf
        else:
            self.area = self.W / (self.rho * self.V)


    def T_from_H(self, h, FAR, Thi, Tlo):
        """ Finding temperature from enthalpy polynomial using the Bisection Method """
        h_lo = self.get_ht(Tlo, FAR)
        h_hi = self.get_ht(Thi, FAR)
        if not (h_lo <= h <= h_hi or h_hi <= h <= h_lo): raise ValueError("Bisection bounds do not bracket the solution")

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

    def get_properties(self):
        if self.M != None: station_data = [self.idx, self.W, self.Tt, self.Pt/1000, self.ht/1000, self.Wc, self.Wf, self.FAR, self.M, self.T, self.h/1000, self.P/1000, self.V, self.rho, self.area]
        else: station_data = [self.idx, self.W, self.Tt, self.Pt/1000, self.ht/1000, self.Wc, self.Wf, self.FAR, None, None, None, None, None, None, None]
        return station_data


class Inlet:
    def __init__(self, cycle_parameters=None, component_parameters=None):
        if cycle_parameters != None:
            # CYCLE ANALYSIS
            M_inlet = cycle_parameters["inlet M"]
            M_exit = cycle_parameters["exit M"]
            Pt_recovery = cycle_parameters["total pressure recovery"]
            engine = cycle_parameters["engine"]
            freestream = engine.ambient

            if hasattr(engine, "fan_diameter"):
                D1 = engine.fan_diameter * 2.54 / 100 # in to m
                [_, Tt_T, Pt_P, rhot_rho, _] = isentropic(M_inlet, freestream.gamma, lookup_key="M")
                T1 = (1/Tt_T) * freestream.Tt
                a1 = numpy.sqrt(freestream.gamma * freestream.R * T1)
                V1 = M_inlet * a1
                A1 = numpy.pi * (D1 / 2)**2
                rho1 = freestream.rhot * (1/rhot_rho)
                W1 = V1 * rho1 * A1
            elif hasattr(engine, "W"):
                W1 = engine.W

            FAR = 0
            self.freestream = Station(W1, freestream.Tt, freestream.Pt, FAR, M=freestream.Minf, idx=0) # Station 0
            self.inlet = Station(W1, freestream.Tt, freestream.Pt, FAR, M=M_inlet, idx=1) # Station 1
            
            # Check if it's a BypassEngine
            if isinstance(engine, BypassEngine):
                bypass_W = W1*engine.B / (engine.B+1)
                root_W = W1 - bypass_W
                self.bypass_exit = Station(bypass_W, freestream.Tt, freestream.Pt, FAR, M=M_exit, idx=1.2)
                self.root_exit = Station(root_W, freestream.Tt, freestream.Pt, FAR, M=M_exit, idx=2)
            else:
                self.exit = Station(W1, freestream.Tt, freestream.Pt*Pt_recovery, FAR, M=M_exit, idx=2) # Station 1

        if component_parameters != None: 
            # COMPONENT DESIGN
            self.design_component(component_parameters)

    def design_component(self, component_parameters):
        self.component_parameters = component_parameters
        W0 = self.component_parameters["freestream"]["W"]
        Tt0 = self.component_parameters["freestream"]["Tt"]
        Pt0 = self.component_parameters["freestream"]["Pt"]
        gamma0 = self.component_parameters["freestream"]["gamma"]
        R0 = self.component_parameters["freestream"]["R"]
        cp0 = self.component_parameters["freestream"]["cp"]
        M0 = self.component_parameters["freestream"]["M"]
        M1 = self.component_parameters["M1"]
        Mth = self.component_parameters["Mth"]
        M_exit = self.component_parameters["M exit"]
        
        """ NOTE: N/R comes from a figure in Farohki. The user must read this value off of the graph based on the pressure coefficients """
        self.N_R = self.component_parameters["N/R"]
        self.CD = self.component_parameters["CD"]
        self.Zth = self.component_parameters["Zth"]
        self.freestream = Station(W0, Tt0, Pt0, M=M0, idx=0)
        self.inlet = Station(W0, Tt0, Pt0, M=M1, idx=1)
        self.throat = Station(W0, Tt0, Pt0, M=Mth, idx=1.2)
        self.exit = Station(W0, Tt0, Pt0, M=M_exit, idx=2)
        self.MFR = self.freestream.area / self.inlet.area
        self.AR =  self.exit.area / self.inlet.area
        self.CPR = (self.exit.P - self.inlet.P) / (self.inlet.Pt - self.inlet.P)
        self.CPR_ideal = 1 - (1 / self.AR**2)
        self.radii = self.get_radii()
        self.N = self.radii[0] * self.N_R
        self.z_coords = self.get_z_coords()


    def get_radii(self):
        r1 = numpy.sqrt(self.inlet.area / numpy.pi)
        rth = numpy.sqrt(self.throat.area / numpy.pi)
        r2 = numpy.sqrt(self.exit.area / numpy.pi)
        return [r1, rth, r2]
    
    def get_z_coords(self):
        z1 = 0
        zth = z1 + self.Zth
        z2 = zth + self.N
        return [z1, zth, z2]
    
    def plot_contour(self):
        plt.figure()
        plt.plot(self.z_coords, self.radii)
        plt.xlabel("Axial Distance [m]")
        plt.ylabel("Radius [m]")
        plt.title("Inlet Contour")
        plt.grid()
        plt.xlim([0, 0.1])
        plt.ylim([0, 0.1])
        #plt.legend()
        plt.show()


class Compressor:
    def __init__(self, upstream:Station, cycle_parameters=None, component_parameters=None, root_upstream=None):
        # CYCLE ANALYSIS
        if cycle_parameters != None:
            engine = cycle_parameters["engine"]
            self.e_t = cycle_parameters["e"]
            if "fan" in cycle_parameters:
                # Fan
                self.is_fan = True
                self.cooling = 0
                self.packing = 0
                self.B = engine.B
                self.rootPR = cycle_parameters["root PR"]
                self.tipPR = cycle_parameters["tip PR"]
                Mroot_exit = cycle_parameters["root exit M"]
                Mtip_exit = cycle_parameters["tip exit M"]
                self.root_inlet = copy.deepcopy(root_upstream)
                self.tip_inlet = copy.deepcopy(upstream)
                self.root_inlet.idx = "2"
                self.tip_inlet.idx = "1.2"
                self.root_exit = self.solve_exit(self.root_inlet, "2.05", self.rootPR, 0, 0, Mroot_exit)
                self.tip_exit = self.solve_exit(self.tip_inlet, "1.3", self.tipPR, 0, 0, Mtip_exit)
                root_delta_ht = self.root_exit.ht - self.root_inlet.ht
                tip_delta_ht = self.tip_exit.ht - self.tip_inlet.ht
                self.delta_ht = root_delta_ht + tip_delta_ht
                self.power = (self.root_inlet.W*root_delta_ht) + (self.tip_inlet.W*tip_delta_ht)
            else:
                # Compressor
                M_exit = cycle_parameters["exit M"]
                inlet_idx = cycle_parameters["inlet idx"]
                exit_idx = cycle_parameters["exit idx"]
                self.cooling = cycle_parameters["cooling"]
                self.packing = cycle_parameters["packing"]
                self.PR = cycle_parameters["PR"]
                self.inlet = upstream
                self.inlet.idx = inlet_idx
                self.exit = self.solve_exit(self.inlet, exit_idx, self.PR, self.cooling, self.packing, M_exit)
                self.delta_ht =  (self.exit.cp*self.exit.Tt) - (self.inlet.cp*self.inlet.Tt)
                self.power = self.inlet.W * self.delta_ht

        # COMPONENT DESIGN
        if component_parameters != None:
            self.design_component(component_parameters)

    def solve_exit(self, inlet, exit_idx, PR, cooling, packing, M_exit):
        exit = copy.deepcopy(inlet)
        exit.idx = exit_idx
        exit.M = M_exit
        exit.W = inlet.W * (1 - cooling - packing)
        exit.Pt = inlet.Pt * PR
        exit.Tt = inlet.Tt * (PR)**((inlet.gamma - 1) / (inlet.gamma*self.e_t))
        exit.set_statics(exit.M)
        return exit

    # Cycle method for performing preliminary component design ("specification" portion of the component design parameters)
    # Only call this after performing cycle analysis, NOT component design
    # NOTE: Missing the RPM / Nmech (% design point)
    def get_specification(self):
        specification = {
                "PR": self.PR,
                "power": self.power,
                "polytropic efficiency": self.e_t,
                "W": self.inlet.W,
                "Wc": self.inlet.Wc,
                "Tt": self.inlet.Tt,
                "Pt": self.inlet.Pt,
                "R": self.inlet.R,
                "gamma": self.inlet.gamma,
                "Cp": self.inlet.cp,
                "FAR": self.inlet.FAR
        },
        return specification

    def design_component(self, component_parameters):
        self.component_parameters = component_parameters
        flow = component_parameters["flow"]
        self.toggle_solution = self.component_parameters["toggle solution"]
        self.PR = self.component_parameters["specification"]["PR"]
        self.power = self.component_parameters["specification"]["power"]
        self.e_t = self.component_parameters["specification"]["polytropic efficiency"]
        self.rpm = self.component_parameters["rpm"]
        self.omega = self.rpm * ((2*numpy.pi/60))
        self.delta_ht = self.power / self.component_parameters["specification"]["W"]
        self.stages = list()
        match flow:
            case "axial":
                self.num_radii = self.component_parameters["number of radii"]
                self.U1t = self.component_parameters["U1t"]
                self.e_t = self.component_parameters["specification"]["polytropic efficiency"]
                self.average_dTt = self.component_parameters["average dTt"]
                # Estimate number of stages
                inlet_Tt = self.component_parameters["specification"]["Tt"]
                inlet_gamma = self.component_parameters["specification"]["gamma"]
                TR = self.PR**((inlet_gamma - 1)/inlet_gamma)
                dTt = inlet_Tt * (TR - 1)
                self.estimated_num_stages = round(dTt / self.average_dTt)
                # Ensure odd number of radii (must have a "mid-radius")
                if self.num_radii % 2 == 0 or self.num_radii < 3: raise ValueError("Number of radii per blade must be odd and greater than 3.")
                self.solve_axial()
            case "radial":
                self.integration = self.component_parameters["integration"]
                self.solve_radial()


    def solve_axial(self):
        if self.toggle_solution:
            # Handle stage inlets
            for idx, parameters in enumerate(self.component_parameters["stages"]):
                if idx == 0:
                    W = self.component_parameters["specification"]["W"]
                    Tt = self.component_parameters["specification"]["Tt"]
                    Pt = self.component_parameters["specification"]["Pt"]
                    R = self.component_parameters["specification"]["R"]
                    gamma = self.component_parameters["specification"]["gamma"]
                    Cp = self.component_parameters["specification"]["Cp"]
                    FAR = self.component_parameters["specification"]["FAR"]
                    M = self.component_parameters["specification"]["M"]
                    parameters["upstream"] = Station(W, Tt, Pt, FAR=FAR, M=M)
                else:
                    # Subsequent stages
                    parameters["upstream"] = self.stages[idx - 1].stations[3]
                self.stages.append(AxialStage(idx, "compressor", self))


    def solve_radial(self):
        if self.toggle_solution:
            # Handle stage inlets
            for idx, parameters in enumerate(self.component_parameters["stages"]):
                if idx == 0:
                    W = self.component_parameters["specification"]["W"]
                    Tt = self.component_parameters["specification"]["Tt"]
                    Pt = self.component_parameters["specification"]["Pt"]
                    R = self.component_parameters["specification"]["R"]
                    gamma = self.component_parameters["specification"]["gamma"]
                    Cp = self.component_parameters["specification"]["Cp"]
                    FAR = self.component_parameters["specification"]["FAR"]
                    M = self.component_parameters["specification"]["M"]
                    parameters["upstream"] = Station(W, Tt, Pt, FAR=FAR, M=M)
                else:
                    # Subsequent stages
                    parameters["upstream"] = self.stages[idx - 1].stations[3]
            self.stages.append(RadialStage(0, "compressor", self))

    def get_results(self):
        match self.integration.lower():
            case "cfturbo":
                pass
                # Outputs (CFTurbo)
                # print(self.compressor.PR)
                # print(s0.W)
                # print(D2)
                # print(self.b2)
                # isentropic_power = W * cp * Tt * (self.compressor.PR**((gamma - 1)/(gamma)) - 1)
                # eta_tt = (self.compressor.PR**((gamma-1)/gamma) - 1) / (self.compressor.PR**((gamma-1)/(gamma*self.compressor.e_t)) - 1) 
                # Dh1 = Rh1 * 2
                # Dc1 = Rc1 * 2
                # Uc1 = Rc1 * self.compressor.omega
                # diameter_ratio = (Rc1 * 2) / D2
            case "compal":
                pass
                # Outputs (Compal)
                # print(self.compressor.PR)
                # print(s0.W, s0.Tt, s0.Pt)
                # print(self.compressor.rpm)
                # print(self.exit_channel_type)
                # print(self.diffuser_type)
                # print(self.beta_blade)
                # print(self.diffuser_Ptloss)
            case "bladegen":
                z_coords = numpy.array([-self.stages[0].L, -self.stages[0].L, -self.stages[0].b2, 0])
                r_coords, _ = self.stages[0].get_meridional_coordinates()
                return r_coords, z_coords
                


class Burner:
    def __init__(self, upstream:Station, cycle_parameters=None, component_parameters=None):
        # CYCLE ANALYSIS
        if cycle_parameters != None:
            self.Ptloss_b = cycle_parameters["total pressure loss"]
            self.eta_b = cycle_parameters["efficiency"]
            self.M_exit = cycle_parameters["exit M"]
            self.engine = cycle_parameters["engine"]
            self.LHV = self.engine.LHV
            self.solve_exit(upstream)

        # COMPONENT DESIGN
        if component_parameters != None:
            self.design_component(component_parameters)
    
    # For Cycle Analaysis (not component design)
    def solve_exit(self, upstream):
        self.inlet = upstream
        self.exit = copy.deepcopy(self.inlet)
        self.exit.idx = 4
        self.exit.M = self.M_exit
        # Handle different choices for user input parameters
        if hasattr(self.engine, "TET"):
            TET = self.engine.TET
            FAR = self.get_FAR(TET, self.inlet.Tt, self.inlet.FAR, self.LHV, self.eta_b)
            self.exit.Wf = self.inlet.W * FAR
        elif hasattr(self.engine, "Wf"):
            self.exit.Wf = self.engine.Wf
            FAR = self.exit.Wf / upstream.W
            TET = bisection(self.get_FAR, FAR, 1800, 100, self.inlet.Tt, 0, self.LHV, self.eta_b)
        self.exit.W = self.inlet.W * (1 + self.exit.FAR)
        self.exit.Pt = self.inlet.Pt * (1 - self.Ptloss_b)
        self.exit.Tt = TET
        self.exit.set_statics(self.exit.M)

    # Cycle method
    def get_FAR(self, T2, T1, FAR1, LHV, eta):
        FARnew = 0.02
        FAR = -1 
        h1 = self.inlet.get_ht(T1, FAR1)
        tolerance = 0.00001
        error = (abs(FAR - FARnew) / FARnew) 
        while error > tolerance:
            FAR = FARnew
            h2 = self.exit.get_ht(T2, FAR)
            FARnew = (h2 - h1) / (LHV * eta)
            error = (abs(FAR - FARnew) / FARnew) 
        return FARnew
    
    # Component Design method
    def design_component(self, component_parameters):
        self.component_parameters = component_parameters
        self.machine = component_parameters["machine"]
        match self.machine.lower():
            case "annular":
                self.solve_annular()
            case "canannular":
                self.solve_canannular()
            case "can":
                self.solve_can()
    

    # Component Design method
    def solve_annular(self):
        # Design Parameters
        self.Pt_loss = self.component_parameters["Pt loss"]
        self.omega_cold = self.component_parameters["omega cold"]
        self.K_OTDF = self.component_parameters["K OTDF"]
        self.K_hot = self.component_parameters["K hot"]
        self.liner_area_frac = self.component_parameters["liner area fraction"]
        self.r_tip = self.component_parameters["r tip"]
        self.length = self.component_parameters["combustor length"]
        self.phi_PRZ = self.component_parameters["phi PRZ"]
        self.phi_SEC = self.component_parameters["phi SEC"]
        self.Tt_in = self.component_parameters["inlet Tt"]
        self.Pt_in = self.component_parameters["inlet Pt"]
        self.Tt_exit = self.component_parameters["exit Tt"]
        self.Pt_exit = self.component_parameters["exit Pt"]
        W = self.component_parameters["W"]
        Wf = self.component_parameters["Wf"]
        LHV = self.component_parameters["LHV"]
        R = self.component_parameters["R"]
        gamma = self.component_parameters["gamma"]
        FAR = Wf / (W - Wf)
        self.inlet = Station(W, self.Tt_in, self.Pt_in, FAR=FAR)
        self.exit = copy.deepcopy(self.inlet)
        self.exit.W += Wf
        self.exit.Tt = self.Tt_exit
        self.exit.Pt = self.Pt_exit

        # Calculations
        self.TR = self.Tt_exit / self.Tt_in
        self.omega_hot = self.K_hot * (self.TR - 1)
        self.omega_ref = self.omega_cold + self.omega_hot
        self.Aref = numpy.sqrt( (R/2) * (W * numpy.sqrt(self.Tt_in) / self.Pt_in)**2 * (self.omega_ref / self.Pt_loss) )
        if numpy.pi * self.r_tip**2 <= self.Aref:
            # Check if tip radius leads causes area to be smaller than reference area
            raise ValueError(f"Tip radius is too small. Increase tip radius so that area is less than reference area ({self.Aref}m^2)")
        FAR_stoichiometric = 1 / 15
        self.FAR_overall = Wf / W
        self.phi_overall = self.FAR_overall / FAR_stoichiometric
        self.dPt_check = self.omega_ref * (R/2) * (W * numpy.sqrt(self.Tt_in) / (self.Aref * self.Pt_in))**2
        self.rho_t3 = self.Pt_in / (R * self.Tt_in)
        self.Vref   = W / (self.rho_t3 * self.Aref)
        self.q_ref  = 0.5 * self.rho_t3 * self.Vref**2
        self.Mref = self.Vref / numpy.sqrt(gamma * R * self.Tt_in)
        self.Aliner = self.liner_area_frac * self.Aref
        self.r_hub = numpy.sqrt(self.r_tip**2 - self.Aref / numpy.pi)
        self.dl = (self.r_tip - self.r_hub)
        self.liner_r_tip = numpy.sqrt(self.Aliner / (2 * numpy.pi)) + (self.r_hub + self.dl/2)
        self.liner_r_hub = self.liner_r_tip - numpy.sqrt(self.Aliner / (2 * numpy.pi))
        self.Wa_PRZ = Wf / (self.phi_PRZ * FAR_stoichiometric)
        self.Wa_to_SEC = Wf / (self.phi_SEC * FAR_stoichiometric)
        self.Wa_SEC = self.Wa_to_SEC - self.Wa_PRZ
        self.Wa_DIL = W - self.Wa_to_SEC
        self.PRZ_Wa_W = (self.Wa_PRZ / W) * 100
        self.SEC_Wa_W = (self.Wa_SEC / W) * 100
        self.DIL_Wa_W = (self.Wa_DIL / W) * 100
        self.OTDF = 1 - numpy.exp(1 / (-self.K_OTDF * (self.length / self.dl) * self.omega_cold))
        self.volume = self.Aref * self.length
        self.tau_res = self.length / self.Vref
        self.tau_res_ms = self.tau_res * 1000
        Pt3_atm = self.Pt_in / 101325
        self.theta_i = (Wf * LHV) / (self.volume * Pt3_atm)
        self.theta_L = W / (self.volume * (Pt3_atm**1.8) * 10**(0.00145 * (self.Tt_in - 400)))
        self.eta_comb = (-5.46974e-10*self.theta_L**5) + (3.97923e-8*self.theta_L**4) - (8.73718e-6*self.theta_L**3) + (3.00007e-4*self.theta_L**2) - (4.568246e-3*self.theta_L) + 99.7
        self.Aheff = self.Aref / numpy.sqrt(self.omega_cold)
        # Check Results 
        all_pass = True
        if self.OTDF > 0.25:
            print('FLAG: OTDF = %.4f exceeds 0.25. Increase Lcomb or adjust liner sizing.\n', self.OTDF)
            all_pass = False
        if (self.theta_i / 1e6) > 60:
            print('FLAG: theta_i = %.2f MW/(m^3*atm) exceeds 60 SLS limit. Increase volume.\n', self.theta_i / 1e6)
            all_pass = False
        if self.theta_L > 5:
            print('FLAG: theta_L = %.4f exceeds 5 kg/(s*atm^1.8*m^3) SLS stability limit.\n', self.theta_L)
            all_pass = False
        if self.tau_res_ms < 3:
            print('FLAG: Residence time = %.3f ms is below 3 ms minimum. Increase Lcomb.\n', self.tau_res_ms)
            all_pass = False
        if self.eta_comb < 95:
            print('FLAG: Combustion efficiency = %.2f %% below 95%%. Reduce loading or increase volume.\n', self.eta_comb)
            all_pass = False
        if self.phi_PRZ < 0.8 or self.phi_PRZ > 1.3:
            print('FLAG: phi_PRZ = %.4f outside stable ignition range 0.8-1.3.\n', self.phi_PRZ)
            all_pass = False
        if all_pass == False:
            print("One or more checks failed, so the results have not been saved.")
    
    def solve_canannular(self):
        pass

    def solve_can(self):
        pass
    
    def get_results(self):
        thermo = {
            "Inlet Tt": self.Tt_in,
            "Exit Tt": self.Tt_exit,
            "Inlet Pt": self.Pt_in / 1e3,
            "Exit Pt": self.Pt_exit / 1e3,
            "Inlet Mass Flow": self.inlet.W,
            "Fuel Mass Flow": self.exit.Wf,
            "overall FAR": self.FAR_overall,
            "overall phi": self.phi_overall,
            "TR": self.TR,
            "omega cold": self.omega_cold,
            "omega hot": self.omega_hot,
            "omega ref": self.omega_ref,
            "dPt/Pt (target)": self.Pt_loss,
            "dPt/Pt (check from Aref)": self.dPt_check,
            "PRZ": [self.phi_PRZ, self.PRZ_Wa_W],
            "SEC": [self.phi_SEC, self.SEC_Wa_W],
            "DIL":  self.DIL_Wa_W,
            "OTDF": self.OTDF,
            "Residence time": self.tau_res_ms,
            "intensity": self.theta_i / 1e6,
            "loading": self.theta_L,
            "Combustion efficiency": self.eta_comb
        }
        geometry = {
            "Aref": self.Aref,
            "Vref":  self.Vref,
            "qref":  self.q_ref,
            "Mref": self.Mref,
            "Aliner": self.Aliner,
            "casing tip radius": self.r_tip,
            "casing hub radius": self.r_hub,
            "liner tip radius": self.liner_r_tip,
            "liner hub radius": self.liner_r_hub,
            "combustor length": self.length,
            "dl": self.dl,
            "Volume": self.volume,
            "Aheff": self.Aheff
        }
        return thermo, geometry


class Turbine:
    def __init__(self, upstream, compressor, cycle_parameters=None, component_parameters=None):
        if cycle_parameters != None:
            # CYCLE ANALYSIS
            self.cycle_parameters = cycle_parameters
            self.engine = self.cycle_parameters["engine"]
            self.upstream = upstream
            self.compressor = compressor
            self.e_t = self.cycle_parameters["e"]
            self.mdot_cool = self.compressor.cooling
            self.eta_m = self.cycle_parameters["mechanical efficiency"]
            self.power = abs(self.compressor.power) / self.eta_m
            self.exit_idx = self.cycle_parameters["exit idx"]
            self.M_exit = self.cycle_parameters["exit M"]
            self.cycle_parameters = self.cycle_parameters
            self.solve_exit()
        elif component_parameters != None:
            # COMPONENT DESIGN
            self.design_component(component_parameters)

    # Turbine Cycle Analysis (design-point)
    def solve_exit(self, upstream=None):
        if upstream != None:
            self.inlet = upstream
        else:
            self.inlet = self.upstream
        self.exit = copy.deepcopy(self.inlet)
        self.exit.idx = self.exit_idx
        self.exit.M = self.M_exit
        self.exit.W = self.inlet.W + self.mdot_cool
        exit_ht = (self.inlet.ht * self.inlet.W - (self.power/self.eta_m)) / self.exit.W
        FAR = self.inlet.W * self.inlet.FAR / (self.exit.W - (self.inlet.W * self.inlet.FAR))
        self.exit.Wf = (self.inlet.W - self.inlet.Wf) * FAR
        self.exit.Tt = self.inlet.T_from_H(exit_ht, self.exit.FAR, self.engine.TET, 100)
        self.ER = (self.exit.Tt  / self.inlet.Tt)**(-self.inlet.gamma / ((self.inlet.gamma - 1)*self.e_t))
        self.exit.Pt = self.inlet.Pt / self.ER
        self.exit.set_statics(self.exit.M)

    def design_component(self, component_parameters):
        self.component_parameters = component_parameters
        flow = self.component_parameters["flow"]
        self.rpm = self.component_parameters["rpm"] # Must ensure this matches the compressor rpm at each operating condition
        self.component_parameters = self.component_parameters
        self.specification = self.component_parameters["specification"]
        self.ER = self.specification["ER"]
        self.power = self.specification["power"]
        self.delta_ht = self.power / self.specification["W"] * 1000
        self.stages = list()
        match flow:
            case "axial":
                self.solve_axial()
            case "radial": 
                self.solve_radial()

    def solve_axial(self):
        for idx, parameters in enumerate(self.component_parameters["stages"]):
            # Handle the upstream station of each stage
            if idx == 0:
                # First Stage
                W = self.specification["W"]
                Wc = self.specification["ER"]
                Tt = self.specification["Tt"]
                Pt = self.specification["Pt"]
                R = self.specification["R"]
                gamma = self.specification["gamma"]
                Cp = self.specification["Cp"]
                FAR = self.specification["FAR"]
                #mid_radius = self.component_parameters["stages"][0]["U3m"] * self.omega
                #Vu = self.component_parameters["stages"][0]["Vu1m"]
                #mid = VelocityTriangle(label="0.5", radius=mid_radius, omega=self.omega, Vu=Vu, Vax, alpha, station=None, Mabs=None, Mrel=None)
                parameters["upstream"] = Station(W, Tt, Pt, FAR=FAR)
            else:
                # Subsequent stages
                parameters["upstream"] = self.stages[idx - 1].stations[3]
            self.stages.append(AxialStage(idx, "turbine", self))
    
    def solve_radial(self):
        return

    # Cycle method for performing preliminary component design ("specification" portion of the component design parameters)
    # Only call this after performing cycle analysis, NOT component design
    # NOTE: Missing the RPM / Nmech (% design point)
    def get_specification(self):
        specification = {
                "ER": self.ER,
                "power": self.power,
                "polytropic efficiency": self.e_t,
                "W": self.inlet.W,
                "Wc": self.inlet.Wc,
                "Tt": self.inlet.Tt,
                "Pt": self.inlet.Pt,
                "R": self.inlet.R,
                "gamma": self.inlet.gamma,
                "Cp": self.inlet.cp,
                "FAR": self.inlet.FAR
        }
        return specification


    # Report Component Design Results (plots and data)
    def get_results(self, flags):
        if flags["plots"]:
            # Velocity Triangles
            plt_idx = 1
            for stage in self.stages:
                for radius_idx in range(stage.num_radii):
                    for station_idx in range(3):
                        stage.stations[station_idx+1].triangles[radius_idx].plot_triangle(stage.num_radii+1, 3, plt_idx)
                        plt_idx += 1
                plt.savefig(f"Stage {stage.idx+1} Velocity Triangles.png")
                plt.clf()
            # Meridional View
            self.r_coords = numpy.array([])
            self.z_coords = numpy.array([])
            for stage_idx, stage in enumerate(self.stages):
                r_coords, z_coords = stage.get_meridional_coordinates()
                if stage.idx > 0:
                    # All axial coordinates of each blade start at 0 --> must shift each new blade row by stage axial spacing and the last axial position of the previous stage's rotor
                    previous_rotor_z, _ = self.stages[stage_idx - 1].rotor.get_axial_coords()
                    z_coords = [z + previous_rotor_z[-2] + stage.stator.axial_spacing for z in z_coords]
                self.r_coords = numpy.append(self.r_coords, r_coords)
                self.z_coords = numpy.append(self.z_coords, z_coords)
            
            plt.scatter(self.r_coords, self.z_coords)
            plt.savefig("Meridional View.png")
            plt.clf()
            
            # Blade Sections
            # NOTE: doing blade sections only makes sense if a 2D Euler solver was used to interate the curvature of each upper surface, lower surface, and camberline of each section at every radius --> skip for now
            '''
            plt.savefig("Blade Sections.png")
            plt.clf()
            '''
        if flags["data"]:
            velocities = {key: list() for key in ["V", "Vax", "Vu", "W", "Wu", "U", "Mabs", "Mrel", "T", "P"]}
            thermo = {key: list() for key in ["mdot", "Tt", "T", "Pt", "P"]}
            geometry = {key: list() for key in ["rm", "rt", "rh", "area"]}
            for stage in self.stages:
                # Stage Data
                _, stage_thermo, stage_geometry = stage.get_data(int((stage.num_radii-1) / 2))
                # Thermodynamics
                thermo["mdot"].append(stage_thermo["mdot"])
                thermo["Tt"].append(stage_thermo["Tt"])
                thermo["T"].append(stage_thermo["T"])
                thermo["Pt"].append(stage_thermo["Pt"])
                thermo["P"].append(stage_thermo["P"])
                # Geometry
                geometry["rm"].append(stage_geometry["rm"])
                geometry["rt"].append(stage_geometry["rt"])
                geometry["rh"].append(stage_geometry["rh"])
                geometry["area"].append(stage_geometry["area"])
                # Velocities
                for radius_idx in range(stage.num_radii):
                    radius_velocities, _, _ = stage.get_data(radius_idx)
                    velocities["V"].append(radius_velocities["V"])
                    velocities["Vax"].append(radius_velocities["Vax"])
                    velocities["Vu"].append(radius_velocities["Vu"])
                    velocities["W"].append(radius_velocities["W"])
                    velocities["Wu"].append(radius_velocities["Wu"])
                    velocities["U"].append(radius_velocities["U"])
                    velocities["Mabs"].append(radius_velocities["Mabs"])
                    velocities["Mrel"].append(radius_velocities["Mrel"])
            return velocities, thermo, geometry

    # Report Componenet Design Performance
    def display_performance(self):
        pass


# Velocity Triangle for axial turbomachines
class VelocityTriangle:
    def __init__(self, label, radius, omega, Vu, Vm, alpha, flow:str="axial", station=None, Mabs=None, Mrel=None):
        self.label = label
        self.radius = radius
        self.omega = omega
        self.U = self.radius * self.omega
        self.station = station
        self.alpha = alpha
        self.Vu = Vu

        # Absolute & Relative Velocities
        match flow:
            case "axial":
                self.Vr = 0
                self.Vax = Vm
                self.V = numpy.sqrt(self.Vu**2 + self.Vax**2)
                self.Wax = self.Vax
                self.Wu = self.Vu - self.U
                if self.Wu == 0:
                    self.W = 0
                else:
                    self.W = numpy.sqrt(self.Wu**2 + self.Vax**2)
                self.beta = numpy.atan(self.Wu / self.Wax)
            case "radial":
                self.Vax = 0
                self.Vr = Vm
                self.V = numpy.sqrt(self.Vu**2 + self.Vr**2)
                self.Wr = self.Vr
                self.Wu = self.Vu - self.U
                if self.Wu == 0:
                    self.W = 0
                else:
                    self.W = numpy.sqrt(self.Wu**2 + self.Wr**2)
                self.beta = numpy.atan(self.Wu / self.Wr)

        # Mach Numbers
        self.Mabs = Mabs
        self.Mrel = Mrel

        if self.station is not None:
            self.set_station(self.station)

    @property
    def Mabs(self): return self.get_Mabs()

    @Mabs.setter
    def Mabs(self, value): self._Mabs = value

    def get_Mabs(self): 
        if self.station is None:
            return None
        else:
            return self.V / numpy.sqrt(self.station.gamma * self.station.R * self.station.T)

    @property
    def Mrel(self): return self.get_Mrel()

    @Mrel.setter
    def Mrel(self, value): self._Mrel = value

    def get_Mrel(self): 
        if self.station is None:
            return None
        else:
            return self.W / numpy.sqrt(self.station.gamma * self.station.R * self.station.T)

    
    def plot_triangle(self, num_rows, num_columns, plt_idx):
        # Absolute Velocities
        plt.subplot(num_rows, num_columns, plt_idx)
        plt.arrow(0, 0, self.Vax, 0, width=2, fc="blue") # Vax
        plt.arrow(0, 0, self.Vax, self.Vu, width=2, fc='red') # V
        plt.arrow(self.Vax, 0, 0, self.Vu, width=1.5, linestyle='--', fc='red') # Vu
        
        # Relative Velocities
        plt.arrow(0, 0, self.Vax, self.Wu, width=2, fc='blue') # W
        plt.arrow(self.Vax, 0, 0, self.Wu, width=1.5, linestyle='--', fc='blue') # Wu
        
        # Wheel Velocity
        plt.arrow(self.Vax, self.Wu, 0, self.U, width=2.5, fc="brown")
   
        plt.title("Velocity Triangle", fontsize=12)
        plt.xlabel('Axial Velocity [m/s]')
        plt.ylabel('Tangential Velocity [m/s]')
        plt.legend(loc="best")

    def set_station(self, station):
        self.station = station
        self.T = self.station.Tt - self.V**2/self.station.cp
        self.P = self.station.Pt * (self.T / self.station.Tt)**(self.station.gamma/(self.station.gamma - 1))


# Stations for axial turbomachines (adds radii and velocity triangles to the engine-flow stations)
class AxialStation(Station):
    def __init__(self, idx, W, Tt, Pt, omega, mid:VelocityTriangle, num_radii=None, FAR=None, M=None):
        self.idx = idx
        self.mid = mid
        self.omega = omega
        self.num_radii = num_radii
        self.triangles = list()
        super().__init__(W, Tt, Pt, FAR=FAR, M=M, idx=idx)
        self.mid.set_station(self)
        

    def set_statics(self, M=None):
        if M is not None:
            self.M = M
        [_, Tt_T, Pt_P, _, _] = isentropic(self.M, self.gamma, lookup_key="M")
        self.T = (1 / Tt_T) * self.Tt
        self.h = self.T * self.cp
        self.P = (1 / Pt_P) * self.Pt
        self.V = self.M * numpy.sqrt(self.gamma * self.R * self.T)
        self.rho = self.P / (self.R * self.T)
        self.area = self.W / (self.rho * self.mid.Vax) # Calculate area based on the axial velocity since the velocity has a tangential component
        self.rhub = self.mid.radius - self.area/(4*numpy.pi*self.mid.radius)
        self.rtip = self.mid.radius + self.area/(4*numpy.pi*self.mid.radius)
        
    
    # Solve velocity triangles at every radius along the blade
    def define_velocity_triangles(self, num_radii=None):
        self.triangles.clear()
        if num_radii is not None:
            self.num_radii = num_radii
        self.radii = list(numpy.linspace(self.rhub, self.rtip, self.num_radii, endpoint=True, dtype=float))
        for idx, radius in enumerate(self.radii):
            U = radius * self.omega
            Vu = (self.mid.Vu * self.mid.radius) / radius # Free Vortex Equation
            alpha = numpy.atan(Vu / self.mid.Vax)
            label = f"{idx+1} / {num_radii}"
            self.triangles.append(VelocityTriangle(label, radius, self.omega, Vu, self.mid.Vax, alpha, station=self, flow="axial"))

            
    def plot_triangles(self):
        for triangle in self.triangles: 
            triangle.plot_triangle


class RadialStation(Station):
    def __init__(self, idx, W, Tt, Pt, triangle:VelocityTriangle=None, FAR=None, M=None):
        super().__init__(W, Tt, Pt, FAR=FAR, M=M, idx=idx)
        self.triangle = triangle
    

class BladeGeometry:
    def __init__(self, machine, flow, stage, blade:str, parameters:dict=None):
        self.machine = machine
        self.stage = stage
        self.blade = blade
        self.parameters = parameters

        if self.parameters is not None:
            match machine:
                case "turbine":
                    match flow:
                        case "axial":
                            self.AR = parameters["AR"]
                            self.zweiffel = parameters["zweiffel"]
                            match self.blade.lower():
                                case "stator":
                                    # AXIAL TURBINE STATOR 
                                    self.axial_turbine_stator()
                                case "rotor":
                                    # AXIAL TURBINE ROTOR 
                                    self.axial_turbine_rotor()
                        case "radial":
                            pass
                case "compressor":
                    match flow:
                        case "axial":
                            self.solidity = parameters["solidity"]
                            self.NOB = parameters["NOB"]
                            match self.blade.lower():
                                case "rotor":
                                    # AXIAL COMPRESSOR ROTOR 
                                    self.axial_compressor_rotor()
                                case "stator":
                                    # AXIAL COMPRESSOR STATOR 
                                    self.axial_compressor_stator()
                        case "radial":
                            match self.blade.lower():
                                case "rotor":
                                    # RADIAL COMPRESSOR ROTOR 
                                    self.radial_compressor_rotor()
                                case "stator":
                                    # RADIAL COMPRESSOR STATOR 
                                    self.radial_compressor_stator()


    def axial_turbine_stator(self):
        # Axial Turbine Stator
        s1 = self.stage.stations[1]
        s2 = self.stage.stations[2]
        s3 = self.stage.stations[3]
        S_rh_avg = (s1.rhub + s2.rhub) / 2
        S_rt_avg = (s1.rtip + s2.rtip) / 2
        self.HT_ratio = S_rh_avg / S_rt_avg
        self.h = S_rt_avg - S_rh_avg
        self.chord = self.h / self.AR
        
        # Stagger and Axial Chord
        self.stagger = [(s1.triangles[radius_idx].alpha + s2.triangles[radius_idx].alpha) / 2 for radius_idx in range(self.stage.num_radii)]
        self.cax = [self.chord * numpy.cos(self.stagger[radius_idx]) for radius_idx in range(self.stage.num_radii)]
        self.deflections = [s2.triangles[radius_idx].alpha - s1.triangles[radius_idx].alpha for radius_idx in range(self.stage.num_radii)]
        
        # Taper Ratio & Axial Spacing (between Stator and Rotor)
        self.taper_ratio = self.cax[-1] / self.cax[0]
        self.axial_spacing = 0.25 * self.cax[int((self.stage.num_radii + 1) / 2)]
        
        # Solidity and Pitch
        self.solidity = (2/self.zweiffel) * numpy.cos(s2.mid.alpha)**2 * (numpy.tan(s2.mid.alpha) - numpy.tan(s1.mid.alpha))
        self.pitch = self.chord / self.solidity

        # Number of Blades
        self.NOB = numpy.ceil((2 * numpy.pi * s1.mid.radius) / self.pitch)
        
        # Opening
        self.os = numpy.cos(s2.mid.alpha)
        self.opening = self.os * self.pitch


    def axial_turbine_rotor(self):
        # Axial Turbine Rotor 
        s1 = self.stage.stations[1]
        s2 = self.stage.stations[2]
        s3 = self.stage.stations[3]
        R_rh_avg = (s2.rhub + s3.rhub) / 2
        R_rt_avg = (s2.rtip + s3.rtip) / 2
        self.HT_ratio = R_rh_avg / R_rt_avg
        self.h = R_rt_avg - R_rh_avg
        self.chord = self.h / self.AR
        
        # Stagger and Axial Chord
        self.stagger = [(s2.triangles[radius_idx].beta + s3.triangles[radius_idx].beta) / 2 for radius_idx in range(self.stage.num_radii)]
        self.cax = [self.chord * numpy.cos(self.stagger[radius_idx]) for radius_idx in range(self.stage.num_radii)]
        self.deflections = [s3.triangles[radius_idx].beta - s2.triangles[radius_idx].beta for radius_idx in range(self.stage.num_radii)]
        
        # Taper Ratio
        self.taper_ratio = self.cax[-1] / self.cax[0]
        
        # Solidity and Pitch
        self.solidity = (2/self.zweiffel) * numpy.cos(s3.mid.beta)**2 * (numpy.tan(s3.mid.beta) - numpy.tan(s2.mid.beta))
        self.pitch = self.chord / self.solidity

        # Number of Blades
        self.NOB = numpy.ceil((2 * numpy.pi * s1.mid.radius) / self.pitch)
        
        # Opening
        self.os = numpy.cos(s2.mid.beta)
        self.opening = self.os * self.pitch

    def axial_compressor_rotor(self):
        # Axial Compressor Rotor 
        s1 = self.stage.stations[1]
        s2 = self.stage.stations[2]
        s3 = self.stage.stations[3]
        rh_avg = (s1.rhub + s2.rhub) / 2
        rt_avg = (s1.rtip + s2.rtip) / 2
        self.HT_ratio = rh_avg / rt_avg
        self.h = rt_avg - rh_avg

        # Pitch, chord, & AR
        self.pitch = numpy.ceil((2 * numpy.pi * s1.mid.radius) / self.NOB)
        self.chord = self.pitch * self.solidity
        self.AR = self.h / self.chord
        
        # Stagger & Axial Chord
        self.stagger = [(s1.triangles[radius_idx].beta + s2.triangles[radius_idx].beta) / 2 for radius_idx in range(self.stage.compressor.num_radii)]
        self.cax = [self.chord * numpy.cos(self.stagger[radius_idx]) for radius_idx in range(self.stage.compressor.num_radii)]
        self.deflections = [s2.triangles[radius_idx].beta - s1.triangles[radius_idx].beta for radius_idx in range(self.stage.compressor.num_radii)]
        
        # Taper Ratio
        self.taper_ratio = self.cax[-1] / self.cax[0]
        
        # Solidity and Pitch
        self.zweiffel = (2/self.solidity) * numpy.cos(s2.mid.beta)**2 * (numpy.tan(s2.mid.beta) - numpy.tan(s1.mid.beta))
        
    def axial_compressor_stator(self):
        # Axial Compressor Rotor 
        s1 = self.stage.stations[1]
        s2 = self.stage.stations[2]
        s3 = self.stage.stations[3]
        rh_avg = (s2.rhub + s3.rhub) / 2
        rt_avg = (s2.rtip + s3.rtip) / 2
        self.HT_ratio = rh_avg / rt_avg
        self.h = rt_avg - rh_avg

        # Pitch, chord, & AR
        self.pitch = numpy.ceil((2 * numpy.pi * s2.mid.radius) / self.NOB)
        self.chord = self.pitch * self.solidity
        self.AR = self.h / self.chord
        
        # Stagger & Axial Chord
        self.stagger = [(s2.triangles[radius_idx].alpha + s3.triangles[radius_idx].alpha) / 2 for radius_idx in range(self.stage.compressor.num_radii)]
        self.cax = [self.chord * numpy.cos(self.stagger[radius_idx]) for radius_idx in range(self.stage.compressor.num_radii)]
        self.deflections = [s3.triangles[radius_idx].alpha - s2.triangles[radius_idx].alpha for radius_idx in range(self.stage.compressor.num_radii)]
        
        # Taper Ratio
        self.taper_ratio = self.cax[-1] / self.cax[0]
        
        # Solidity and Pitch
        self.zweiffel = (2/self.solidity) * numpy.cos(s3.mid.alpha)**2 * (numpy.tan(s3.mid.alpha) - numpy.tan(s2.mid.alpha))

    def get_axial_coords(self):
        # This is for an axial blade (needs revision to distinguish between axial and radial blades)
        z_coords = []
        for idx in range(self.stage.num_radii * 2):
            if idx == 0:
                # Hub Leading Edge (Corner)
                z_coords.append(0)
            elif idx < self.stage.num_radii:
                # Leading Edge
                z_coords.append((self.cax[0] - self.cax[idx]) /  2)
            else:
                # Trailing Edge
                # Formula: leading Edge z-coordinate + axial chord length (leading edge coordinate from a radius index is the number of radii minus the given radius index)
                radius_idx = (2*self.stage.num_radii) - 1 - idx
                z_coords.append(z_coords[radius_idx] + self.cax[radius_idx])
        return z_coords

    def radial_compressor_rotor(self):
        pass

    def radial_compressor_stator(self):
        pass


# General axial-stage object to be used for turbines & compressors (handles stage calculations, not be used outside of component classes)
class AxialStage:
    def __init__(self, idx, machine, component):
        self.idx = idx
        match machine.lower():
            case "turbine":
                # AXIAL TURBINE
                self.turbine = component
                self.turbine_parameters = self.turbine.component_parameters
                self.stage_parameters = self.turbine.component_parameters["stages"][idx]
                self.solve_turbine()
            case "compressor": 
                # AXIAL COMPRESSOR
                self.compressor = component
                self.compressor_parameters = self.compressor.component_parameters
                self.stage_parameters = self.compressor.component_parameters["stages"][idx]
                self.solve_compressor()

    def solve_turbine(self):
        # Load in stage parameters
        self.aerodynamics = self.stage_parameters["aerodynamics"]
        self.geometry = self.stage_parameters["geometry"]
        self.cooling = self.stage_parameters["cooling"]
        self.upstream = self.stage_parameters["upstream"]
        self.num_radii = self.turbine_parameters["number of radii"]
        self.omega = self.turbine.rpm * (2*numpy.pi / 60)
        efficiency = self.turbine_parameters["specification"]["polytropic efficiency"]
        phi = self.aerodynamics["flow coefficient"]
        psi = self.aerodynamics["loading coefficient"]
        loss_coefficient = self.aerodynamics["loss coefficient"]
        M2m = self.aerodynamics["M2m"]
        U3m = self.aerodynamics["U3m"]
        Vu1m = self.aerodynamics["Vu1m"]
        Vax2_Vax1 = self.aerodynamics["Vax climb stator"]
        Vax3_Vax2 = self.aerodynamics["Vax climb rotor"]
        Rm2_Rm1 = self.geometry["rm climb rate stator"]
        Rm3_Rm2 = self.geometry["rm climb rate rotor"]
        AR_stator = self.geometry["AR stator"]
        AR_rotor = self.geometry["AR rotor"]
        zweiffel = self.geometry["Zweiffel"]

        # Ensure odd number of radii (must have a "mid-radius")
        if self.num_radii % 2 == 0 or self.num_radii < 3: raise ValueError("Number of radii per blade must be odd and greater than 3.")

        # Radii & Axial Velocities
        Rm3 = U3m / self.omega
        Rm2 = Rm3 / Rm3_Rm2
        Vax3 = U3m * phi
        Vax2 = Vax3 / Vax3_Vax2
        if self.idx == 0:
            # First Stage
            Rm1 = Rm2 / Rm2_Rm1
            Vax1 = Vax2 / Vax2_Vax1
        else:
            # Subsequent Stages
            Rm1 = self.upstream.radii[-1]
            Vax1 = self.upstream.mid.Vax
        
        # Stage Quantites
        self.delta_ht = psi * (Rm3*self.omega)**2
        self.work_split = self.delta_ht / self.turbine.delta_ht
        self.capacity = self.upstream.W * numpy.sqrt(self.upstream.Tt) / (self.upstream.Pt / 101.325)

        # Meanline and Radial Calculations (Velocity Triangles)
        # Station 1
        W1 = self.upstream.W
        Tt1 = self.upstream.Tt
        Pt1 = self.upstream.Pt
        alpha1 = numpy.atan(Vu1m / Vax1)
        mid1 = VelocityTriangle("station 1 mid", Rm1, self.omega, Vu1m, Vax1, alpha1, flow="axial")
        s1 = AxialStation(1, W1, Tt1, Pt1, self.omega, mid=mid1, num_radii=self.num_radii)
        s1.T = s1.Tt - s1.mid.V**2/(2*s1.cp)
        s1.M = numpy.sqrt((2/(s1.gamma - 1)) * (s1.Tt/s1.T) - 1)
        s1.set_statics()
        s1.define_velocity_triangles(self.num_radii)

        # Station 2
        s2 = copy.deepcopy(s1)
        s2.idx = 2
        s2.Pt = s1.Pt - loss_coefficient*(s1.Pt - s1.P)
        s2.set_statics(M2m)
        V2 = M2m * numpy.sqrt(s2.gamma * s2.R * s2.T)
        alpha2 = numpy.acos(Vax2 / V2)
        Vu2 = V2 * numpy.sin(alpha2)
        s2.mid = VelocityTriangle("station 2 mid", Rm2, self.omega, Vu2, Vax2, alpha2, flow="axial")
        s2.mid.set_station(s2)
        s2.define_velocity_triangles()

        # Station 3
        s3 = copy.deepcopy(s2)
        s3.idx = 3
        s3.Tt = s2.Tt - self.delta_ht/s2.cp
        s3.Pt = s2.Pt * (s3.Tt/s2.Tt)**(s3.gamma/(efficiency*(s3.gamma - 1)))
        Vu3 = (self.delta_ht - Rm2*self.omega*s2.mid.Vu ) / (Rm3*self.omega) # Euler Turbine Equation
        alpha3 = numpy.atan(Vu3/Vax3)
        s3.mid = VelocityTriangle("station 3 mid", Rm3, self.omega, Vu3, Vax3, alpha3, flow="axial")
        s3.mid.set_station(s3)
        s3.T = s3.Tt - s3.mid.V**2/(2*s3.cp)
        s3.M = numpy.sqrt((2/(s3.gamma - 1)) * (s3.Tt/s3.T) - 1)
        s3.set_statics()
        s3.define_velocity_triangles()

        # Store stations in dictionary object
        self.stations = {1: s1, 2: s2, 3: s3}

        # Degrees of Reaction (DoR)
        self.DoR = list()
        for radius_idx in range(self.num_radii):
            self.DoR.append(self.get_DoR(radius_idx, machine="turbine"))

        # Deflections
        self.deflections = list()
        for radius_idx in range(self.num_radii):
            self.deflections.append(self.get_deflection(radius_idx, machine="turbine"))

        # Expansion Ratio (Inlet Pt / Exit Pt)
        self.ER = self.stations[1].Pt / self.stations[3].Pt

        # Blade Geometries
        self.stator = BladeGeometry(machine="turbine", flow="axial", stage=self, blade="stator", parameters={"AR": AR_stator, "zweiffel": zweiffel})
        self.rotor = BladeGeometry(machine="turbine", flow="axial", stage=self, blade="rotor", parameters={"AR": AR_rotor, "zweiffel": zweiffel})

        # Cooling
        self.solve_turbine_cooling()

    
    def get_DoR(self, radius_idx, machine):
        ht1 = self.stations[1].ht
        ht3 = self.stations[3].ht
        match machine.lower():
            case "turbine":
                h2 = self.stations[2].ht - self.stations[2].triangles[radius_idx].V**2/2
                h3 = self.stations[3].ht - self.stations[3].triangles[radius_idx].V**2/2
                return (h2 - h3) / abs(ht3 - ht1)
            case "compressor":
                h2 = self.stations[2].ht - self.stations[2].triangles[radius_idx].V**2/2
                h1 = self.stations[1].ht - self.stations[1].triangles[radius_idx].V**2/2
                return (h2 - h1) / abs(ht3 - ht1)


    def get_deflection(self, radius_idx, machine):
        match machine.lower():
            case "turbine":
                stator_deflection = self.stations[2].triangles[radius_idx].alpha - self.stations[1].triangles[radius_idx].alpha
                rotor_deflection = self.stations[3].triangles[radius_idx].beta - self.stations[2].triangles[radius_idx].beta
                return {"stator": stator_deflection, "rotor": rotor_deflection}
            case "compressor":
                rotor_deflection = self.stations[2].triangles[radius_idx].beta - self.stations[1].triangles[radius_idx].beta
                stator_deflection = self.stations[3].triangles[radius_idx].alpha - self.stations[2].triangles[radius_idx].alpha
                return {"rotor": rotor_deflection, "stator": stator_deflection}


    def solve_turbine_cooling(self):
        s1 = self.stations[1]
        s2 = self.stations[2]
        s3 = self.stations[3]
        OTDF = self.cooling["OTDF"]
        RTDF = self.cooling["RTDF"]
        CDT = self.cooling["CDT"]
        cp_coolant = self.cooling["cp coolant"]
        metal_TtMax = self.cooling["metal TtMax"]
        cooling_efficiency = self.cooling["cooling efficiency"]

        # STATOR
        # Peak gas temperature seen by stator (using Inlet Total Temp)
        T_peak_S = OTDF * (s1.Tt - CDT) + s1.Tt
        
        # Corrected temperature based on recovery factor (0.15 logic from original)
        T_peak_S_corr = T_peak_S - 0.15 * (s2.mid.V)**2 / (2 * cp_coolant)
        
        # Required cooling effectiveness
        eps_req_S = (T_peak_S_corr - metal_TtMax) / (T_peak_S_corr - CDT)
        
        # Non-dimensional coolant parameter
        K_cool_S = eps_req_S / (cooling_efficiency * (1 - eps_req_S))
        
        # Heat Transfer: Stator
        Area_S = 2 * self.stator.chord * self.stator.h * self.stator.NOB
        mu_S = 1.458e-6 * T_peak_S_corr**(3/2) / (T_peak_S_corr + 110.4)
        Re_S = (s1.rho * s1.mid.V * self.stator.chord) / mu_S
        k_S = 0.000053983 * T_peak_S_corr + 0.013568
        Nu_S = 0.488 * Re_S**0.592
        h_conv_S = Nu_S * k_S / self.stator.chord
        
        # Resulting Stator Coolant Flow
        self.stator.mdot_cool = K_cool_S * h_conv_S * Area_S / cp_coolant

        # ROTOR
        # Peak gas temperature seen by rotor (using RTDF)
        T_peak_R = RTDF * (s1.Tt - CDT) + s1.Tt
        
        # Corrected temperature based on relative velocity (W) for the rotor
        T_peak_R_corr = T_peak_R - 0.15 * (s2.mid.W)**2 / (2 * cp_coolant)
        
        # Required cooling effectiveness
        eps_req_R = (T_peak_R_corr - metal_TtMax) / (T_peak_R_corr - CDT)
        
        # Non-dimensional coolant parameter
        K_cool_R = eps_req_R / (cooling_efficiency * (1 - eps_req_R))
        
        # Heat Transfer: Rotor
        Area_R = 2 * self.rotor.chord * self.rotor.h * self.rotor.NOB
        mu_R = 1.458e-6 * T_peak_R_corr**(3/2) / (T_peak_R_corr + 110.4)
        Re_R = (s2.rho * s2.mid.W * self.rotor.chord) / mu_R
        k_R = 0.000053983 * T_peak_R_corr + 0.013568
        Nu_R = 0.488 * Re_R**0.592
        h_conv_R = Nu_R * k_R / self.rotor.chord
        
        # Resulting Rotor Coolant Flow
        self.rotor.mdot_cool = K_cool_R * h_conv_R * Area_R / cp_coolant
        
        # Total Stage Cooling
        self.total_mdot_cool = self.stator.mdot_cool + self.rotor.mdot_cool

    def get_dehaller(self):
        pass

    def get_diffusion(self):
        pass
    
    def solve_compressor(self):
        # Load in stage parameters
        efficiency = self.compressor_parameters["specification"]["polytropic efficiency"]
        loss_coefficient = self.stage_parameters["aerodynamics"]["loss coefficient"]
        upstream = self.stage_parameters["upstream"]
        M3m = self.stage_parameters["aerodynamics"]["exit M"]
        alpha3m = self.stage_parameters["aerodynamics"]["exit alpha"]
        Vax2_Vax1 = self.stage_parameters["aerodynamics"]["Vax climb stator"]
        Vax3_Vax2 = self.stage_parameters["aerodynamics"]["Vax climb rotor"]
        Rm2_Rm1 = self.stage_parameters["geometry"]["rm climb rate stator"]
        Rm3_Rm2 = self.stage_parameters["geometry"]["rm climb rate rotor"]
        solidity_rotor = self.stage_parameters["geometry"]["stator solidity"]
        solidity_stator = self.stage_parameters["geometry"]["rotor solidity"]
        NOB_stator = self.stage_parameters["geometry"]["stator NOB"]
        NOB_rotor = self.stage_parameters["geometry"]["rotor NOB"]
        dTt = self.stage_parameters["dTt"]

        # Radii & Axial Velocities
        if self.idx == 0:
            # First Stage
            Rt1 = self.compressor.U1t / self.compressor.omega
            Rh1 = numpy.sqrt(Rt1**2 - upstream.area/numpy.pi)
            Rm1 = (Rt1 - Rh1) / 2 + Rh1
            Vax1 = upstream.V
            alpha1 = self.compressor_parameters["specification"]["alpha"]
        else:
            # Subsequent Stages
            Rm1 = upstream.radii[-1]
            Vax1 = upstream.mid.Vax
            alpha1 = upstream.mid.alpha
        Rm2 = Rm2_Rm1 * Rm1
        Rm3 = Rm3_Rm2 * Rm2
        Vax2 = Vax2_Vax1 * Vax1
        Vax3 = Vax3_Vax2 * Vax2
        
        # Stage Quantites
        self.delta_ht = dTt * upstream.cp
        self.work_split = self.delta_ht / self.compressor.delta_ht
        self.capacity = upstream.W * numpy.sqrt(upstream.Tt) / (upstream.Pt / 101.325)

        # Meanline and Radial Calculations (Velocity Triangles)
        # Station 1
        W1 = upstream.W
        Tt1 = upstream.Tt
        Pt1 = upstream.Pt
        Vu1m = upstream.V * alpha1
        mid1 = VelocityTriangle("station 1 mid", Rm1, self.compressor.omega, Vu1m, Vax1, alpha1, flow="axial")
        s1 = AxialStation(1, W1, Tt1, Pt1, self.compressor.omega, mid=mid1, num_radii=self.compressor.num_radii)
        s1.T = s1.Tt - s1.mid.V**2/(2*s1.cp)
        s1.M = numpy.sqrt((2/(s1.gamma - 1)) * (s1.Tt/s1.T) - 1)
        s1.set_statics()
        s1.define_velocity_triangles(self.compressor.num_radii)

        # Station 2
        s2 = copy.deepcopy(s1)
        s2.idx = 2
        s2.Tt = s1.Tt + dTt
        s2.Pt = s1.Pt * (s2.Tt/s1.Tt)**(s2.gamma/(efficiency*(s2.gamma - 1)))
        Vu2 = (Rm1*self.compressor.omega*s1.mid.Vu - self.delta_ht) / (Rm3*self.compressor.omega) # Euler Turbine Equation
        alpha2 = numpy.atan(Vu2/Vax2)
        s2.mid = VelocityTriangle("station 2 mid", Rm2, self.compressor.omega, Vu2, Vax2, alpha2, flow="axial")
        s2.mid.set_station(s2)
        s2.T = s2.Tt - s2.mid.V**2/(2*s2.cp)
        s2.M = numpy.sqrt((2/(s2.gamma - 1)) * (s2.Tt/s2.T) - 1)
        s2.set_statics()
        s2.define_velocity_triangles()

        # Station 3
        s3 = copy.deepcopy(s2)
        s3.idx = 3
        s3.Pt = s2.Pt - loss_coefficient*(s2.Pt - s2.P)
        s3.set_statics(M3m)
        V3 = M3m * numpy.sqrt(s2.gamma * s2.R * s2.T)
        Vu3 = V3 * numpy.sin(alpha3m)
        s3.mid = VelocityTriangle("station 3 mid", Rm3, self.compressor.omega, Vu3, Vax3, alpha3m, flow="axial")
        s3.mid.set_station(s3)
        s3.define_velocity_triangles()

        # Store stations in dictionary object
        self.stations = {1: s1, 2: s2, 3: s3}

        # Degrees of Reaction (DoR)
        self.DoR = list()
        for radius_idx in range(self.compressor.num_radii):
            self.DoR.append(self.get_DoR(radius_idx, machine="compressor"))

        # Deflections
        self.deflections = list()
        for radius_idx in range(self.compressor.num_radii):
            self.deflections.append(self.get_deflection(radius_idx, machine="compressor"))

        # Expansion Ratio (Inlet Pt / Exit Pt)
        self.ER = self.stations[1].Pt / self.stations[3].Pt

        # Blade Geometries
        self.rotor = BladeGeometry(machine="compressor", flow="axial", stage=self, blade="rotor", parameters={"solidity": solidity_rotor, "NOB": NOB_rotor})
        self.stator = BladeGeometry(machine="compressor", flow="axial", stage=self, blade="stator", parameters={"solidity": solidity_stator, "NOB": NOB_stator})


    def get_meridional_coordinates(self):
        # Radial Coordinates
        r_coords = self.stations[1].radii + self.stations[2].radii + self.stations[2].radii + self.stations[3].radii
        # Axial Coordinates
        stator_z = self.stator.get_axial_coords()
        rotor_z = [z + self.stator.axial_spacing for z in self.rotor.get_axial_coords()]
        z_coords = stator_z + rotor_z
        return r_coords, z_coords


    def get_data(self, radius_idx):
        velocities = {key: list() for key in ["V", "Vax", "Vu", "W", "Wu", "U", "Mabs", "Mrel", "T", "P"]}
        thermo = {key: list() for key in ["mdot", "Tt", "T", "Pt", "P"]}
        geometry = {key: list() for key in ["rm", "rt", "rh", "area"]}
        for station in self.stations.values():
            velocities["V"].append(station.triangles[radius_idx].V)
            velocities["Vax"].append(station.triangles[radius_idx].Vax)
            velocities["Vu"].append(station.triangles[radius_idx].Vu)
            velocities["W"].append(station.triangles[radius_idx].W)
            velocities["Wu"].append(station.triangles[radius_idx].Wu)
            velocities["U"].append(station.triangles[radius_idx].U)
            velocities["Mabs"].append(station.triangles[radius_idx].Mabs)
            velocities["Mrel"].append(station.triangles[radius_idx].Mrel)
            velocities["T"].append(station.triangles[radius_idx].T)
            velocities["P"].append(station.triangles[radius_idx].P)
            thermo["mdot"].append(station.W)
            thermo["Tt"].append(station.Tt)
            thermo["Pt"].append(station.Pt)
            thermo["T"].append(station.mid.T)
            thermo["P"].append(station.mid.P)
            geometry["rm"].append(station.mid.radius)
            geometry["rt"].append(station.radii[-1])
            geometry["rh"].append(station.radii[0])
            geometry["area"].append(station.area)
        return velocities, thermo, geometry


# General axial-stage object to be used for turbines & compressors (handles stage calculations, not be used outside of component classes)
class RadialStage:
    def __init__(self, idx, machine, component):
        self.idx = idx
        match machine.lower():
            case "turbine":
                # RADIAL TURBINE
                self.turbine = component
                self.turbine_parameters = self.turbine.component_parameters
                self.stage_parameters = self.turbine.component_parameters["stages"][idx]
                self.solve_turbine()
            case "compressor": 
                # RADIAL COMPRESSOR
                self.compressor = component
                self.compressor_parameters = self.compressor.component_parameters
                self.stage_parameters = self.compressor.component_parameters["stages"][idx]
                self.solve_compressor()
    
    def solve_turbine(self):
        return
    
    def solve_compressor(self):
        W = self.compressor_parameters["specification"]["W"]
        Wc = self.compressor_parameters["specification"]["Wc"]
        Tt = self.compressor_parameters["specification"]["Tt"]
        Pt = self.compressor_parameters["specification"]["Pt"]
        R = self.compressor_parameters["specification"]["R"]
        gamma = self.compressor_parameters["specification"]["gamma"]
        cp = self.compressor_parameters["specification"]["Cp"]
        FAR = self.compressor_parameters["specification"]["FAR"]
        M = self.compressor_parameters["specification"]["M"]
        self.diffuser_type = self.compressor_parameters["diffuser type"]
        self.exit_channel_type = self.compressor_parameters["exit channel type"]
        self.slip_model = self.compressor_parameters["slip model"]
        self.psi = self.stage_parameters["aerodynamics"]["loading coefficient"]
        self.beta_blade = numpy.deg2rad(self.stage_parameters["aerodynamics"]["exit beta"])
        self.blockage_model = self.stage_parameters["aerodynamics"]["blockage model"]
        self.diffuser_Ptloss = self.stage_parameters["aerodynamics"]["diffuser Pt loss"]
        self.L = self.stage_parameters["geometry"]["impeller height"]
        self.impeller_NOB = self.stage_parameters["geometry"]["impeller NOB"]
        self.diffuser_NOB = self.stage_parameters["geometry"]["diffuser NOB"]
        self.nu_h1 = self.stage_parameters["geometry"]["hub diameter ratio"]
        self.nu_c1 = self.stage_parameters["geometry"]["casing diameter ratio"]
        self.width_ratio = self.stage_parameters["geometry"]["exit width ratio"]
        self.t2 = self.stage_parameters["geometry"]["trailing blade hub thickness"]
        self.diffuser_width = self.stage_parameters["geometry"]["diffuser width"]
        self.U2 = numpy.sqrt(self.compressor.delta_ht / self.psi)
        self.R2 = self.U2 / self.compressor.omega

        # Inducer (Station 1)
        D2 = self.R2 * 2
        s0 = Station(W, Tt, Pt, FAR=FAR, M=M) # Exit of inlet (psuedo station "0")
        [_, _, _, _, A0_Astar] = isentropic(s0.M, s0.gamma)
        Rh1 = (1/2) * self.nu_h1 * D2
        Rc1 = (1/2) * self.nu_c1 * D2
        R1_mid = Rh1 + (Rc1 - Rh1)/2
        Astar = (1/A0_Astar) * s0.area
        A1 = numpy.pi * (Rc1**2 - Rh1**2)
        if A1 < Astar:
            raise ValueError("Impeller inlet area cannot be smaller than critical area --> increase diameter ratios.")
        A1_Astar = A1 * A0_Astar / s0.area
        [M1, Tt_T, Pt_P, _, _] = isentropic(A1_Astar, s0.gamma, lookup_key="area ratio")
        rho1 = (s0.Pt / Pt_P) / (R * (s0.Tt * Tt_T))
        Vax1 = s0.W / (rho1 * A1)
        mid1 = VelocityTriangle("inlet mid radius", R1_mid, self.compressor.omega, Vu=0, Vm=Vax1, alpha=0, flow="axial")
        s1 = AxialStation(1, W, Tt, Pt, self.compressor.omega, mid=mid1, num_radii=3, FAR=FAR, M=M1)
        s1.define_velocity_triangles()

        # Impeller (Station 2)
        s2 = RadialStation(2, W, Tt, Pt, FAR=FAR)
        s2.Tt = s2.T_from_H(s1.ht+self.compressor.delta_ht, s2.FAR, 1800, 0)
        s2.Pt = self.compressor.PR * s1.Pt
        self.b2 = self.width_ratio * D2
        A_effective = self.calculate_blockage()
        Vu2 = self.compressor.delta_ht / self.U2
        Vu2_real = self.calculate_slip(Vu2, self.beta_blade)
        Vr2 = self.iterate_radial_velocity(A_effective, 100, Vu2_real, s2)
        Wu2 = Vu2 - self.U2
        beta_tangential = abs(numpy.atan(Vr2 / Wu2)) # FROM TANGENTIAL DIRECTION --> this is CFTurbo's convention for alpha and beta
        beta_radial = abs(numpy.atan(Wu2 / Vr2)) # FROM RADIAL (or Meridional) DIRECTION --> this is standard convention for alpha and beta
        alpha_tangential = numpy.atan(Vr2 / Vu2_real)
        alpha_radial = numpy.atan(Vu2_real / Vr2)
        triangle2 = VelocityTriangle("impeller exit", self.R2, self.compressor.omega, Vu=Vu2, Vm=Vr2, alpha=alpha_radial, flow="radial") # this was made using the radial alpha 
        s2.triangle = triangle2

        rho2 = s2.W / (A_effective * Vr2)
        rhot2 = s2.Pt / (R * s2.Tt)
        P2 = s2.Pt * (rho2 / rhot2)**(s2.gamma)

        # Diffuser (Station 3)
        s3 = copy.deepcopy(s2)
        s3.idx = 3
        s3.Pt = s2.Pt * (1 - self.diffuser_Ptloss)
        match self.diffuser_type:
            case "vaneless":
                R3 = self.R2 + self.diffuser_width
                A3 = 2 * numpy.pi * R3 * self.b2
                Vr3 = Vr2 * self.R2 / R3
                Vu3 = Vu2 * self.R2 / R3
                alpha3 = numpy.atan(Vr3/Vu3)
                triangle3 = VelocityTriangle("diffuser exit", R3, 0, Vu=Vu3, Vm=Vr3, alpha=alpha3, flow="radial")
                s3.triangle = triangle3
            case "vaned":
                pass
            case "both":
                # Assumes a vaneless diffuser and then a vaned diffuser
                pass
        
        # Exit Channel (Station 4)
        match self.exit_channel_type:
            case "return channel":
                pass
            case "volute":
                pass

        self.stations = {1: s1, 2: s2, 3: s3}
        self.impeller = BladeGeometry(machine="compressor", flow="radial", stage=self, blade="rotor", parameters=self.stage_parameters["geometry"])
        self.diffuser = BladeGeometry(machine="compressor", flow="radial", stage=self, blade="stator", parameters=self.stage_parameters["geometry"])

    def calculate_blockage(self):
        match self.blockage_model:
            case "tangential":
                effective_area = self.b2 * (2*numpy.pi*self.R2 -  (self.impeller_NOB*self.t2/numpy.sin(self.beta_blade))) 
        return effective_area

    def iterate_radial_velocity(self, A, Vr_guess, Vu, station):
        tolerance = 0.0001
        W, ht, Tt, Pt, gamma, FAR, R = station.W, station.ht, station.Tt, station.Pt, station.gamma, station.FAR, station.R
        V = numpy.sqrt(Vr_guess**2 + Vu**2)
        h = ht - V**2/2
        T = station.T_from_H(h, FAR, 1800, 0)
        P = Pt * (T/Tt)**(gamma / (gamma - 1))
        rho = P / (R*T)
        Vr_new = W / (rho * A)
        error = abs(Vr_new - Vr_guess)
        while error >= tolerance:
            Vr_guess = copy.deepcopy(Vr_new)
            V = numpy.sqrt(Vr_guess**2 + Vu**2)
            h = ht - V**2/2
            T = station.T_from_H(h, FAR, 1800, 0)
            P = Pt * (T/Tt)**(gamma / (gamma - 1))
            rho = P / (R*T)
            Vr_new = W / (rho * A)
            error = abs(Vr_new - Vr_guess)
        return Vr_new
    
    def calculate_slip(self, *args):
        match self.slip_model:
            case "radial":
                Vu_real = self.radial_slip(args)
            case "Aungier-Wiesner":
                Vu_real = self.aungier_wiesner_slip(*args)
        return Vu_real
    
    def radial_slip(self, Vu_ideal):
        epsilon = 1 - (1.98 / self.impeller_NOB)
        Vu_real = Vu_ideal * epsilon
        return Vu_real

    def aungier_wiesner_slip(self, Vu_ideal, beta_blade):
        outflow_coefficient = 1 - numpy.sqrt(numpy.sin(beta_blade) / (self.impeller_NOB**0.7))
        Vu_real = Vu_ideal - (1 - outflow_coefficient)*self.U2
        return Vu_real

    def get_meridional_coordinates(self):
        # Currently only for the impeller (for BladeGen integration)
        r_coords = list()
        z_coords = list()
        # Radial Coordinates
        r_coords.extend(self.stations[1].radii + [self.stations[2].triangle.radius, self.stations[2].triangle.radius])
        # Axial Coordinates
        z_coords.extend([0, 0, self.L, self.L + self.b2])
        return r_coords, z_coords



class Mixer:
    def __init__(self, hot_inlet:Station, cold_inlet:Station, component_parameters=None): 
        # CYCLE ANALYSIS
        # NOTE The hot stream Mach Number is chosen, but the cold stream Mach Number can equivalently be chosen
        # Hot Stream (static properties)
        hot_inlet.T = hot_inlet.Tt * (1 + (hot_inlet.gamma - 1)/2 * hot_inlet.M**2)**(-1)
        hot_inlet.P = hot_inlet.Pt * (1 + (hot_inlet.gamma - 1)/2 * hot_inlet.M**2)**(-hot_inlet.gamma / (hot_inlet.gamma - 1))
        hot_inlet.rho = hot_inlet.P / (hot_inlet.R*hot_inlet.T)
        hot_inlet.V = hot_inlet.M * numpy.sqrt(hot_inlet.gamma * hot_inlet.R * hot_inlet.T)
        hot_inlet.area = hot_inlet.W / (hot_inlet.rho*hot_inlet.V)
        self.hot_inlet = hot_inlet
        momentum_hot = hot_inlet.W*hot_inlet.V + hot_inlet.area*hot_inlet.P

        # Cold Stream (static properties)
        cold_inlet.P = hot_inlet.P
        cold_inlet.M = numpy.sqrt(2 / (cold_inlet.gamma - 1) * ((cold_inlet.Pt/cold_inlet.P)**((cold_inlet.gamma - 1)/cold_inlet.gamma) - 1))
        cold_inlet.T = cold_inlet.Tt * (1 + (cold_inlet.gamma - 1)/2 * cold_inlet.M**2)**(-1)
        cold_inlet.rho = cold_inlet.P / (cold_inlet.R*cold_inlet.T)
        cold_inlet.V = cold_inlet.M * numpy.sqrt(cold_inlet.gamma * cold_inlet.R * cold_inlet.T)
        cold_inlet.area = cold_inlet.W / (cold_inlet.rho*cold_inlet.V)
        self.cold_inlet = cold_inlet
        momentum_cold = cold_inlet.W*cold_inlet.V + cold_inlet.area*cold_inlet.P

        # Design Parameter (hot stream Mach Number) Restriction
        if cold_inlet.P > cold_inlet.Pt: raise ValueError("Static pressure is larger than the total pressure (not physically possible). A larger Mhot must be chosen.")

        # Mixer Exit
        idx = "6"
        W = hot_inlet.W + cold_inlet.W
        Wf = hot_inlet.Wf
        A = hot_inlet.area + cold_inlet.area
        Tt = (hot_inlet.W*hot_inlet.Tt + cold_inlet.W*cold_inlet.Tt) / W
        FAR = Wf / (W - Wf)
        R = hot_inlet.get_R(FAR)
        Cp = hot_inlet.get_cp(Tt, FAR)
        gamma = Cp / (Cp - R)
        K = (momentum_hot + momentum_cold)**2 * (gamma / (W**2 * Tt * R))
        M = numpy.sqrt(((2*gamma - K) + numpy.sqrt((K - 2*gamma)**2 + 4*(((gamma-1)/2)*K-gamma**2))) / (2*((gamma-1)/2*K-gamma**2)))
        T = Tt * (1 + (gamma - 1)/2 * M**2)**(-1)
        P = (W*numpy.sqrt(R*T)/(A*M*numpy.sqrt(gamma)))
        Pt = P * (Tt/T)**(gamma / (gamma - 1))
        self.exit = Station(W, Tt, Pt, FAR=FAR, idx=idx, M=M)

        # COMPONENT DESIGN
        if component_parameters != None: pass


class Afterburner:
    def __init__(self, upstream:Station=None, cycle_parameters=None, component_parameters=None):
        if cycle_parameters != None:
            # CYCLE ANALYSIS
            engine = cycle_parameters["engine"]
            self.toggle = cycle_parameters["toggle"]
            self.Ttmax = cycle_parameters["AET"]
            self.pi_ab_on = cycle_parameters["pi_hot"]
            self.pi_ab_off = cycle_parameters["pi_cold"]
            self.eta = cycle_parameters["efficiency"]
            self.LHV = engine.LHV
            self.solve_exit(upstream)

        if component_parameters != None: 
            # COMPONENT DESIGN
            self.design_component(component_parameters)

    def solve_exit(self, upstream=None):
        if upstream != None:
            self.inlet = upstream
        else:
            self.inlet = self.upstream
        self.exit = copy.deepcopy(self.inlet)
        self.exit.idx = "7"
        # Handle different engine modes (afterburner on or off)
        match self.toggle: 
            case False:
                self.exit.Pt = self.inlet.Pt * self.pi_ab_off # Afterburner off
            case True: # Afterburner on
                self.exit.Pt = self.inlet.Pt * self.pi_ab_on
                Tt_Ttstar1 = get_Tt_Ttstar(self.inlet.M, self.inlet.gamma)
                Tt_Ttstar2 = self.Ttmax * (Tt_Ttstar1) / self.inlet.Tt
                Ttstar = 1 / Tt_Ttstar1 * self.inlet.Tt

                # Check if the afteburner chokes before reaching the given total
                # temperature (based on Rayleigh flow theory)
                if Tt_Ttstar2 > 1:
                    print(f"Given Afterburner Exit Temperature exceeds maximum allowed amount {Ttstar}.\n"
                        "This value has been set as the exit total temperatue and a lower value is recommended to be chosen.\n\n")
                    self.exit.Tt = Ttstar
                    self.exit.M = 1
                elif Tt_Ttstar2 <= 1:
                    self.exit.Tt = self.Ttmax
                    self.exit.M = bisection(get_Tt_Ttstar, Tt_Ttstar2, 1, 0, self.inlet.gamma)

                self.exit.FAR = Burner.get_FAR(self, self.exit.Tt, self.inlet.Tt, self.inlet.FAR, self.LHV, self.eta) 
                self.exit.Wf = (self.inlet.W - self.inlet.Wf) * self.exit.FAR
                self.exit.W = self.inlet.W + self.exit.Wf
                self.exit.set_statics(self.exit.M)

    def design_component(self, component_parameters):
        self.component_parameters = component_parameters
        match self.component_parameters["flameholder geometry"].lower():
            case "vee gutter":
                self.solve_vee_gutter()
            case "assymetric":
                self.solve_assymetric()


    def solve_vee_gutter(self):
        # "Dry Mode" is when no fuel is added --> exit gamma is the same as inlet gamma
        # "Wet Mode" is when fuel is added --> exit gamma is not the same as inlet gamma (Tt increases from heat)
        # Equations from Farohki (FOURTH EDITION)
        self.N = self.component_parameters["num rings"]
        self.CD = self.component_parameters["drag coefficient"]
        self.eta = self.component_parameters["efficiency"]
        self.LHV = self.component_parameters["lower heating value"]
        self.d = self.component_parameters["base thickness"]
        self.di = self.component_parameters["inner diameter"]
        self.alpha = self.component_parameters["apex angle"]
        W = self.component_parameters["inlet"]["W"]
        Tt = self.component_parameters["inlet"]["Tt"]
        Pt = self.component_parameters["inlet"]["Pt"]
        FAR = self.component_parameters["inlet"]["FAR"]
        M = self.component_parameters["inlet"]["M"]
        self.inlet = Station(W, Tt, Pt, FAR=FAR, M=M)
        self.exit = copy.deepcopy(self.inlet)
        self.exit.Tt = self.component_parameters["Tt exit"]
        self.exit.FAR = Burner.get_FAR(self, self.exit.Tt, self.inlet.Tt, self.inlet.FAR, self.LHV, self.eta)
        self.q = self.exit.FAR * self.LHV * self.eta
        # Calculations
        Adry = self.get_A(0)
        Awet = self.get_A(self.q)
        Me_dry = self.get_Me(Adry, self.exit.gamma)
        Me_wet = self.get_Me(Awet, self.exit.gamma)
        self.Pt_loss_dry = self.get_loss(self.inlet.M, Me_dry, self.inlet.gamma, self.inlet.gamma)
        self.Pt_loss_wet = self.get_loss(self.inlet.M, Me_wet, self.inlet.gamma, self.exit.gamma)

    # Vee-Gutter Method
    def get_Me(self, A, gamma):
        return numpy.sqrt(((2*(gamma - 1) - gamma*A**2 + numpy.sqrt((gamma*A**2 - 2*gamma + 2)**2 + 2*((gamma - 1)**2)*(A**2 - 2)))) / ((gamma - 1)*(A**2 - 2)))

    # Vee-Gutter Method
    def get_A(self, q):
        Mi = self.inlet.M
        gamma_i = self.inlet.gamma
        return ((1 + gamma_i*Mi**2*(1-self.CD/2)) / (gamma_i*Mi)) * numpy.sqrt(((gamma_i - 1)/(1 + (gamma_i - 1)/2*Mi**2)) * (1 / (1 + (q/self.inlet.ht))))

    # Vee-Gutter Method --> Eq. 5.99
    def get_loss(self, Mi, Me, gamma_i, gamma_e):
        return (1 + gamma_i*Mi**2*(1 - self.CD/2)*((1 + (gamma_e - 1)/2*Me**2)**(gamma_e/(gamma_e-1)))) / ((1 + gamma_e*Me**2)*((1 + (gamma_i - 1)/2*Mi**2)**(gamma_i/(gamma_i-1))))

    def solve_assymetric(self):
        pass


class Nozzle:
    def __init__(self, upstream:Station, cycle_parameters=None, component_parameters=None):
        if cycle_parameters != None:
            # CYCLE ANALYSIS
            self.cycle_parameters = cycle_parameters
            self.geometry = self.cycle_parameters["C/CD"]
            self.engine = self.cycle_parameters["engine"]
            self.Pinf = self.engine.ambient.P
            self.CD = self.cycle_parameters["discharge coefficient"]
            self.CV = self.cycle_parameters["velocity coefficient"]
            self.solve_exit(upstream)
            
            # COMPONENT DESIGN
            if component_parameters != None:
                self.design_component(component_parameters)

    # Cycle method
    def solve_exit(self, upstream):
        # Handle Converging or Converging-Diverging nozzle geometries
        match self.geometry:
            case "C": self.solve_converging(upstream)
            case "CD": self.solve_CD(upstream)

    # Cycle method
    def solve_converging(self, upstream):
        self.inlet = upstream
        gamma = self.inlet.gamma
        cp = self.inlet.cp
        R = self.inlet.R
        critical_NPR = (1 + ((gamma-1) / 2))**(gamma / (gamma-1))
        self.NPR = self.inlet.Pt / self.Pinf
        self.exit = copy.deepcopy(self.inlet)
        self.exit.idx = 8 # Throat station
        if self.NPR >= critical_NPR:
            # Choked
            self.exit.M = 1
            self.exit.T = (2 / (gamma + 1)) * self.exit.Tt
            self.exit.V = numpy.sqrt(gamma * R * self.exit.T)
            self.exit.P = self.exit.Pt * (1 / critical_NPR)
        else:
            # Unchoked
            self.exit.P = self.Pinf
            self.exit.T = self.exit.Tt * (self.exit.P/self.exit.Pt)**((gamma - 1) / gamma)
            self.exit.V = numpy.sqrt(2 * cp * (self.exit.Tt - self.exit.T))
            self.exit.M = self.exit.V / numpy.sqrt(gamma * R * self.exit.T)
        self.exit.set_statics(self.exit.M)

    # Cycle method
    def solve_CD(self, upstream):
        # This function forces the throat to always be choked (CHECK THIS)
        self.inlet = upstream
        gamma = self.inlet.gamma
        cp = self.inlet.cp
        R = self.inlet.R
        self.NPR = self.inlet.Pt / self.Pinf
        self.exit = copy.deepcopy(self.inlet)
        self.throat = copy.deepcopy(self.inlet)
        self.throat.idx = 8
        self.exit.idx = 9
        # Throat Calculations
        if self.inlet.M < 1:
            self.throat.M = 1
            self.throat.Tt = self.inlet.Tt
            self.throat.ht = self.throat.Tt * cp
            self.throat.Pt = self.inlet.Pt
            self.throat.mdot = self.inlet.mdot
            self.throat.set_statics(self.throat.M)
        # Exit Plane Calculations
        self.exit.P = self.Pinf; # Assuming perfectly expanded flow
        self.exit.M = numpy.sqrt((2/(gamma - 1)) * ((self.exit.Pt/self.exit.P)^((gamma-1)/gamma) - 1)); # Isentropic Relation
        self.exit.set_statics(self.exit.M)

    # Cycle method for performing preliminary component design ("specification" portion of the component design parameters)
    # Only call this after performing cycle analysis, NOT component design
    def get_specification(self):
        specification = {
            "W": self.inlet.W,
            "Pt": self.inlet.Pt / 1000,
            "Tt": self.inlet.Tt,
            "FAR": self.inlet.FAR,
            "M": self.inlet.M
        }
        return specification

    def design_component(self, component_parameters):
        self.component_parameters = component_parameters
        self.type = self.component_parameters["type"]
        W_in = self.component_parameters["inlet station"]["W"]
        Pt_in = self.component_parameters["inlet station"]["Pt"] * 1000
        Tt_in = self.component_parameters["inlet station"]["Tt"]
        FAR_in = self.component_parameters["inlet station"]["FAR"]
        M_in = self.component_parameters["inlet station"]["M"]
        self.inlet = Station(W_in, Tt_in, Pt_in, FAR=FAR_in, M=M_in)
        if not hasattr(self, "geometry"):
            self.geometry = self.component_parameters["geometry"]
        if not hasattr(self, "Pinf"):
            self.Pinf = self.component_parameters["Pinf"] * 1000
        match self.geometry:
            case "C": 
                self.solve_converging(self.inlet)
                self.alpha = numpy.deg2rad(self.component_parameters["half angle"])
                self.r_inlet = numpy.sqrt(self.inlet.area / numpy.pi)
                self.r_exit = numpy.sqrt(self.exit.area / numpy.pi)
                self.length = (self.r_inlet - self.r_exit) / numpy.tan(self.alpha)
                match self.type:
                    case "conical":
                        self.CA = (1 + numpy.cos(self.alpha)) / 2
                    case "rectangular":
                        self.CA = numpy.sin(self.alpha) / self.alpha
            case "CD": 
                self.solve_CD(self.inlet)
                self.C_alpha = numpy.deg2rad(self.component_parameters["converging half angle"])
                self.CD_alpha = numpy.deg2rad(self.component_parameters["diverging half angle"])
                self.r_inlet = numpy.sqrt(self.inlet.area / numpy.pi)
                self.r_th = numpy.sqrt(self.throat.area / numpy.pi)
                self.r_exit = numpy.sqrt(self.exit.area / numpy.pi)
                self.converging_length = (self.r_inlet - self.r_th) / numpy.tan(self.C_alpha)
                self.diverging_length = (self.r_exit - self.r_th) / numpy.tan(self.CD_alpha)
                self.total_length = self.converging_length + self.diverging_length
                match self.type:
                    case "conical":
                        self.CA = (1 + numpy.cos(self.CD_alpha)) / 2
                    case "rectangular":
                        self.CA = numpy.sin(self.CD_alpha) / self.CD_alpha


class Recuperator:
    def __init__(self, cold_inlet:Station, cycle_parameters, component_parameters=None):
        self.cyce_parameters = cycle_parameters
        self.engine = cycle_parameters["engine"]
        self.Mexit_cold = cycle_parameters["cold exit M"]
        self.Mexit_hot = cycle_parameters["hot exit M"]
        self.pi_cold = cycle_parameters["pi_cold"]
        self.pi_hot = cycle_parameters["pi_hot"]
        self.delta_ht = cycle_parameters["delta_ht"] * 10**3 # kJ/kg to J/kg

        self.cold_inlet = copy.deepcopy(cold_inlet)
        self.cold_inlet.idx = "3.06"
        self.cold_exit = copy.deepcopy(self.cold_inlet)
        self.cold_exit.idx = "3.07"
        self.cold_exit.Tt = self.cold_exit.T_from_H(self.cold_inlet.ht+self.delta_ht, self.cold_exit.FAR, 1800, self.cold_inlet.Tt)
        self.cold_exit.Pt *= self.pi_cold
        self.cold_exit.M = self.Mexit_cold
        self.cold_exit.set_statics(self.cold_exit.M)
    
        # COMPONENT DESIGN
        if component_parameters != None: pass

    def pass_hot_stream(self, hot_upstream):
        self.hot_inlet = hot_upstream
        self.hot_inlet.idx = "6.07"
        self.hot_exit = copy.deepcopy(self.hot_inlet)
        self.hot_exit.idx = "6.08"
        self.hot_exit.Tt = self.cold_exit.T_from_H(self.hot_inlet.ht-self.delta_ht, self.hot_exit.FAR, self.hot_inlet.Tt, 0)
        self.hot_exit.Pt *= self.pi_hot
        self.hot_exit.M = self.Mexit_hot
        self.hot_exit.set_statics(self.hot_exit.M)

        return self.hot_exit


class Bleed:
    def __init__(self, cooling, packing):
        self.cooling = cooling
        self.packing = packing


def format_axes(ax, title, ylabel):
    ax.set_title(title, fontsize=14, weight='bold')
    ax.set_ylabel(ylabel, fontsize=12)
    ax.legend(frameon=False, loc="best")
    # Grid control
    ax.minorticks_on()
    ax.grid(True, which='major', linestyle='--', linewidth=0.7, alpha=0.7)
    ax.grid(True, which='minor', linestyle=':', linewidth=0.5, alpha=0.5)
    # Clean spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


'''
Essentially a turbojet core that will serve as a parent class to other architectures like turbofan, turboshaft, ramjet, etc.
Handles a recuperator and afterburner as well, but the default is just a core
'''
class Engine:
    def __init__(self, engine_parameters):
        self.spools = engine_parameters["spools"]
        self.LHV = engine_parameters["LHV"]
        self.Minf = engine_parameters["Minf"]
        self.altitude = engine_parameters["altitude"]
        self.ambient = Ambient(self.altitude, Minf=self.Minf)

        # Handle user errors in engine parameters
        if "mdot" in engine_parameters and "fan diameter" in engine_parameters:
            raise ValueError("Error. Cannot have two sizing parameters. Must choose either inlet mass flow or fan face diameter.")
        else:
            if "mdot" in engine_parameters: self.W = engine_parameters["mdot"] # kg/sec
            elif "fan diameter" in engine_parameters: self.fan_diameter = engine_parameters["front diameter"] # inches
        if "TET" in engine_parameters and "mdotf" in engine_parameters:
            raise ValueError("Error. Cannot have two burner parameters. Must choose either TET or fuel mass flow rate.")
        else:
            if "TET" in engine_parameters: self.TET = engine_parameters["TET"] # K
            if "mdotf" in engine_parameters: self.Wf = engine_parameters["mdotf"] # kg/sec

    @property
    def components(self): return self._components

    @components.setter
    def components(self, value): self._components = value

    def set_components(self, parameters):
        self.parameters = parameters

        # Include the engine object in each set of component parameters
        for component_name in parameters:
            if component_name != "engine":
                parameters[component_name]["engine"] = self

        # Handle spool counts
        match self.spools:
            case 1:
                self.PR = parameters["engine"]["PR"]
                parameters["compressor"]["PR"] = self.PR
                self.components = [
                    inlet := Inlet(parameters["intake"]),
                    compressor := Compressor(inlet.exit, parameters["compressor"]),
                    burner := Burner(compressor.exit, parameters["burner"]),
                    turbine := Turbine(burner.exit, compressor, parameters["turbine"]),
                    exhaust := Nozzle(turbine.exit, parameters["nozzle"])
                ]
                '''
                self.components = [inlet := Inlet(parameters["intake"])]
                parameters["compressor"]["upstream"] = inlet.exit
                parameters["burner"]["upstream"] = compressor.exit
                parameters["turbine"]["upstream"] = burner.exit
                parameters["turbine"]["compressor"] = compressor
                parameters["nozzle"]["upstream"] = turbine.exit
                self.components.append(compressor := Compressor(parameters["compressor"]))
                self.components.append(burner := Burner(parameters["burner"]))
                self.components.append(turbine := Turbine(parameters["turbine"]))
                self.components.append(exhaust := Nozzle(turbine.exit, parameters["nozzle"]))
                self.compressor = compressor
                self.turbine = turbine
                '''
                self.compressor = compressor
                self.turbine = turbine
            case 2:
                self.OPR = parameters["engine"]["OPR"]
                self.HPR = parameters["engine"]["HPR"]
                self.LPR = self.OPR / self.HPR
                parameters["hpc"]["PR"] = self.HPR
                parameters["lpc"]["PR"] = self.LPR

                self.components = [
                    inlet := Inlet(parameters["intake"]),
                    lpc := Compressor(inlet.exit, parameters["lpc"]),
                    hpc := Compressor(lpc.exit, parameters["hpc"]),
                    burner := Burner(hpc.exit, parameters["burner"]),
                    hpt := Turbine(burner.exit, hpc, parameters["hpt"]),
                    lpt := Turbine(hpt.exit, lpc, parameters["lpt"]),
                    exhaust := Nozzle(lpt.exit, parameters["nozzle"])
                ]
                self.lpc = lpc
                self.hpc = hpc
                self.hpt = hpt
                self.lpt = lpt
            case 3:
                self.OPR = parameters["engine"]["OPR"]
                self.HPR = parameters["engine"]["HPR"]
                self.IPR = parameters["engine"]["IPR"]
                self.LPR = self.OPR / self.HPR / self.IPR
                parameters["hpc"]["PR"] = self.HPR
                parameters["ipc"]["PR"] = self.IPR
                parameters["lpc"]["PR"] = self.LPR

                self.components = [
                    inlet := Inlet(parameters["intake"]),
                    lpc := Compressor(inlet.exit, parameters["lpc"]),
                    ipc := Compressor(lpc.exit, parameters["hpc"]),
                    hpc := Compressor(ipc.exit, parameters["hpc"]),
                    burner := Burner(hpc.exit, parameters["burner"]),
                    hpt := Turbine(burner.exit, hpc, parameters["hpt"]),
                    ipt := Turbine(hpt.exit, ipc, parameters["ipt"]),
                    lpt := Turbine(ipt.exit, lpc, parameters["lpt"]),
                    exhaust := Nozzle(lpt.exit, parameters["nozzle"])
                ]
                self.lpc = lpc
                self.ipc = ipc
                self.hpc = hpc
                self.hpt = hpt
                self.ipt = ipt
                self.lpt = lpt
            case _:
                raise ValueError("Number of spools must be either 1, 2, or 3.")

        self.inlet = inlet
        self.burner = burner
        self.exhaust = exhaust
        # Check for a recuperator in the user input parameters
        if "recuperator" in parameters: self.add_recuperator(parameters["recuperator"])
        # Check for an afterburner in the user input parameters
        if "afterburner" in parameters: self.add_afterburner(parameters["afterburner"])


    # Add an afterburner to the engine
    def add_afterburner(self, parameters):
        insert_idx = len(self.components) - 1
        inlet_station = self.components[insert_idx - 1].exit
        afterburner = Afterburner(inlet_station, parameters)
        self.afterburner = afterburner
        self.components[insert_idx].solve_exit(afterburner.exit)
        self.components.insert(insert_idx, afterburner) 


    # Add a recuperator to the engine
    def add_recuperator(self, parameters):
        for idx, component in enumerate(self.components):
            if isinstance(component, Burner): 
                insert_idx = idx
                recuperator = Recuperator(self.components[idx - 1].exit, parameters)
                component.solve_exit(recuperator.cold_exit)
            elif isinstance(component, Turbine):
                component.solve_exit(self.components[idx-1].exit)
            elif not isinstance(component, Turbine) and isinstance(self.components[idx - 1], Turbine):
                component.solve_exit(recuperator.pass_hot_stream(self.components[idx-1].exit))
            elif idx > 2:
                component.solve_exit(self.components[idx - 1].exit)
        self.components.insert(insert_idx, recuperator)
        self.recuperator = recuperator


    # Retrieve flow properties at every station
    def get_station_data(self): 
        raw_data = list()
        # Handle Recuperator station data (afterburner is treated like the other components)
        if "recuperator" in self.parameters:
            for component in self.components:
                if isinstance(component, Inlet):
                    raw_data.append(component.freestream.get_properties())
                    raw_data.append(component.inlet.get_properties())
                    raw_data.append(component.exit.get_properties())
                elif isinstance(component, Recuperator):
                    raw_data.append(component.cold_exit.get_properties())
                elif isinstance(component, Turbine):
                    raw_data.append(component.exit.get_properties())
                    raw_data.append(self.recuperator.pass_hot_stream(component.exit).get_properties())
                else:
                    raw_data.append(component.exit.get_properties())
        else:
            for component in self.components:
                if isinstance(component, Inlet):
                    raw_data.append(component.freestream.get_properties())
                    raw_data.append(component.inlet.get_properties())
                    raw_data.append(component.exit.get_properties())
                else:
                    raw_data.append(component.exit.get_properties())
                    
        rounded_data = numpy.round(numpy.array(raw_data, dtype=float), 3).tolist()
        station_data = pandas.DataFrame(rounded_data, columns=Station.column_names)
        
        return station_data

    # Plot the temperatures and pressures throughout the whole engine
    def plot_thermo(self):
        Tt = self.get_station_data()["Tt [K]"]
        Ts = self.get_station_data()["Ts [K]"]
        Pt = self.get_station_data()["Pt [kPa]"]
        Ps = self.get_station_data()["Ps [kPa]"]
        stations = self.get_station_data()["Station"]
        marker = "o"
        size = 4
        color = "white"
        edge = 1.5
        width = 1
        plt.close("all")
        plt.style.use('seaborn-v0_8-whitegrid')
        fig, axes = plt.subplots(2, 1, figsize=(7, 7), dpi=120, sharex=True)
        # Temperatures
        axes[0].plot(stations, Tt, marker=marker, markersize=size, markerfacecolor=color, markeredgewidth=edge, linewidth=width, label="Total Temperature")
        axes[0].plot(stations, Ts, marker=marker, markersize=size, markerfacecolor=color, markeredgewidth=edge, linewidth=width, label="Static Temperature")
        format_axes(axes[0], "Temperature", "Temperature [K]")
        # Pressures
        axes[1].plot(stations, Pt, marker=marker, markersize=size, markerfacecolor=color, markeredgewidth=edge, linewidth=width, label="Total Pressure")
        axes[1].plot(stations, Ps, marker=marker, markersize=size, markerfacecolor=color, markeredgewidth=edge, linewidth=width, label="Static Pressure")
        format_axes(axes[1], "Pressure", "Pressure [kPa]")
        axes[1].set_xlabel("Station", fontsize=11)
        plt.tight_layout()
        plt.show()

    # Display the full engine performance
    def get_performance(self):
        """ Performance Parameters """
        station_data = self.get_station_data()
        W = station_data["W [kg/sec]"].values
        V = station_data["V [m/sec]"].values
        P = station_data["Ps [kPa]"].values
        A = station_data["Area [m^2]"].values
        FAR = station_data["FAR"].values
        Wf_in = self.burner.exit.Wf
        FARexit = self.exhaust.exit.FAR
        pressure_thrust = A[-1] * ((P[-1]*1000) - self.ambient.P)
        pressure_power = pressure_thrust * V[-1]
        # Specific Thrust
        T_ma = V[-1]*(1 + FAR[-1]) - V[0] + pressure_thrust/W[0] 
        # Thrust Specific Fuel Consumption (TSFC)
        TSFC = (FAR[-1] / T_ma) * 10**6 
        # Thrust
        thrust = T_ma * W[0] 
        # Propulsive Efficiency
        eta_p = thrust*V[0] / (0.5*W[-1]*V[-1]**2 - 0.5*W[0]*V[0]**2 + pressure_power) 
        # Thermal Efficiency
        eta_th = (0.5*W[-1]*V[-1]**2 - 0.5*W[0]*V[0]**2 + pressure_power) / (Wf_in*self.LHV) 
        # Overall Efficiency
        eta_o = eta_p * eta_th 
        performance = pandas.DataFrame({
                                         "Specific Thrust [m/sec]": T_ma,
                                         "TSFC [g/(sec*kN)]": TSFC,
                                         "Thrust [kN]": thrust/1000,
                                         "Propulsive Efficiency": eta_p,
                                         "Thermal Efficiency": eta_th,
                                         "Overall Efficiency": eta_o
                                        }, index=[0])
        return station_data, performance

    def get_specifications(self):
        specifications = dict()
        for component in self.components:
            if isinstance(component, Compressor):
                specifications["compressor"] = component.get_specification()
            elif isinstance(component, Turbine):
                specifications["turbine"] = component.get_specification()
            elif isinstance(component, Nozzle):
                specifications["nozzle"] = component.get_specification()
        return specifications

    def optimize(self, performance_parameter): pass
    def sensitivity_study(self): pass
    def off_design(self): pass


''' handles turbofan (mixed/unmixed) and variable bypass ramjet '''
class BypassEngine(Engine):
    def __init__(self, engine_parameters):
        self.spools = engine_parameters["spools"]
        self.B = engine_parameters["bypass ratio"]
        self.mixed = engine_parameters["mixed"]
        self.is_fan = engine_parameters["fan"]
        self.TET = engine_parameters["TET"]
        self.LHV = engine_parameters["LHV"]
        self.Minf = engine_parameters["Minf"]
        self.altitude = engine_parameters["altitude"]
        self.ambient = Ambient(self.altitude, Minf=self.Minf)

        # Handle user errors in engine parameters
        if "mdot" in engine_parameters and "fan diameter" in engine_parameters:
            raise ValueError("Error. Cannot have two sizing parameters. Must choose either inlet mass flow or fan face diameter.")
        else:
            if "mdot" in engine_parameters: self.W = engine_parameters["mdot"] # kg/sec
            elif "fan diameter" in engine_parameters: self.fan_diameter = engine_parameters["fan diameter"] # inches


    def set_components(self, parameters):
        self.parameters = parameters

        # Include the engine object in each set of component parameters
        for component_name in parameters:
            if component_name != "engine":
                parameters[component_name]["engine"] = self

        # Handle spool counts
        match self.spools:
            case 1:
                self.PR = parameters["engine"]["PR"]
                parameters["compressor"]["PR"] = self.PR

                # Check for user errors
                if self.is_fan: raise ValueError("Cannot have a fan with only 1 spool. Must be at least two spools.")
                if self.mixed == False: raise ValueError("Must have a mixed exhaust with only 1 spool.")
                self.components = [
                    inlet := Inlet(parameters["intake"]),
                    compressor := Compressor(inlet.root_exit, parameters["compressor"]),
                    burner := Burner(compressor.exit, parameters["burner"]),
                    turbine := Turbine(burner.exit, compressor, parameters["turbine"]),
                    mixer := Mixer(turbine.exit, inlet.bypass_exit),
                    exhaust := Nozzle(mixer.exit, parameters["nozzle"])
                ]
                self.compressor = compressor
                self.turbine = turbine
            case 2:
                # Turbofan or Turbojet
                if self.is_fan: 
                    parameters["lpc"]["fan"] = True
                    self.OPR = parameters["engine"]["OPR"]
                    self.LPR = parameters["lpc"]["root PR"]
                    self.HPR = self.OPR / self.LPR
                    parameters["hpc"]["PR"] = self.HPR

                    # Mixed Exhaust
                    match self.mixed:
                        case True: 
                            self.components = [
                                inlet := Inlet(parameters["intake"]),
                                lpc := Compressor(inlet.bypass_exit, parameters["lpc"], root_upstream=inlet.root_exit),
                                hpc := Compressor(lpc.root_exit, parameters["hpc"]),
                                burner := Burner(hpc.exit, parameters["burner"]),
                                hpt := Turbine(burner.exit, hpc, parameters["hpt"]),
                                lpt := Turbine(hpt.exit, lpc, parameters["lpt"]),
                                mixer := Mixer(lpt.exit, lpc.tip_exit),
                                exhaust := Nozzle(mixer.exit, parameters["nozzle"])
                            ]
                            self.mixer = mixer
                        case False: 
                            self.components = [
                                inlet := Inlet(parameters["intake"]),
                                lpc := Compressor(inlet.bypass_exit, parameters["lpc"], root_upstream=inlet.root_exit),
                                cold_nozzle := Nozzle(lpc.tip_exit, parameters["cold nozzle"]),
                                hpc := Compressor(lpc.exit, parameters["hpc"]),
                                burner := Burner(hpc.exit, parameters["burner"]),
                                hpt := Turbine(burner.exit, hpc, parameters["hpt"]),
                                lpt := Turbine(hpt.exit, lpc, parameters["lpt"]),
                                exhaust := Nozzle(lpt.exit, parameters["nozzle"])
                            ]
                            self.cold_nozzle = cold_nozzle
                else: 
                    parameters["lpc"]["fan"] = False
                    self.OPR = parameters["engine"]["OPR"]
                    self.LPR = parameters["engine"]["LPR"]
                    self.HPR = self.OPR / self.LPR
                    parameters["hpc"]["PR"] = self.HPR
                    parameters["lpc"]["PR"] = self.LPR

                    # Mixed Exhaust
                    if self.mixed == False: raise ValueError("Must have a mixed exhaust if there's no fan.")

                    self.components = [
                        inlet := Inlet(parameters["intake"]),
                        lpc := Compressor(inlet.root_exit, parameters["lpc"]),
                        hpc := Compressor(lpc.exit, parameters["hpc"]),
                        burner := Burner(hpc.exit, parameters["burner"]),
                        hpt := Turbine(burner.exit, hpc, parameters["hpt"]),
                        lpt := Turbine(hpt.exit, lpc, parameters["lpt"]),
                        mixer := Mixer(lpt.exit, inlet.bypass_exit),
                        exhaust := Nozzle(mixer.exit, parameters["nozzle"])
                    ]
                    self.mixer = mixer
                self.lpc = lpc
                self.hpc = hpc
                self.hpt = hpt
                self.lpt = lpt
            case 3:
                # Check if the engine is a turbofan
                if self.is_fan:
                    parameters["lpc"]["fan"] = True
                    self.OPR = parameters["engine"]["OPR"]
                    self.LPR = parameters["lpc"]["root PR"]
                    self.HPR = parameters["engine"]["HPR"]
                    self.IPR = self.OPR / self.HPR / self.LPR
                    parameters["hpc"]["PR"] = self.HPR
                    parameters["ipc"]["PR"] = self.IPR
                    # Mixed Exhaust
                    match self.mixed:
                        case True: 
                            self.components = [
                                inlet := Inlet(parameters["intake"]),
                                lpc := Compressor(inlet.exit, parameters["lpc"]),
                                ipc := Compressor(lpc.exit, parameters["hpc"]),
                                hpc := Compressor(ipc.exit, parameters["hpc"]),
                                burner := Burner(hpc.exit, parameters["burner"]),
                                hpt := Turbine(burner.exit, hpc, parameters["hpt"]),
                                ipt := Turbine(hpt.exit, ipc, parameters["ipt"]),
                                lpt := Turbine(ipt.exit, lpc, parameters["lpt"]),
                                mixer := Mixer(lpt.exit, lpc.tip_exit),
                                exhaust := Nozzle(mixer.exit, parameters["nozzle"])
                            ]
                            self.mixer = mixer
                        case False:
                            self.components = [
                                inlet := Inlet(parameters["intake"]),
                                lpc := Compressor(inlet.exit, parameters["lpc"]),
                                cold_nozzle := Nozzle(lpc.tip_exit, parameters["cold nozzle"]),
                                ipc := Compressor(lpc.exit, parameters["hpc"]),
                                hpc := Compressor(ipc.exit, parameters["hpc"]),
                                burner := Burner(hpc.exit, parameters["burner"]),
                                hpt := Turbine(burner.exit, hpc, parameters["hpt"]),
                                ipt := Turbine(hpt.exit, ipc, parameters["ipt"]),
                                lpt := Turbine(ipt.exit, lpc, parameters["lpt"]),
                                exhaust := Nozzle(lpt.exit, parameters["nozzle"])
                            ]
                            self.cold_nozzle = cold_nozzle
                else:
                    parameters["lpc"]["fan"] = False
                    self.OPR = parameters["engine"]["OPR"]
                    self.HPR = parameters["engine"]["HPR"]
                    self.IPR = parameters["engine"]["IPR"]
                    self.LPR = self.OPR / self.HPR / self.IPR
                    parameters["hpc"]["PR"] = self.HPR
                    parameters["ipc"]["PR"] = self.IPR
                    parameters["lpc"]["PR"] = self.LPR
                    # Mixed Exhaust
                    if self.mixed == False: raise ValueError("Must have a mixed exhaust if there's no fan.")
                    self.components = [
                        inlet := Inlet(parameters["intake"]),
                        lpc := Compressor(inlet.root_exit, parameters["lpc"]),
                        ipc := Compressor(lpc.exit, parameters["hpc"]),
                        hpc := Compressor(ipc.exit, parameters["hpc"]),
                        burner := Burner(hpc.exit, parameters["burner"]),
                        hpt := Turbine(burner.exit, hpc, parameters["hpt"]),
                        ipt := Turbine(hpt.exit, ipc, parameters["ipt"]),
                        lpt := Turbine(ipt.exit, lpc, parameters["lpt"]),
                        mixer := Mixer(lpt.exit, inlet.bypass_exit),
                        exhaust := Nozzle(mixer.exit, parameters["nozzle"])
                    ]
                    self.mixer = mixer
                self.lpc = lpc
                self.ipc = ipc
                self.hpc = hpc
                self.hpt = hpt
                self.ipt = ipt
                self.lpt = lpt
            case _:
                raise ValueError("Number of spools must be either 1, 2, or 3.")

        self.inlet = inlet
        self.burner = burner
        self.exhaust = exhaust
        # Check for an afterburner in the user input parameters
        if "afterburner" in parameters: self.add_afterburner(parameters["afterburner"])

    # Retrieve flow properties at every station
    def get_station_data(self): 
        raw_data = list()
        
        for component in self.components:
            if isinstance(component, Inlet):
                raw_data.append(component.freestream.get_properties())
                raw_data.append(component.inlet.get_properties())
                raw_data.append(component.bypass_exit.get_properties())
                raw_data.append(component.root_exit.get_properties())
            elif isinstance(component, Compressor):
                if hasattr(component, "is_fan"):
                    raw_data.append(component.root_exit.get_properties())
                    raw_data.insert(-2, component.tip_exit.get_properties())
                else:
                    raw_data.append(component.exit.get_properties())
            else:
                raw_data.append(component.exit.get_properties())
                    
        rounded_data = numpy.round(numpy.array(raw_data, dtype=float), 3).tolist()
        station_data = pandas.DataFrame(rounded_data, columns=Station.column_names)
        return station_data


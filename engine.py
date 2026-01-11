"""
class for cycle analysis
"""

from lookup_tables import *
import numpy

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

    def __init__(self, W, Tt, Pt, ht, name:dict=None):
        self.W = W
        self.Tt = Tt
        self.Pt = Pt
        self.ht = ht
        self.Wc = self.get_Wc()
        self.name = name

        # Statics
        self.M, self.T, self.P, self.V, self.rho, self.A = 0

        # Gas Model
        self.FAR = 0
        self.cp = self.get_cp(self.Tt, self.FAR)
        self.R = self.get_R(self.FAR)
        self.gamma = self.cp / (self.cp - self.R)
        
    def set_statics(self, M):
        [_, Tt_T, Pt_P, _, _] = isentropic(self.M, lookup_key="M")
        self.T = (1 / Tt_T) * self.Tt
        self.P = (1 / Pt_P) * self.Pt
        self.V = self.M * numpy.sqrt(self.gamma * self.R * self.T)
        self.rho = self.P / (self.R * self.T)
        self.A = self.W / (self.rho * self.V)

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
    
    def get_cp(self, T, FAR):
        """ Formula 3.24 """
        TZ = T/1000
        cp = (
              self.A["A0"] + self.A["A1"] * TZ + (self.A["A2"] * TZ**2) + (self.A["A3"] * TZ**3) +
              (self.A["A4"] * TZ**4) + (self.A["A5"] * TZ**5) + (self.A["A6"] * TZ**6) +
              (self.A["A7"] * TZ**7) + (self.A["A8"] * TZ**8) + (self.B["B0"] + self.B["B1"] * TZ +
              (self.B["B2"] * TZ**2) + (self.B["B3"] * TZ**3) + (self.B["B4"] * TZ**4) +
              (self.B["B5"] * TZ**5) + (self.B["B6"] * TZ**6) + (self.B["B7"] * TZ**7)) * (FAR / (1 + FAR))
             )
        return cp
    
    def get_R(self, FAR):
        """ Formula 3.22 for Kerosene """
        return 287.05 - 0.0099 * FAR + 0.0000001 * FAR**2
    
    def get_Wc(self):
        Wc = self.W * math.sqrt(self.Tt / self.Tref) / (self.Pt / self.Tref)
        return Wc

class Inlet:
    def __init__(self, upstream:Station, Pt_recovery):
        self.exit_W = upstream.W
        self.exit_Pt = upstream.Pt * Pt_recovery
        self.exit_Tt = upstream.Tt
        self.exit_ht = upstream.get_ht(self.exit_Tt, 0)
        self.exit = Station(self.exit_W, self.exit_Tt, self.exit_Pt, self.exit_ht)

class Compressor:
    def __init__(self, upstream:Station, PR):
        self.exit_W = upstream.W
        self.exit_Pt = upstream.Pt * PR
        self.exit_Tt = upstream.Tt
        self.exit_ht = upstream.get_ht(self.exit_Tt, 0)
        self.exit = Station(self.exit_W, self.exit_Tt, self.exit_Pt, self.exit_ht)

class Burner:
    def __init__(self, upstream:Station, TIT, LHV, pressure_loss):
        self.eta_b = 1
        self.exit_FAR = self.get_FAR(TIT, self.inlet.Tt, LHV)
        self.exit_W = upstream.W * (1 + self.exit_FAR)
        self.exit_Pt = upstream.Pt * (1 - pressure_loss)
        self.exit_Tt = TIT
        self.exit_ht = upstream.get_ht(self.exit_Tt, self.exit_FAR)
        self.exit = Station(self.exit_W, self.exit_Tt, self.exit_Pt, self.exit_ht)
    
    def get_FAR(self, T2, T1, LHV):
        """ Iterating Fuel-to-Air ratio until the exit total temperature matches the TIT """
        # Initialize variables
        FARnew = 0.02
        FAR = -1 # Start with a value that ensures the loop condition is met
        h1 = self.outlet.get_ht(T1, 0)
        
        # Perform iterations
        while (abs(FAR - FARnew) / FARnew) > 0.00001:
            FAR = FARnew
            h2 = self.outlet.get_ht(T2, FAR)
            FARnew = (h2 - h1) / (LHV * self.eta_b)

        return FARnew

class Turbine:
    def __init__(self, upstream:Station):
        pass
        '''
        self.exit_W = upstream.W * (1 + self.exit_FAR)
        self.exit_Pt = upstream.Pt * (1 - pressure_loss)
        self.exit_Tt = self.TfromH(self.exit_ht, self.FAR, TIT, 0)
        self.exit_ht = upstream.get_ht(self.exit_Tt, self.exit_FAR)
        self.exit = Station(self.exit_W, self.exit_Tt, self.exit_Pt, self.exit_ht)
        '''

    def TfromH(self, h, FAR, Thi, Tlo):
        """ Finding temperature from enthalpy polynomial using the Bisection Method """
        # Initial calculation of hmid to prevent uninitialized use
        Tmid = (Thi + Tlo) / 2
        hmid = self.get_ht(Tmid, FAR)
        
        # Perform bisection
        while (abs(hmid - h) / h) > 0.001:
            Tmid = (Thi + Tlo) / 2
            hmid = self.get_ht(Tmid, FAR)
            
            if hmid < h:
                Tlo = Tmid
            elif hmid > h:
                Thi = Tmid
            
            iterations = iterations + 1
            TfromH = Tmid
        return TfromH

class Nozzle:
    def __init__(self, upstream):
        self.exit_W = upstream.W
        self.exit_Pt = upstream.Pt
        self.exit_Tt = upstream.Tt
        self.exit_ht = upstream.get_ht(self.exit_Tt, self.exit_FAR)
        self.exit = Station(self.exit_W, self.exit_Tt, self.exit_Pt, self.exit_ht)

        def set_statics(self,  Pa):
            gamma = self.exit.gamma
            critical_NPR = (1 + ((gamma-1) / 2))**(gamma / (gamma-1))
            NPR = upstream.Pt / Pa
            if NPR > critical_NPR: self.exit.M = 1

class Diffuser:
    def __init__(self, upstream):
        self.exit_W = upstream.W
        self.exit_Pt = upstream.Pt
        self.exit_Tt = upstream.Tt
        self.exit_ht = upstream.get_ht(self.exit_Tt, self.exit_FAR)
        self.exit = Station(self.exit_W, self.exit_Tt, self.exit_Pt, self.exit_ht)

class Engine:
    def __init__(self):
        self.OPR = 20
        self.TIT = 1600
        self.BPR = 4
        self.FPR = 2.1
        self.CPR = self.OPR / self.FPR
        self.LHV = 43.15*10**6
        self.inlet = Inlet()
        self.compressor = Compressor()
        self.burner = Burner()
        self.turbine =  Turbine()
        self.exhaust =  Nozzle()

    def validate_cycle(self):
        pass
    def optimize(self, perf_param):
        pass
    def sensitivity_study(self):
        pass
    def off_design(self):
        pass
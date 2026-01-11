import matplotlib.pyplot as plt
import pandas
import numpy
import math
import copy
import os

# Change the current working directory to the file location
filepath = os.path.abspath(__file__)
directory = os.path.dirname(filepath)
os.chdir(directory)

from atmosphere import atmosphere

def Cdo_correlation(Cdo):
    if Cdo > 0.8:
        return 0.02 + ((Cdo - 0.8)*0.286) + (0.02*(0.2**2))
    else:
        return 0.02 + (0.02*(0.2**2))

def excess_power():
    thrust_loading = 0.45
    wing_loading = 45 * 47.88
    n = 1
    beta = 1
    num_points = 100

    M = numpy.linspace(0.000001, 1, num_points)
    altitude = numpy.linspace(0, 60000, num_points) * 0.3048
    ambient = list(map(atmosphere, list(altitude)))
    temperature, pressure, rho, a = zip(*ambient)
    temperature, pressure = numpy.meshgrid(temperature, pressure)
    rho, a = numpy.meshgrid(rho, a)
    temperature = temperature.T
    rho = rho.T
    M, altitude = numpy.meshgrid(M, altitude)

    Cdo = numpy.vectorize(Cdo_correlation)(copy.deepcopy(M))
    Pt_Ps = (1 + ((gamma - 1)/2) * M**2)**(gamma / (gamma - 1))
    alpha = Pt_Ps * pressure / Pref
    V = M * a
    q = (1/2)* rho * V**2
    Ps = V * ((alpha/beta)*thrust_loading - K1*wing_loading*(n**2)*beta/q - K2*n - (Cdo/(beta*wing_loading/q)))
    TSFC = (0.9 + 0.3 * M) * numpy.vectorize(math.sqrt)(temperature/Tref) 
    fs = Ps / alpha / TSFC
    fig, ax = plt.subplots()
    contour = ax.contourf(M, altitude/0.3048, fs)
    cbar = fig.colorbar(contour, ax=ax)
    plt.show()

def get_thrust_loading(M, altitude, beta, wing_loading, n, K1, K2, CDR, VS, acceleration):
    temperature, pressure, rho, a = atmosphere(altitude)
    Cdo = numpy.vectorize(Cdo_correlation)(copy.deepcopy(M))
    Pt_Ps = (1 + ((gamma - 1)/2) * M**2)**(gamma / (gamma - 1))
    alpha = Pt_Ps * pressure / Pref
    V = M * a
    q = (1/2)* rho * V**2
    thrust_loading = (beta/alpha) * ((q/beta/wing_loading) * (K1 * (n*beta*wing_loading/q)**2 + K2*(n*beta*wing_loading/q) + Cdo + CDR) + (VS/V) * (acceleration/g))

    return thrust_loading

def constraint(data):
    for idx, row in data.iterrows():
        M = row["M"]
        altitude = row["altitude"] * 0.3048 # m
        beta = row["beta"]
        n = row["n"]
        VS = row["VS"]
        sample_wing_loadings = numpy.linspace(20, 120, 200) * 47.88
        thrust_loadings = numpy.vectorize(get_thrust_loading)(M, altitude, beta, sample_wing_loadings, n, K1, K2, CDR, VS, acceleration)
        plt.scatter(sample_wing_loadings, thrust_loadings)
    plt.show()

def mission(data):
    for idx, row in data.iterrows():
        M = row["M"]
        altitude = row["altitude"] * 0.3048 # m
        beta = row["beta"]
        n = row["n"]
        VS = row["VS"]
        temperature, pressure, rho, a = atmosphere(altitude)
        Pt_Ps = (1 + ((gamma - 1)/2) * M**2)**(gamma / (gamma - 1))
        alpha = Pt_Ps * pressure / Pref
        TSFC = (0.9 + 0.3 * M)*math.sqrt(temperature/Tref)
        thrust_loading = get_thrust_loading(M, altitude, beta, wing_loading, n, K1, K2, CDR, VS, acceleration)

# INPUTS
gamma = 1.4
R = 287
g = 9.8
Tref, Pref, rhoref, aref = atmosphere(0)
AR = 10
e = 0.8
K1 = 1/(math.pi*AR*e)+0.02
K2 = -0.008
CDR = 0
acceleration = 0
wing_loading = 45 * 47.88

data = pandas.read_excel("inputs.xlsx")

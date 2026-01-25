import math

# Only metric units
def atmosphere(height, unit = "metric"):
    #if unit == "imperial": height = height * 0.3048

    # assuming height is in m
    troposphere = 11000 # m
    stratosphere = 25000 # m
    if height < troposphere:
        temperature = 15.04 - 0.00649 * height
        pressure = 101.29 * ((temperature + 273.1) / 288.08)**5.256
    elif troposphere < height < stratosphere:
        temperature = -56.46
        pressure = 22.65 * math.exp(1.73 - 0.000157 * height)
    elif height > stratosphere:
        temperature = -131.21 + 0.00299 * height
        pressure = 2.488 * ((temperature + 273.1) / 216.6)**-11.388
    
    pressure = pressure * 10**3 # Pa
    temperature = temperature + 273.15 # K
    R = 287
    rho = pressure / (R * temperature)

    return [temperature, pressure, rho]


class Ambient:
    R = 287
    gamma = 1.4

    def __init__(self, Minf, altitude, alpha):
        self.T, self.P, self.rho = atmosphere(altitude)
        self.rhoinf = self.P / (self.R * self.T)
        self.Minf = Minf
        self.a = math.sqrt(self.gamma * self.R * self.T)
        self.Vinf = self.Minf * self.a
        self.alpha = alpha

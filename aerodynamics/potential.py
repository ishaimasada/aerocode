# Module for potential flow applications
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import numpy
import copy
import os
import sys

# Change the current working directory to the file location
filepath = os.path.abspath(__file__)
directory = os.path.dirname(filepath)
os.chdir(directory)

from curves import Point
from atmosphere import Ambient

''' classes for thin airfoil theory '''
''' Numerical Thin Airfoil Theory using Weissinger's Approximation '''
class Panel:
    def __init__(self, start_point:Point, end_point:Point):
        self.start_point = copy.deepcopy(start_point)
        self.end_point = copy.deepcopy(end_point)

    def get_position(self, parameter):
        v = self.end_point - self.start_point
        position = self.start_point + parameter * v
        return position
    
    @property
    def length(self):
        return numpy.sqrt((self.end_point.x_coord - self.start_point.x_coord)**2 + (self.end_point.y_coord - self.start_point.y_coord)**2 + (self.end_point.z_coord - self.start_point.z_coord)**2)
    
    @property
    def tangent_vector(self):
        return (self.end_point - self.start_point) / self.length

    @property
    def normal_vector(self):
        pass

    @property
    def phi(self):
        if self.end_point.x_coord - self.start_point.x_coord != 0: 
            phi = numpy.atan2((self.end_point.y_coord - self.start_point.y_coord), (self.end_point.x_coord - self.start_point.x_coord)) 
        else: phi = 0
        return phi
    
    @property
    def beta(self):
        return self.phi + (numpy.pi / 2)
    
class Thin_Airfoil:
    def __init__(self, designation:str, chord, num_panels, flap_angle=0, flap_start=None):
        self.designation = designation
        self.num_panels = num_panels
        self.chord = chord
        self.M = int(self.designation[0]) / 100
        self.P = int(self.designation[1]) / 10
        self.flap_angle = flap_angle
        self.flap_start = flap_start
        self.panels = self.generate_panels()

    def generate_panels(self):
        x_positions = numpy.linspace(0, 1, self.num_panels + 1)
        z_positions = numpy.array([self.get_camber(x)[0] for x in x_positions])
        panels = [Panel(Point(x_positions[i], z_positions[i], 0), Point(x_positions[i + 1], z_positions[i + 1], 0)) for i in range(self.num_panels)]

        common_length = 1 / (self.num_panels + 1)

        precision = 5
        error_threshold = 10**-precision
        panel_lengths = set(round(panel.length, precision) for panel in panels)
        while len(panel_lengths) > 1:
            for idx, panel in enumerate(panels):
                if idx > 0:
                    panel.start_point = copy.deepcopy(panels[idx - 1].end_point)
                if idx == len(panels) - 1:
                    continue
                length_error = abs(common_length - panel.length)
                while length_error > error_threshold:
                    #print(panel.get_length(), common_length, panel.start_point.x_coord, panel.end_point.x_coord)
                    if panel.length < common_length:
                        panel.end_point.set_x(panel.end_point.x_coord + error_threshold/2)
                    elif panel.length > common_length and panel.end_point.x_coord > panel.start_point.x_coord:
                        panel.end_point.set_x(panel.end_point.x_coord - error_threshold/2)
                    if panel.end_point.x_coord > 1: panel.end_point.set(1)
                    else: break
                    panel.end_point.set_y(self.get_camber(panel.end_point.x_coord)[0])
                    length_error = abs(common_length - panel.length)
            if panels[-1].length > common_length:
                common_length += abs(common_length - panels[-1].length) / self.num_panels
            elif panels[-1].length < common_length:
                common_length -= abs(common_length - panels[-1].length) / self.num_panels
            panel_lengths = set(round(panel.length, precision - 2) for panel in panels)

        # Each panel's vortex influences every panel's 3/4 flow tangency point
        self.quarter_points = numpy.tile(numpy.array([panel.get_position(0.25) for panel in panels]).reshape(self.num_panels, 1), (1, self.num_panels)) # vortex position vectors
        self.three_quarter_points = numpy.tile(numpy.array([panel.get_position(0.75) for panel in panels]).reshape(1, self.num_panels), (self.num_panels, 1)) # flow tangency position vectors

        return panels

    def get_camber(self, x_c):
        ''' Equations referenced from Airfoil Tools '''
        # input x/c instead of x-coordinates to get the chord-normalized z-coordinates
        if 0 <= x_c < self.P:
            z_c = (self.M/self.P**2) * (2*self.P*x_c - x_c**2)
            dz_dx = (2*self.M/self.P**2) * (self.P - x_c)
        elif self.flap_start == None and self.P <= x_c <= 1:
            z_c = (self.M/(1 - self.P)**2) * (1 - 2*self.P + 2*self.P*x_c - x_c**2)
            dz_dx = (2*self.M/(1 - self.P)**2) * (self.P - x_c)
        elif self.flap_start != None and self.P <= x_c <= self.flap_start:
            z_c = (self.M/(1 - self.P)**2) * (1 - 2*self.P + 2*self.P*x_c - x_c**2)
            dz_dx = (2*self.M/(1 - self.P)**2) * (self.P - x_c)
        elif self.flap_start != None and self.flap_start < x_c <= 1:
            flap_start_height = (self.M/(1 - self.P)**2) * (1 - 2*self.P + 2*self.P*self.flap_start - self.flap_start**2) # b
            dz_dx = numpy.tan(-self.flap_angle) # m
            z_c = dz_dx*(x_c - self.flap_start) + flap_start_height # mx + b
        else: return None
        
        return z_c, dz_dx
    
    @property
    def lift(self, ambient:Ambient):

        # Camber slopes at the points of tangency (three quarter point on panel)
        slopes = numpy.array([self.get_camber(position.x_coord)[1] for position in self.three_quarter_points[0]]) 
        
        # Connect every quarter-panel point to every three-quarter-panel point of every panel
        panel_connectors = numpy.vectorize(Panel)(self.quarter_points, self.three_quarter_points)
        r = numpy.vectorize(lambda panel: panel.length)(panel_connectors) # distances between control points and flow tangency points

        # Induced velocity with an absolute angle determined by a clockwise circulation
        theta = numpy.vectorize(lambda panel: panel.phi + 3*numpy.pi/2)(panel_connectors) # direction of induced velocity
        panel_lengths = numpy.tile(numpy.array([panel.length for panel in self.panels]), (self.num_panels, 1))

        # Apply Weissinger's Approximation in matrix form
        A = panel_lengths / (2*numpy.pi*r) * (numpy.sin(theta) - numpy.cos(theta)*slopes) # vertical velocities minus horizontal velocities
        b = numpy.transpose(ambient.Vinf*(numpy.cos(ambient.alpha)*slopes - numpy.sin(ambient.alpha)))

        # Lift
        gamma_distribution = numpy.linalg.solve(A, b) # distribution of length-specific vortex strengths
        Vz = ambient.Vinf*numpy.sin(ambient.alpha) + numpy.sum(((panel_lengths*gamma_distribution) / (2*numpy.pi*r)) * numpy.sin(slopes), axis=0)
        Vx = ambient.Vinf*numpy.cos(ambient.alpha) + numpy.sum(((panel_lengths*gamma_distribution) / (2*numpy.pi*r)) * numpy.cos(slopes), axis=0)
        V = numpy.vectorize(numpy.sqrt)(Vx**2 + Vz**2)
        Cp = 1 - (V/ambient.Vinf)**2
        L = float(ambient.rhoinf * ambient.Vinf * sum(gamma_distribution * panel_lengths[0]))
        L_prime = ambient.rhoinf * ambient.Vinf * gamma_distribution * panel_lengths[0]
        Cl = float(L / ((1/2)*ambient.rhoinf*ambient.Vinf**2))

        return L, Cl, L_prime, Cp
    
    @property
    def moment(self):
        if self.flap_start != None:
            theta_f = numpy.acos(1 - (2*self.flap_start))
            A1 = 4*self.M / (numpy.pi * (1 - self.P)**2) * ((1 - self.chord/2)*numpy.sin(theta_f) + (theta_f*self.chord/2) + (numpy.sin(2*theta_f)/4)) - ((2/numpy.pi)*numpy.tan(-self.flap_angle)*numpy.sin(theta_f))
            A2 = 4*self.M / (numpy.pi * (1 - self.P)**2) * ((1 - self.chord/2)*(theta_f/2 + numpy.sin(2*theta_f)/4) + (self.chord/2)*(numpy.sin(theta_f) - (1/3)*(numpy.sin(theta_f)**3))) + (2/numpy.pi*numpy.tan(-self.flap_angle)*((numpy.pi - theta_f)/2 - numpy.sin(2*theta_f)/4))
        else:
            A1 = self.chord * self.M / (1 - self.P)**2
            A2 = - self.chord * self.M / (2 * (1 - self.P)**2)
        cmc4 = (numpy.pi/4) * (A2 - A1)
        return cmc4

    def plot_results(self):
        fig, [ax0, ax1, ax2, ax3] = plt.subplots(2, 2, figsize=(12, 5))

        # Camberline
        x_positions = list()
        y_positions = list()
        for panel in self.panels:
            x_positions.append(float(panel.start_point.x_coord))
            x_positions.append(float(panel.end_point.x_coord))
            y_positions.append(float(panel.start_point.y_coord))
            y_positions.append(float(panel.end_point.y_coord))
        
        control_x = [point.x_coord for point in self.quarter_points.flatten()]
        control_y = [point.y_coord for point in self.quarter_points.flatten()]
        vortex_x = [point.x_coord for point in self.three_quarter_points.flatten()]
        vortex_y = [point.y_coord for point in self.three_quarter_points.flatten()]
        ax0.scatter(control_x, control_y, color="blue", label="Control Points")
        ax0.scatter(vortex_x, vortex_y, color="red", label="Flow Tangency Points")
        ax0.scatter(x_positions, y_positions, color='green', label="Endpoints")
        ax0.plot(x_positions, y_positions)
        ax0.xlim(0, 1)
        ax0.ylim(-0.2, 0.2)
        ax0.legend()

        _, _, L_prime, Cp = self.lift
        x_c = numpy.transpose(numpy.array([point.x_coord for point in self.three_quarter_points[0]]))
        
        # Lift Distribution
        ax1.plot(x_c, L_prime, label="L_prime")
        if self.flap_start != None: ax1.title("Flapped " + self.designation + " Lift Distribution")
        else: ax1.title(self.designation + " Lift Distribution")

        # Pressure Distribution
        ax2.plot(x_c, Cp, label="Cp")
        if self.flap_start != None: ax2.title("Flapped " + self.designation + " Pressure Distribution")
        else: ax2.title(self.designation + " Pressure Distribution")

        plt.show()


# (2D) Source Panel Method uniform flow + octagon circle, 
# J. Anderson, Ch 3.17, pg 284 - 294

import matplotlib.pyplot as plt
import warnings
import math
import numpy
import sys
import os

warnings.filterwarnings("ignore")

# Change the current working directory to the file location
filepath = os.path.abspath(__file__)
directory = os.path.dirname(filepath)
os.chdir(directory)

# Add custom classes to search locations
sys.path.append(r'..\curves')

# Import the Position Vector class
from Point import Point
from BSpline import BSpline

degree = 3
num_points = 50
control_points = [Point(0,0,0), Point(2,2,0), Point(5,4,0), Point(7,4,0)]
spline = BSpline(control_points, degree, num_points)
spline_points = spline.get_positions() # airfoil points

V = 1 # freestream

# 2D body defintion
x_values = [-0.9239, -0.9239, -0.3827, 0.3827, 0.9239, 0.9239, 0.3827, -0.3827]
y_values = [-0.3827, 0.3827, 0.9239, 0.9239, 0.3827, -0.3827, -0.9239, -0.9239]
#x_values = [point.x_coord for point in spline_points]
#y_values = [point.y_coord for point in spline_points]

x_values.append(x_values[0])
y_values.append(y_values[0])
X = numpy.array(x_values)
Y = numpy.array(y_values)

# panel length, mid-point and orientation angle
n = len(X) - 1
Sj = numpy.zeros((n,1))
phi = numpy.zeros((n,1))
xmp = numpy.zeros((n,1))
ymp = numpy.zeros((n,1))
for i in range(1, n):
    # Length of each panel
    Sj[i] = math.sqrt((X[i+1] - X[i])**2 + (Y[i+1] - Y[i])**2)
    # mid-point of each panel
    xmp[i] = 0.5*(X[i+1] + X[i])
    ymp[i] = 0.5*(Y[i+1] + Y[i])
    # bottom of panel orientation angle wrt positive x-axis, CCW positive
    dx = X[i+1] - X[i]
    dy = Y[i+1] - Y[i]
    angle = math.atan(numpy.where(dx != 0, dy / dx, math.inf)) # conditional division
    phi[i] = math.degrees(angle) % 360

plt.plot(X,Y,'o-')
plt.plot(xmp,ymp,'*')
#plt.show()

# build the In(i,j) matrix  
In = numpy.identity(n) * math.pi
Bi = numpy.zeros((n,1))
for i in range(1, n):
    for j in range(1, n):
        if i != j:
            phi_i = math.radians(phi[i])
            phi_j = math.radians(phi[j])
            A = -(xmp[i] - X[j])*math.cos(phi_j) -(ymp[i] - Y[j])*math.sin(phi_j)
            B = (xmp[i] - X[j])**2 + (ymp[i] - Y[j])**2
            C = math.sin(phi_i - phi_j)
            D = (ymp[i] - Y[j])*math.cos(phi_i) - (xmp[i] - X[j])*math.sin(phi_i)
            #E = (xmp[i] - X[j])*math.sin(phi_j) - (ymp[i] - Y[j])*math.cos(phi_j) 
            E = math.sqrt(B - A**2)

            # normal component
            In[i,j] = (C/2)*math.log((Sj[i]**2 + 2*A*Sj[i] + B)/B) + ((D - A*C)/E)*(math.atan((Sj[i] + A)/E) - math.atan(A/E))

    # Build the Bi vector
    Bi[i] = -2 * math.pi * V * math.cos(phi_i + math.pi/2)

# Determine source/sink strength lamda per panel
lower_lambda = numpy.linalg.inv(In)*Bi

# Verify sum of lamda*Lj = check = 0
check = 0
for j in range(1, n):
    check += lower_lambda[j]*Sj[j]

print(sum(check))


'''

# Build the It(i,j) matrix
It = zeros(n,n);  Vi = zeros(n,1); Cp = zeros(n,1)
for i = 1:n
    sum_lamda = 0
    for j = 1:n
        if i ~= j
            A = -(xmp(i) - X(j))*cosd(phi(j)) ...
                -(ymp(i) - Y(j))*sind(phi(j))
            B = (xmp(i) - X(j))^2 + (ymp(i) - Y(j))^2
            C = sind(phi(i) - phi(j))
            D = (ymp(i) - Y(j))*cosd(phi(i)) ...
                - (xmp(i) - X(j))*sind(phi(i))
            E = sqrt(B - A^2)
            # tangential component
            It(i,j) = (D - A*C)/(2*E)*log((Sj(i)^2 + 2*A*Sj(i) + B)/B) ...
                      - C*(atan((Sj(i) + A)/E) - atan(A/E))
        sum_lamda = sum_lamda + lamda(j)/(2*pi)*It(i,j)
    # Determine Vi to find Cp
    Vi(i) = V*sind(phi(i) + 90) + sum_lamda
    Cp(i) = (1 - (Vi(i)/V)^2)


beta = zeros(n,1)
for i = 1:n
    beta(i) = phi(i) + 90
    if beta(i) >= 360:
        beta(i) = beta(i) - 360
[theta,idx] = sort(beta)
theta_p1 = [theta;360]
Cp = Cp(idx)
Cp_p1 = [Cp;Cp(1)]

subplot(2,1,2)
cs = spline(theta_p1,Cp_p1)
xq = linspace(theta_p1(1),theta_p1(n+1),100)
plot(theta_p1,Cp_p1,'r*',xq,ppval(cs,xq),'k-')
label_1 = string([5,4,3,2,1,8,7,6,5])
text(theta_p1,Cp_p1,label_1,'VerticalAlignment','bottom')
'''
""" Module for shape design """
import math
import matplotlib.pyplot as plt
import copy
import numpy

# This is a point, but the truer representation of the object is a position vector. Will manipulate it as such.
class Point:
    def __repr__(self):
            return f"Point: {round(float(self.x_coord), 5)}, {round(float(self.y_coord), 5)}, {round(float(self.z_coord), 5)}"

    def __init__(self, x, y, z):
        self.x_coord = x
        self.y_coord = y
        self.z_coord = z

    def __add__(self, other_point):
        return Point(self.x_coord + other_point.x_coord, self.y_coord + other_point.y_coord, self.z_coord + other_point.z_coord)

    def __sub__(self, other_point):
        return Point(self.x_coord - other_point.x_coord, self.y_coord - other_point.y_coord, self.z_coord - other_point.z_coord)

    def __mul__(self, num):
        return Point(self.x_coord * num, self.y_coord * num, self.z_coord * num)

    def __rmul__(self, num):
        return self.__mul__(num)
    
    def __truediv__(self, scalar):
        return Point(self.x_coord / scalar, self.y_coord / scalar, self.z_coord / scalar)
    
    def make_copy(self):
        return Point(copy.deepcopy(self.x_coord), copy.deepcopy(self.y_coord), copy.deepcopy(self.z_coord))

    def set_x(self, new_x):
        self.x_coord = new_x

    def set_y(self, new_y):
        self.y_coord = new_y

    def set_z(self, new_z):
        self.z_coord = new_z
    
    @property
    def magnitude(self):
        return numpy.sqrt(self.x_coord**2 + self.y_coord**2 + self.z_coord**2)

'''
Bezier Curve class for curve design
'''
class BezierCurve():
    def __init__(self, control_points, ps_ss=None):
        self.control_points = control_points
        self.positions = []
        self.ps_ss = ps_ss
        self.degree = len(self.control_points) - 1


    def get_positions(self, parameters):
        '''
        Returns positions based on the parameters given
        '''
        # Returns the bernstein polynomial coefficient
        def basis_polynomial(parameter, iterator):
            binomial_coefficient = math.factorial(self.degree) / (math.factorial(iterator) * math.factorial(self.degree - iterator))

            return binomial_coefficient * (parameter**iterator) * ((1 - parameter)**(self.degree - iterator))

        # Clear any existing positions
        self.positions.clear()

        # Iterate Through Each Parameter Step
        for t in parameters:
            position = Point(0, 0)

            # Apply the effects of each control point to the parameter
            for idx, point in enumerate(self.control_points):
                position += point.scalar_mul(basis_polynomial(t, self.degree, idx))

            # Store the position
            self.positions.append(position)

        return self.positions


    def plot_points(self):
        '''
        Plots the positions on the curve
        '''
        if len(self.positions) == 0:
            print("No positions have been generated. Parameter(s) must be provided")
            return

        x_positions = [point.x_coord for point in self.positions]
        y_positions = [point.y_coord for point in self.positions]
        plt.scatter(x_positions, y_positions, 'o')
        plt.plot(x_positions, y_positions)

'''
Create a BSpline curve 
Uses the Cox-de Boor recursion formula to determine the basis coefficients

NOTE: This program currently assumes a uniform knot vector. Different types of knot vectors must be accounted for
'''
class BSpline():
    def __init__(self, control_points, degree, num_points):
        self.control_points = control_points
        self.degree = degree
        self.positions = []

        # Uniform knot vector 
        num_knots = degree + len(self.control_points) + 1
        self.knots = numpy.linspace(0, num_knots, num=num_knots)

        # Parameter values
        min_parameter = self.knots[0]
        max_parameter = self.knots[-1]
        self.parameters = numpy.linspace(min_parameter, max_parameter, num_points)

    def get_positions(self):
        '''
        Returns the positions of each point based on the control points
        '''
        def basis_function(parameter, degree, i):
            # Check for curve simplicity
            if degree == 0:
                if self.knots[i] <= parameter < self.knots[i + 1]:
                    return 1
                return 0

            # First Term
            if self.knots[i + degree] == self.knots[i]:
                term1 = 0
            else:
                term1 = ((parameter - self.knots[i]) / (self.knots[i + degree] - self.knots[i])) * basis_function(parameter, degree - 1, i)
            
            # Second Term
            if self.knots[i + degree + 1] == self.knots[i + 1]:
                term2 = 0
            else:
                term2 = ((self.knots[i + degree + 1] - parameter) / (self.knots[i + degree + 1] - self.knots[i + 1])) * basis_function(parameter, degree - 1, i + 1)

            return term1 + term2
        
        for parameter in self.parameters:
            point = Point(0,0,0)
            for i in range(len(self.control_points)):
                basis_coefficient = basis_function(parameter, self.degree, i)
                point += self.control_points[i] * basis_coefficient
            self.positions.append(point)

        return self.positions

    def plot_points(self):
        '''
        Plots the spline points.
        '''
        if len(self.positions) == 0:
            print("No positions have been generated. Parameter(s) must be provided")
            return

        # Retrieve the x-components and y-components of the spline points and control points
        spline_x = [point.x_coord for point in self.positions]
        spline_y = [point.y_coord for point in self.positions]
        control_x = [point.x_coord for point in self.control_points]
        control_y = [point.y_coord for point in self.control_points]

        # Plot all points
        plt.scatter(spline_x, spline_y, label = "Resultant BSpline Curve")
        plt.plot(spline_x, spline_y, label = "Resultant BSpline Curve")
        plt.scatter(control_x, control_y)
        plt.show()

'''
""" BSpline code from ChatGPT """
def find_knot_span(u, degree, knots, n):
    """Find the knot span index for parameter u."""
    if u >= knots[n + 1]:  # Clamp to last span
        return n
    if u <= knots[degree]:  # Clamp to first span
        return degree
    for i in range(degree, n + 1):
        if knots[i] <= u < knots[i + 1]:
            return i
    raise ValueError("Parameter u is out of bounds.")

def de_boor(u, degree, knots, control_points):
    """Evaluate B-spline curve at parameter u using De Boor's algorithm."""
    n = len(control_points) - 1
    k = find_knot_span(u, degree, knots, n)

    d = [np.array(control_points[j]) for j in range(k - degree, k + 1)]

    for r in range(1, degree + 1):
        for j in range(degree, r - 1, -1):
            i = k - degree + j
            denom = knots[i + degree - r + 1] - knots[i]
            alpha = (u - knots[i]) / denom if denom != 0 else 0
            d[j] = (1.0 - alpha) * d[j - 1] + alpha * d[j]

    return d[degree]


degree = 3
control_points = [(0, 0), (1, 2), (3, 5), (6, 5), (7, 2), (8, 0)]
n = len(control_points) - 1

# Normalized, clamped uniform knot vector
knots = [0] * (degree + 1) + list(np.linspace(0, 1, n - degree + 1)) + [1] * (degree + 1)

# Evaluate B-spline at 100 points
curve = [de_boor(u, degree, knots, control_points) for u in np.linspace(0, 1, 100)]

# Optional plot
if __name__ == "__main__":
    import matplotlib.pyplot as plt

    x_vals, y_vals = zip(*curve)
    ctrl_x, ctrl_y = zip(*control_points)

    plt.plot(x_vals, y_vals, label="B-spline Curve")
    plt.plot(ctrl_x, ctrl_y, 'o--', label="Control Polygon")
    plt.legend()
    plt.axis('equal')
    plt.title("B-spline Curve (Robust De Boor)")
    plt.show()
'''

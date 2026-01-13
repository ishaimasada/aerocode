from scipy.interpolate import CubicSpline
import numpy, csv, sys, math, os

# Change the current working directory to the file location
filepath = os.path.abspath(__file__)
directory = os.path.dirname(filepath)
os.chdir(directory)

# Add custom classes to search locations
sys.path.append(r'..\curves')

# Import the Position Vector class
from Point import Point

def write_coords(filename, coordinates):
    '''
    Writes the injection locations to a CSV file
    '''
    # Erase the contents of the file
    open(filename, 'w').close()

    # Open the file to use csv module
    with open(filename, 'w', newline='\n') as csvfile:
        coord_writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')

        # File Format Requirements
        csvfile.write('[Name]\nR1 Blade\n\n')
        csvfile.write('[Spatial Fields]\n')
        coord_writer.writerow(['X', 'Y', 'Z'])
        csvfile.write('\n[Data]\n')
        coord_writer.writerow(['X[m]', 'Y[m]', 'Z[m]'])
        csvfile.write('\n')

        theta = 5 * math.pi/180 # assumed angle in radians

        # Write each set of coordinates as a row in the file
        for point in coordinates:
            # Offset
            offset = point[0] * math.tan(theta) # add to y-coordinates

            coord_writer.writerow([point[0], point[1] + offset, point[2]])

def load_data(span, data):
    # All points within airfoil section
    #section_points = [(float(point.strip().split()[0]), float(point.strip().split()[1]), float(point.strip().split()[2])) for point in raw_data[span].splitlines()]
    section_points = [(float(point.strip().split()[1]), float(point.strip().split()[2])) for point in data[span].splitlines()]

    # Span Position
    span_x = int(float(data[span].splitlines()[0].split()[0]))

    # Suction Side Points
    #suction_side_points = [(float(point.strip().split()[0]), float(point.strip().split()[1]), float(point.strip().split()[2])) for point in raw_data[span].split('TE')[0].splitlines()]
    pressure_side_points = [(float(point.strip().split()[1]), float(point.strip().split()[2])) for point in data[span].split('TE')[0].splitlines()]

    # Pressure Side Points
    suction_side_points = data[span].split('TE')[1].splitlines()
    suction_side_points.pop(0)
    suction_side_points = [(float(point.strip().split()[1]), float(point.strip().split()[2])) for point in suction_side_points]

    return span_x, suction_side_points, pressure_side_points

def get_position(chord_position, span_position, pressure_suction_choice):
    match pressure_suction_choice:
        case 'ps':
            parameter_0 = numpy.linspace(0, 1, len(pressure_0))
            span0_section = CubicSpline(parameter_0, pressure_0)
            initial_point = span0_section(chord_position)

            parameter_1 = numpy.linspace(0, 1, len(pressure_1))
            span1_section = CubicSpline(parameter_1, pressure_1)
            final_point = span1_section(chord_position)

        case 'ss':
            parameter_0 = numpy.linspace(0, 1, len(suction_0))
            span0_section = CubicSpline(parameter_0, suction_0)
            initial_point = span0_section(chord_position)

            parameter_1 = numpy.linspace(0, 1, len(suction_1))
            span1_section = CubicSpline(parameter_1, suction_1)
            final_point = span1_section(chord_position)

    desired_x = (span_0 + span_position*(span_1 - span_0))
    initial_point = Point(desired_x, initial_point[0], initial_point[1]).scalar_mul(0.001)
    final_point = Point(desired_x, final_point[0], final_point[1]).scalar_mul(0.001)
    section_point = initial_point + (final_point - initial_point).scalar_mul(span_position)

    return section_point


# Read in user input
with open(r"data/input positions.csv", 'r') as f:
    user_input = f.read().splitlines()
    user_input.pop(0)

points = []

for row in user_input:
    args = row.split(',')

    span_position = float(args[0])
    chord_position = float(args[1])
    ps_ss = args[2]
    rotor_stator = args[3]

    match rotor_stator:
        case "rotor":
            # Read in rotor data
            with open('data/rotor.txt', 'r') as f:
                raw_data = f.read().split('\n\n')

                span_0, suction_0, pressure_0 = load_data(0, raw_data)
                span_1, suction_1, pressure_1 = load_data(4, raw_data)

        case "stator":
            # Read in stator data
            with open('data/stator.txt', 'r') as f:
                raw_data = f.read().split('\n\n')

                span_0, suction_0, pressure_0 = load_data(0, raw_data)
                span_1, suction_1, pressure_1 = load_data(4, raw_data)

    injection_location = get_position(chord_position, span_position, ps_ss)
    points.append([injection_location.x_coord, injection_location.y_coord, injection_location.z_coord])

write_coords('Injection Locations.csv', points) 

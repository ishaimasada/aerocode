''' Example usage code for a Radial Compressor '''
import json
import os
import sys

# Change the current working directory to the file location
filepath = os.path.abspath(__file__)
directory = os.path.dirname(filepath)
os.chdir(directory)

# Load parameters from JSON file
with open("compressor_parameters.json", "r") as file:
    parameters = json.load(file)

sys.path.append(r"../propulsion")

# Import all types from engine module
from engine import Compressor

compressor = Compressor(upstream=None, component_parameters=parameters)

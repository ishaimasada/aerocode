''' Turbojet engine design for XJEP --> Replaces the role of the JetCats and KingTechs by allowing for easy instrumentation '''

import json
import sys
import os

# Change the current working directory to the file location
filepath = os.path.abspath(__file__)
directory = os.path.dirname(filepath)
os.chdir(directory)

# Load engine and component parameters from JSON files
with open("turbojet_parameters.json", "r") as file:
    cycle_parameters = json.load(file)["parameters"]
    engine_parameters = cycle_parameters["engine"]

with open("inlet_parameters.json", "r") as file:
    inlet_parameters = json.load(file)

with open("compressor_parameters.json", "r") as file:
    compressor_parameters = json.load(file)

with open("burner_parameters.json", "r") as file:
    burner_parameters = json.load(file)

with open("turbine_parameters.json", "r") as file:
    turbine_parameters = json.load(file)

with open("nozzle_parameters.json", "r") as file:
    nozzle_parameters = json.load(file)

sys.path.append(r"..\propulsion")

# Import all types from engine module
from engine import * # type: ignore

# CYCLE ANALYSIS
output_filename = "station_data.xlsx"

# Create an instance of the Engine class
engine = Engine(engine_parameters) # type: ignore

# Pass the component parameters to the engine object for cycle analysis
engine.set_components(cycle_parameters)

# Change the current working directory to the file location
os.chdir(directory)

# COMPONENT DESIGN
engine.inlet.design_component(inlet_parameters)
engine.compressor.design_component(compressor_parameters)
engine.burner.design_component(burner_parameters)
engine.turbine.design_component(turbine_parameters)
engine.exhaust.design_component(nozzle_parameters)
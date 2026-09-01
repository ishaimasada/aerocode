''' Turbojet engine design for XJEP --> Replaces the role of the JetCats and KingTechs by allowing for easy instrumentation '''

# Import modules
import json, sys, os

sys.path.append(r"..\propulsion")
from engine import * # type: ignore

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

# CYCLE ANALYSIS
cycle_filename = "cycle.xlsx"
engine = Engine(engine_parameters) # type: ignore
engine.set_components(cycle_parameters)
os.chdir(directory)
engine.write_station_data(cycle_filename) 

# COMPONENT DESIGN
engine.inlet.design_component(inlet_parameters)
engine.compressor.design_component(compressor_parameters)
engine.burner.design_component(burner_parameters)
engine.turbine.design_component(turbine_parameters)
engine.exhaust.design_component(nozzle_parameters)

# Report Results
inlet_coords = engine.inlet.get_results()
compressor_r_coords, compressor_z_coords = engine.compressor.get_results()
engine.turbine.get_results(turbine_parameters["flags"])
burner_thermo, burner_geometry = engine.burner.get_results()
exhaust_coords = engine.exhaust.get_results()
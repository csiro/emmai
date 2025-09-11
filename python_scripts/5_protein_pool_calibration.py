#!/usr/bin/env python
import cobra
import os
import yaml

inputs_path = os.getenv("INPUTS")  # From your env.sh file
if inputs_path is None:
    raise ValueError("The INPUTS environment variable is not set.")
inputs_file = os.path.join(inputs_path, "inputs.yml")
with open(inputs_file, "r") as file:
    data = yaml.safe_load(file)

# Input and output file paths
output_file_path = os.path.join(inputs_path, data["output_file_path"])
sbml_model = os.path.join(inputs_path, data["sbml_model"])
modified_model_name = os.path.splitext(data["sbml_model"])[0]
modified_model_file = os.path.join(
    output_file_path, "output_GEMs", f"patched_ec_{modified_model_name}_final.xml"
)
os.makedirs(output_file_path, exist_ok=True)
transporters = data["transporters"]
media = data["media"]
bounds = tuple(data["bounds"])

ec_model = cobra.io.read_sbml_model(
    os.path.join(
        output_file_path, "output_GEMs", f"patched_ec_{modified_model_name}_fixed.xml"
    )
)

ec_model.medium = media
ec_model.reactions.DM_usage.bounds = bounds
sol = ec_model.optimize()

print(ec_model.summary(sol))
# print(sol.objective_value)

cobra.io.write_sbml_model(
    ec_model, os.path.join(output_file_path, "output_GEMs", modified_model_file)
)

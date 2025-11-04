#!/usr/bin/env python
import cobra
import os
import yaml
import logging

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
    output_file_path, "output_GEMs", f"ec_{modified_model_name}_mod2.xml"
)
os.makedirs(output_file_path, exist_ok=True)
transporters = data["transporters"]
excluded_reactions = data["excluded_reactions"]

logging.getLogger("cobra").setLevel(logging.ERROR)
ec_model = cobra.io.read_sbml_model(
    os.path.join(output_file_path, "output_GEMs", f"ec_{modified_model_name}_mod1.xml")
)

"""Average reaction coefficient"""

usage_coefficients = []

for reaction in ec_model.reactions:
    if "resource_usage_pseudometabolite" in [m.name for m in reaction.metabolites]:
        for m in reaction.metabolites:
            if m.id == "usage":
                coef = reaction.metabolites[m]
                usage_coefficients.append(coef)

average_coef = sum(usage_coefficients) / len(usage_coefficients)
print(abs(average_coef))

patched_model = ec_model.copy()


def contains_keywords(cell):
    return any(keyword.lower() in str(cell).lower() for keyword in transporters)


usage = patched_model.metabolites.get_by_id("usage")

with patched_model:
    for reaction in patched_model.reactions:
        if (
            reaction not in patched_model.boundary
            and not contains_keywords(reaction.name)
            and reaction.name not in excluded_reactions
        ):
            if "resource_usage_pseudometabolite" not in [
                m.name for m in reaction.metabolites
            ]:
                # print(reaction.bounds)
                if reaction.lower_bound < 0 and reaction.upper_bound <= 0:
                    # print('reverse')
                    reaction.add_metabolites({usage: abs(average_coef)})
                if reaction.lower_bound >= 0 and reaction.upper_bound > 0:
                    # print('forward')
                    reaction.add_metabolites({usage: -abs(average_coef)})
                # print(reaction.name, reaction.reaction)

    cobra.io.write_sbml_model(
        patched_model,
        os.path.join(output_file_path, "output_GEMs", modified_model_file),
    )

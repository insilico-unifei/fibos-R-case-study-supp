###############################################################################
# FIBOS: R and Python packages for analyzing protein packing and structure
# 
# PYTHON SCRIPT ACCESSORY FOR: a case study using FIBOS, comparing packing densities 
# between experimental protein structures and AlphaFold predictions.
#
# AUTHORS: 
#
# Herson Soares, João Romanelli, Carlos Silveira
# Federal University of Itajubá (Unifei), Itabira, MG, Brazil.
#
# Patrick Fleming,
# Johns Hopkins University (JHU), Baltimore, MD, USA.
#
# CONTACTS: hersinsoares@gmail.com, carlos.silveira@unifei.edu.br

# '''
# QUERIES RCSB FOR PROTEIN STRUCTURES WITH SPECIFIED LENGTH, COUNT, AND RESOLUTION CRITERIA.
# GROUPS BY SEQUENCE IDENTITY (30% CUTOFF) AND SAVES REPRESENTATIVE ENTITY IDS TO CSV.
# 
# '''

import os
import csv
import sys
from pathlib import Path
from rcsbsearchapi import rcsb_attributes as attrs
from rcsbsearchapi.search import GroupBy

out_path = sys.argv[1]

q1 = attrs.rcsb_entry_info.selected_polymer_entity_types == "Protein (only)"
q2 = attrs.entity_poly.rcsb_sample_sequence_length < 1200
q3 = attrs.entity_poly.rcsb_sample_sequence_length > 50
q4 = attrs.rcsb_entry_info.deposited_polymer_entity_instance_count == 1
q5 = attrs.rcsb_entry_info.assembly_count == 1
q6 = attrs.exptl.method == "X-RAY DIFFRACTION"
q7 = attrs.rcsb_entry_info.resolution_combined <= 1.5

query = q1 & q2 & q3 & q4 & q5 & q6 & q7 


results = list(
    query(
        group_by=GroupBy(aggregation_method="sequence_identity", similarity_cutoff=30),
        group_by_return_type="representatives",
        return_type="polymer_entity"
    )
)

if results:

    total_hits = len(results)
    print(f"Found {total_hits} structures")
    print(results[:10])

    output = results

    print("\nSaving to file: ")
    print(out_path)
    with open(out_path, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows([[x] for x in results])



#!/usr/bin/env python3

import csv
import sys

# Input files
tsv_file = sys.argv[3]
csv_file = sys.argv[2]
output_file = f"{sys.argv[1]}.ninjaMap.abundance.csv"

# Read A.tsv (tab-separated) into a dict keyed by first column
a_dict = {}
with open(tsv_file, newline='') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        key = row[0]
        # store columns 2–4 (index 1–3)
        a_dict[key] = row[1:4]

# Read B.csv (comma-separated)
merged_data = []
with open(csv_file, newline='') as f:
    reader = csv.reader(f)
    for row in reader:
        key = row[0]
        if key in a_dict:
            merged_data.append(row + a_dict[key])
        else:
            # If no match, still include B row (optional)
            merged_data.append(row + ["", "", ""])

# Write output as CSV
with open(output_file, "w", newline='') as f:
    writer = csv.writer(f)
    writer.writerows(merged_data)


'''
# Read A.tsv (tab-separated)
a_data = []
with open(tsv_file, newline='') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        # Take only 2nd, 3rd and 4th columns (index 1 to 3)
        a_data.append(row[1:4])

# Read B.csv (comma-separated)
b_data = []
with open(csv_file, newline='') as f:
    reader = csv.reader(f)
    for row in reader:
        b_data.append(row)

# Merge line by line
merged_data = []
for b_row, a_row in zip(b_data, a_data):
    merged_data.append(b_row + a_row)

# Write output as CSV
with open(output_file, "w", newline='') as f:
    writer = csv.writer(f)
    writer.writerows(merged_data)

print(f"Merged file written to {output_file}")
'''
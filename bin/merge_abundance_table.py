#!/usr/bin/env python3

import csv
import sys

# Input files
tsv_file = sys.argv[3]
csv_file = sys.argv[2]
output_file = f"{sys.argv[1]}.ninjaMap.abundance.csv"

# Read A.tsv (tab-separated)
a_data = []
with open(tsv_file, newline='') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        # Take only 2nd and 3rd columns (index 1 and 2)
        a_data.append(row[1:3])

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

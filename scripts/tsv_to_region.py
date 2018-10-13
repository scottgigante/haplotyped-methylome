import sys
import csv
tsv_fn = sys.argv[1]
output = []
with open(tsv_fn, 'r') as tsv_handle:
    reader = csv.DictReader(tsv_handle, delimiter="\t")
    for row in reader:
        output.append("{}:{}-{}".format(row["chr"], row["start"], row["end"]))

print(" ".join(output))

import sys
import csv
tsv_fn = sys.argv[1]
if len(sys.argv) > 2:
    overhang = int(sys.argv[2])
else:
    overhang = 0
output = []
with open(tsv_fn, 'r') as tsv_handle:
    reader = csv.DictReader(tsv_handle, delimiter="\t")
    for row in reader:
        output.append(
            "{}:{}-{}".format(row["chr"], int(row["start"]) - overhang, int(row["end"]) + overhang))

print(" ".join(output))

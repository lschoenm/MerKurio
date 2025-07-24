# Python Analysis of the JSON Log

[Python](https://www.python.org/) is a general purpose programming language with an easy-to-learn syntax. 

It natively supports parsing JSON files with the `json` module. Replace `log.json` with the path of your log file. It will store the JSON log as a standard Python dictionary:  

```python
import json

# Load the JSON log
with open('log.json') as f:
    data = json.load(f)

# Accessing the "summary_statistics" entry
data["summary_statistics"]
```

## Example: Sorting Records by Match Count

In this example, the JSON log is parsed and the match counts per record are stored in a dictionary. Then, the entries are printed in descending order: 

```python
import json
from collections import defaultdict

# Load the JSON log
with open("log.json") as f:
    data = json.load(f)

# Count matches per distinct record
record_counts = defaultdict(int)
for entry in data["matching_records"]:
    record_id = entry["record_id"]
    record_counts[record_id] += 1

# Sort by match count, descending
sorted_records = sorted(record_counts.items(), key=lambda x: x[1], reverse=True)

# Print results
for record_id, count in sorted_records:
    print(f"{record_id}: {count} match{'es' if count > 1 else ''}")
```

## Example: Plot a Histogram of Match Densities

Using the external dependency Matplotlib, a histogram of match position density can be plotted: 

```python
import json
import matplotlib.pyplot as plt

# Load the JSON log
with open('log.json') as f:
    data = json.load(f)

# Extract and convert positions to integers
positions = [int(entry["position"]) for entry in data["matching_records"]]

# Plot histogram
plt.hist(positions)
plt.xlabel("Pattern match position in record")
plt.ylabel("Frequency")
plt.title("Histogram of Pattern Match Positions")
plt.tight_layout()
plt.show()
```
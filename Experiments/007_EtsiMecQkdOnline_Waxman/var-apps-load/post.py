#!/usr/bin/env python3

import pandas as pd
import sys
from os import getenv

data = pd.read_csv(sys.stdin)

print(f"algorithm    : {data['algorithm'].unique()}")
print(f"arrival_rates: {data['arrival-rate'].unique()}")

data["capacity-res"] = data["capacity-res"] / data["capacity-tot"]
data["processing-res"] = data["processing-res"] / data["processing-tot"]

metrics = [
    "capacity-res",
    "processing-res",
    "blocking-probability",
    "avg-active-apps",
    "path-length-avg",
    "net-rate-tot",
    "signalling-rate",
    "execution-time",
]

for metric in metrics:
    grouped_df = data.groupby(["algorithm", "arrival-rate"])

    # summary statistics
    grouped = grouped_df[metric]
    outdata = dict()
    for (algorithm, arrival_rate), value in grouped.mean().items():
        outdata.setdefault(algorithm, dict())
        outdata[algorithm][arrival_rate] = [value]
    for (algorithm, arrival_rate), value in grouped.std().items():
        outdata[algorithm][arrival_rate].append(value)

    for algorithm, xy_values in outdata.items():
        with open(f"post/{metric}-{algorithm}.dat", "w") as outfile:
            for x_value, y_values in xy_values.items():
                outfile.write(f"{x_value}")
                for y_value in y_values:
                    outfile.write(f" {y_value}")
                outfile.write("\n")

    if getenv("NORAW"):
        continue

    # raw data (e.g., for boxplots or distributions)
    for (algorithm, arrival_rate), values in grouped_df:
        values[metric].to_csv(
            f"post/{metric}-{algorithm}-{arrival_rate}-raw.dat",
            header=False,
            index=False,
        )

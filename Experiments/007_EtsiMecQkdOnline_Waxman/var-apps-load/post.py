#!/usr/bin/env python3

import pandas as pd
import sys
from os import getenv

data = pd.read_csv(sys.stdin)

groups = getenv("KEYS").split(",")
assert len(groups) == 2

for group in groups:
    print(f"{group} values: {data[group].unique()}")

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
    grouped_df = data.groupby(groups)

    # summary statistics
    grouped = grouped_df[metric]
    outdata = dict()
    for (group1, group2), value in grouped.mean().items():
        outdata.setdefault(group1, dict())
        outdata[group1][group2] = [value]
    for (group1, group2), value in grouped.std().items():
        outdata[group1][group2].append(value)

    for group1, xy_values in outdata.items():
        with open(f"post/{metric}-{group1}.dat", "w") as outfile:
            for x_value, y_values in xy_values.items():
                outfile.write(f"{x_value}")
                for y_value in y_values:
                    outfile.write(f" {y_value}")
                outfile.write("\n")

    if getenv("NORAW"):
        continue

    # raw data (e.g., for boxplots or distributions)
    for (group1, group2), values in grouped_df:
        values[metric].to_csv(
            f"post/{metric}-{group1}-{group2}-raw.dat",
            header=False,
            index=False,
        )

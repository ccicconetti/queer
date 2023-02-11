# Experiment 006

Online assignment of paths in a QKD network and tasks to edge nodes.

Topologies generated with the Waxman model:

| Parameter | Value  |
|-----------|--------|
| Alpha     | 0.4    |
| Beta      | 0.4    |
| L         | 100 km |

Allocation algorithms compared:

| Algorithm   | Description                                                                 |
|-------------|-----------------------------------------------------------------------------|
| Random      | Pick an app at random among those feasible                                  |
| Spf         | Select the app with shortest path among those feasible                      |
| Bf          | Select the app with minimum residual on the edge nodes among those feasible |
| RandomBlind | Pick an app at random                                                       |
| SpfBlind    | Select the app with shortest path                                           |
| BfBlind     | Select the app with minimum residual on the edge nodes                      |

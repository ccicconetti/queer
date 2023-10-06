# Experiment 007

Online assignment of paths in a QKD network and tasks to edge nodes, with dynamic applications that enter/leave the system.

Topologies generated with the Waxman model:

| Parameter | Value  |
| --------- | ------ |
| Alpha     | 0.4    |
| Beta      | 0.4    |
| L         | 100 km |

Allocation algorithms compared:

| Algorithm       | Description                                  |
| --------------- | -------------------------------------------- |
| Policy014_k1    | all applications follow the same path        |
| Policy014_k3    | same as above, but keep k=3 alternatives     |
| Policy015       | each application is assigned own path        |
| Policy015_reuse | same as above, but try to re-use allocations |

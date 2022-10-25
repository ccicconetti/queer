# Experiment 005

Peer assignment and resource scheduling of apps in a quantum network.

Topologies generated with a Poisson Point Process in a circle with radius 60 km:

| Parameter      | Dense | Sparse |
|----------------|-------|--------|
| Avg num users  | 50    | 100    |
| Threshold (km) | 20    | 15     |

Scheduling algorithm: DRR

Peer assignment algorithms: random, shortest-path, load-balancing

## Sub-experiment var-num-peers

- Fraction of end users: 20%
- Fraction of data centers: 20%
- Number of peers: from 1 to 5
- Number of apps: 10 or 20

## Sub-experiment var-users

- Fraction of end users: from 10% to 50%
- Fraction of data centers: 20%
- Number of peers: 3
- Number of apps: 15

## Sub-experiment var-dcs

- Fraction of end users: 20%
- Fraction of data centers: from 10% to 50%
- Number of peers: 3
- Number of apps: 15

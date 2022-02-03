# Infection

This package will allow you to create a credible set for the patient zero given a set of infected nodes on a graph.

## Installation

Run the following to install:

```python
pip install infect-net-inference
```

## Usage

```python
from infect_tools_noisy_conv import *

# Simulate an infection on graf. Outputs a list containing the infected nodes in the order they were infected 
order = simulateInfection(graf, first, n_inf, q) 

# Create a dict mapping each node index to the posterior root probability
freq = inferInfection(foo, q, **mcmc_params)
    

```

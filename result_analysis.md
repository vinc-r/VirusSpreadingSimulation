# Virus Spreading Simulation

Here is a visualization of a Virus Spreading Simulation (VSS).

Each point can represent a person. A first point is infect by the virus.
Different parameters affect the propagation. Points can propagate virus between themselves.

Here are some points to understand before seeing :
- The aim is to see differences with or without an isolation of the population.
- The quarantine area represent an area where sick point can't infect any healthy point.
- When the containment threshold is achieve : majority of population don't move


## With containment and quarantine simulation

![](https://github.com/vinc-r/VirusSpreadingSimulation/blob/master/quarantine_and_containment_result.gif?raw=true)

The containment really help to educe the propagation. We also have recovered people propagate very well the immunity.

## With quarantine but without containment simulation

![](https://github.com/vinc-r/VirusSpreadingSimulation/blob/master/only_quarantine_result.gif?raw=true)

Without containment, infected points still continue to propagate too much. Even if there is a quarantine area for sick points.

## No measure simulation

![](https://github.com/vinc-r/VirusSpreadingSimulation/blob/master/no_measure_result.gif?raw=true)

Without any measure to contain, the virus propagate very wel between points. 97% of the population is touch by the virus at the end (only 3% get a chance to be immune before).
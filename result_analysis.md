# Virus Spreading Simulation

-----
author: Vincent Rosset
date: April 03rd 2020
-----

How to visualize the global impact of containment ? To understand, I choose to code in python a **Virus Spreading Simulation (VSS)** considering multiple factors close to the real situation: 
- infection radius and probability, death rate, incubation period, healing duration, asymptomatic rate, containment threshold, % of population staying contain, immunity spreading.



Each dot represents a person. A first point is infect by the virus.
Different parameters affect the propagation. Dots can propagate the virus among themselves.

Here are some points to understand before looking at the simulations:
- The aim is to see the effects of isolation on people.
- The quarantine area represents an area where sick dots cannot infect any healthy dot.
- When the containment threshold is achieved a majority of dots does not move.


## With containment and quarantine simulation

![](https://github.com/vinc-r/VirusSpreadingSimulation/blob/master/quarantine_and_containment_result.gif?raw=true)

The containment really helps to reduce the propagation. Recovered people also propagate the immunity very well.

## With quarantine but without containment simulation

![](https://github.com/vinc-r/VirusSpreadingSimulation/blob/master/only_quarantine_result.gif?raw=true)

Without containment, infected dots still continue to propagate a lot -even if there is a quarantine area for sick dots.

## No measure simulation

![](https://github.com/vinc-r/VirusSpreadingSimulation/blob/master/no_measure_result.gif?raw=true)

Without any measure to contain, the virus propagates very well among dots. 
97% of the population is affected by the virus at the end and only 3% had a chance to be immune before. 


__See code : [https://github.com/vinc-r/VirusSpreadingSimulation/blob/master/VirusSpredingSimulation.py](https://github.com/vinc-r/VirusSpreadingSimulation/blob/master/VirusSpredingSimulation.py)__
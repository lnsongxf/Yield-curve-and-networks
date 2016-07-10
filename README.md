# Yield-curve-and-networks
with Can Gao and Andrea Buraschi

This project is the simulation basis for the paper we are writing.

It contains code to produce graphs such as yield curves, risk premia and others.

Main library contains:
general files, network type + _n + load
load signifies that it loads the states of the world from data library

save: contains files for each topology that save the states of the world
The save file allows us to not compute every time again all the states with Pi matrix
so effectively the Pi matrix is very easy to compute because we accounted for states that are
symmetric.

the next stage would be to simply save the Pi matrices, though that would take a lot of space.

Topologies:
There are several topologies: iid (no network), ring, star, chain.

The main files produce the following graphs:
conditional covariance (rainbow heat-map)

bond premia and yield curves

yield curves -- all in one figure, with x distressed in each subfigure

choice of one sub-figure of the above

In addition there are, for each type of network, a separate file that calculates the spectral gap
and produces an appropriate graph

There are a couple of files that produce pictures of the states of the world:
k_reg_pic for k-regular networks
star_pic for star
chain_pic for chain


Finally, there is file that produces the table of amount of states in relation to the topology
and the number of nodes.


no_load is a primitive version of the code that calculates the Pi matrix without accounting
for symmetric states. It works perfectly fine but is slower.

data libaray contains what save files produce and what load file require to work

figures are saved into figures dir


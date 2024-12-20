# greedy-peering
Code for simulations done as part of the paper "Stability of P2P Networks Under Greedy Peering", SIROCCO'24 [best paper award]

All data is stored under the 'data' directory (~3GBs).

The simulations begin with an initial graph topology (generally a random graph) with $n$ nodes, each making $d$ outgoing connections and 
accepting $d_{in}$ incoming connections ($d_{in} =n$ same as uncapped). In each round each node drops $k$ (generally set to 1) of its 
worst-performing peers and replaces them with $k$ new random peers. Performance is based on avarage distance to some $m$ set of miners. 

The are two models of the simulation. The primary model used in the paper is that all nodes drop their worst edges first, then in a random 
order they each replace the lost edges, and at the end all the distances are updated. For a limited set of simulations (called the 'ideal') 
simulations, in a round nodes take turns dropping their worst edge, replacinc it, and all distances being updated. Unless it is labeled 
'ideal', all simulations follow the former model.

For more information on the simulation setup, check out the paper.

Code to generate the simulation data is found under 'simulation code', and all analysis and plotting code can be found in the included jupyter notebooks.






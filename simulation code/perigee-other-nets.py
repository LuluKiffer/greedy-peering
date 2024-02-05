# run simulations of an initail random graph with nodes following the perigee greedy protocol to get closer to the miners

import sys
import random
import time
#import matplotlib.pyplot as plt
import numpy as np
import os
from library import Graph
from datetime import datetime
import json

global n,m,g,miners,vertices,v_tag,max_in,D
class Node:
    def __init__(self, i, miner=0, frac = 0):
        self.num  = i
        self.in_edges = []
        self.out_edges = []
        self.is_miner = miner
        if miner:
            miners.append(self)
            self.frac = frac
        else:
            self.frac = 0
        self.miner_dist = 1000
        self.last_out = 0
    def update_miners_dist(self,D):
        # gets the average distance to all miners (if a miner, the distance to self is 0)
        self.miner_dist = 0.0
        for m in miners:
            self.miner_dist += m.frac*D[m.num][self.num]
    def remove_worst_out(self):
        # go through out_edges and removes the worst peer (i.e. the one with maximum miner_dist)
        d = 0
        m = 0
        for v in self.out_edges:
            if v.miner_dist > d:
                d = v.miner_dist
                m = v
        self.out_edges.remove(m)
        m.in_edges.remove(self)
        g.remove_edge(self.num, m.num)
        self.last_out = m
        return
    def add_new_edge(self,num):
        # adds a random new out edges from list of all non-peers and not the last one removed
        poss = set(vertices)-set(self.in_edges)-set(self.out_edges)-set([self])
        if self.last_out != 0:
            poss-= set([self.last_out])
        count = 0
        tries = 0
        while(count<num):
            tries += 1
            if(len(poss)==0):
                break
            new = random.choice(list(poss))
            if(len(new.in_edges) < max_in):
                g.add_edge(self.num, new.num,1)
                self.out_edges.append(new)
                new.in_edges.append(self)
                poss.remove(new)
                count += 1
            else:
                poss.remove(new)
                continue
        return tries

def update_miners_dist(m):
    D = g.dijkstra_all().copy()
    diam_n = 0
    diam_m = 0
    max_d = [0]*n # calculate the eccentricity of each vertex

    # get diameter of network and miners
    for i in range(n):
        for j in range(i,n):
            if(D[i][j]>diam_n):
                diam_n = D[i][j]
            if(j<m and D[i][j]>diam_m):
                diam_m = D[i][j]
            if(D[i][j] > max_d[i]):
                max_d[i] = D[i][j]
            if(D[i][j] > max_d[j]):
                max_d[j] = D[i][j]
    return (D, diam_n, diam_m, max_d)
    # return distance matrix, the diameter of the network, the diameter of the miners
    # and a vector of eccentricities

def get_pdf_cdf(list_input):
    counts = {}
    for d in list_input:
        d = round(d, 2)
        try:
            counts[d] += 1
        except:
            counts[d] = 1
    xs = sorted(counts.keys())
    pdf = []
    cdf = []
    total = 0
    for x in xs:
        try:
            v = counts[x]
        except:
            v = 0
        total += v
        pdf.append(v)
        cdf.append(total)
    return((xs,pdf),(xs,cdf))

# creates a random graph with n miners and runs perigee for r rounds
# returns the in_degree for each node per round and the cdf and pdf of the in degrees per round
def run_simulation(num_nodes,num_miners,d_out,r,d_in, graph, beta):
    # number of nodes, miners, out degree, in degree, rounds, and distribution of mining power
    global n,m,g,miners,vertices,v_tag,max_in

    n = num_nodes
    m = num_miners
    d = d_out
    max_in = d_in

    miners = [] # list of miners
    vertices = [] # list of vertices
    g = Graph(n)

    dist = [1.0/m]*m

    # make all the vertices
    for i in range(0,n):
        if(i<m):
            vertices.append(Node(i,1,dist[i]))
        else:
            vertices.append(Node(i))

    # make a small world graph
    if(graph == 'small'):
        # for each node connect it to the d_out nodes after it (mod n)
        for i in range(0,n):
            for j in range(0,d_out):
                m = (i+j)%n
                g.add_edge(i,m,1)
                vertices[i].out_edges.append(vertices[m])
                vertices[m].in_edges.append(vertices[i])
        # for each node, replace each edge with probability beta
        for i in range(0,n):
            for j in range(0,d_out):
                m = (i+j)%n
                if(np.random.binomial(1, beta)):
                    vertices[i].out_edges.remove(vertices[m])
                    vertices[m].in_edges.remove(vertices[i])
                    g.remove_edge(i,m)
                    vertices[i].add_new_edge(1)

    # make a scale-free graph
    else:
        # first create a cyclic graph the first 20 nodes
        for i in range(0,beta):
            for j in range(0,d_out):
                m = (i+j)%beta
                g.add_edge(i,m,1)
                vertices[i].out_edges.append(vertices[m])
                vertices[m].in_edges.append(vertices[i])
        weights = [d_out]*beta
        # rest of nodes connect to random d_out nodes weighted by edges of current nodes
        for i in range(beta, n):
            choices = list(range(i))
            ws = weights.copy()
            for j in range(d_out):
                while(True):
                    choice_i = random.choices(list(range(i-j)),ws,k=1)[0]
                    m = choices[choice_i]
                    if(len(vertices[m].in_edges)==max_in):
                        continue
                    g.add_edge(i,m,1)
                    vertices[i].out_edges.append(vertices[m])
                    vertices[m].in_edges.append(vertices[i])
                    choices.pop(choice_i)
                    ws.pop(choice_i)
                    weights[m] += 1
                    break
            weights += [d_out]

    # get initial stats
    degrees = {} # {r: [in-degree of all nodes]}
    dist_to_miners = {} # {r: [average weighted distance to miners for all nodes]}
    tries_rounds = {} # {r: [tries to make a connection for each node]}
    diameter_rounds = {} # {r: (diameter_network, diameter_miners)}
    eccentricity_rounds = {} # {r:[eccentricities]}
    ds = [] # degrees

    for i in range(0,n):
        ds.append(len(vertices[i].in_edges))
    degrees[0] = ds

    # update miner distances
    dists = []
    (D,diam_all,diam_miners, eccen_n)  = update_miners_dist(m)
    diameter_rounds[0] = (diam_all, diam_miners) #diameter of miners, whole network
    eccentricity_rounds[0] = eccen_n # eccentricity for each node
    for i in range(0,n):
        vertices[i].update_miners_dist(D.copy())
        dists.append(vertices[i].miner_dist)
    dist_to_miners[0] = dists

    # for each round, each miner drops the worst peer and makes a new random peer
    # then updates all distances
    # shuffle the order that nodes drop connections?
    for j in range(1,r+1):
        # remove worst out-edge
        for i in range(0,n):
            vertices[i].remove_worst_out()
        # add a random new edge that is not the one removed
        tries = []
        for i in range(0,n):
            tries.append(vertices[i].add_new_edge(1))
        tries_rounds[j] = tries

        #update miner distances
        dists = []
        (D,diam_all,diam_miners, eccen_n)  = update_miners_dist(m)
        diameter_rounds[j] = (diam_all, diam_miners)
        eccentricity_rounds[j] = eccen_n
        for i in range(0,n):
            vertices[i].update_miners_dist(D.copy())
            dists.append(vertices[i].miner_dist)
        dist_to_miners[j] = dists[:]

        # update degrees
        ds = []
        for i in range(0,n):
            ds.append(len(vertices[i].in_edges))
        degrees[j] = ds[:]

    dists_pdfs = {}
    dists_cdfs = {}
    # get pdf and cdf of dist to miners
    for i in range(r+1):
        dists = sorted(dist_to_miners[i])
        dists_pdfs[i],dists_cdfs[i] = get_pdf_cdf(dists)

    # get the pdf and cdf of the degrees dictionary
    degrees_pdfs = {}
    degrees_cdfs = {}
    for i in range(r+1):
        ds = sorted(degrees[i])
        degrees_pdfs[i], degrees_cdfs[i] = get_pdf_cdf(ds)

    return {'degrees': degrees, 'degrees_pdfs': degrees_pdfs, 'degrees_cdfs': degrees_cdfs,'distance_pdfs':dists_pdfs,'distance_cdfs': dists_cdfs,\
            'connection_tries':tries_rounds,'diameters':diameter_rounds, 'eccentricities': eccentricity_rounds}


# nodes, miners, d_in, d_out, rounds, whether to run multiple miner distributions
n = int(sys.argv[1])
m = int(sys.argv[2])
d_in = int(sys.argv[3])
d_out = int(sys.argv[4])
r = int(sys.argv[5])

sims = 20

if(d_in < n):
        lim = str(d_in)
else:
        lim = 'uncap'

dir = 'perigee-simulations/'

# run simulations with all miners and different models of mining power distributions
# all miners have the same mining power

"""
# small world graph (watts-strogatz graph)
beta = .5
file_name = 'perigee_uni_small_'+str(beta)+'n_'+str(n)+'_m_'+str(m)+'_d_'+str(d_out)+'_din_'+lim+'_r_'+str(r)+'.json'
print('small world')
for i in range(sims):
	s = datetime.now()
	out = run_simulation(n,m,d_out,r,d_in, 'small', beta)
	print(datetime.now()-s)
	with open(dir+file_name, 'a') as outfile:
		outfile.write(json.dumps(out)+'\n')
"""
# scale-free graph (barabasi-albert graph)
beta = 20
file_name = 'perigee_uni_sfree_'+str(beta)+'_n_'+str(n)+'_m_'+str(m)+'_d_'+str(d_out)+'_din_'+lim+'_r_'+str(r)+'.json'
print('scale free')
for i in range(sims):
        s = datetime.now()
        out = run_simulation(n,m,d_out,r,d_in, 'scalefree', beta)
        print(datetime.now()-s)
        with open(dir+file_name, 'a') as outfile:
                outfile.write(json.dumps(out)+'\n')




























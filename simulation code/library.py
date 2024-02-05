from queue import PriorityQueue

class Graph:
	def __init__(self, num_of_vertices):
		self.v = num_of_vertices
		self.edges =[ [-1 for i in range(num_of_vertices)] for j in range(num_of_vertices)]
		self.visited = []
	def add_edge(self, u, v, weight):
		self.edges[u][v] = weight
		self.edges[v][u] = weight
	def remove_edge(self,u,v):
		self.edges[u][v] = -1
		self.edges[v][u] = -1
	def dijkstra(self, start_vertex):
		#find shortest path for all miners, i.e. the first m nodes
		D = {v:float('inf') for v in range(self.v)}
		D[start_vertex] = 0
		pq = PriorityQueue()
		pq.put((0, start_vertex))
		while (not pq.empty()):
			(dist, current_vertex) = pq.get()
			self.visited.append(current_vertex)
			for neighbor in range(self.v):
				if self.edges[current_vertex][neighbor] != -1:
					if neighbor not in self.visited:
						distance = self.edges[current_vertex][neighbor]
						old_cost = D[neighbor]
						new_cost = D[current_vertex] + distance
						if new_cost < old_cost:
							pq.put((new_cost, neighbor))
							D[neighbor] = new_cost
		self.visited = []
		return D
	def dijkstra_all(self):
		# create a n by n matrix to store all minimum distances
		# set initial values to infinity
		n = self.v
		D = {}
		for i in range(n):
			D[i] = {v:float('inf') for v in range(n)}
		# let each node is the source
		for i in range(n):
			D[i][i] = 0  #distance to self is zero
			visited = []
			pq = PriorityQueue()
			pq.put((0,i))
			for j in range(0,i):
				# for j<i we've already found the minimum distance from j to i
				D[i][j] = D[j][i]
				# add their neighbors to the queue
				for neighbor in range(i+1,n):
					if self.edges[j][neighbor] != -1:
						distance = self.edges[j][neighbor]
						old_cost = D[i][neighbor]
						new_cost = D[i][j] + distance
						if new_cost < old_cost:
							pq.put((new_cost, neighbor))
							D[i][neighbor] = new_cost
				visited.append(j)
			while(not pq.empty() and len(visited) != n ):
				(dist, current_vertex) = pq.get()
				visited.append(current_vertex)
				for neighbor in range(n):
						if self.edges[current_vertex][neighbor] != -1:
							if neighbor not in self.visited:
								distance = self.edges[current_vertex][neighbor]
								old_cost = D[i][neighbor]
								new_cost = D[i][current_vertex] + distance
								if new_cost < old_cost:
									pq.put((new_cost, neighbor))
									D[i][neighbor] = new_cost
		return D

	def dijkstra_hops_all(self): # dijkstra but with edge weight set to 1
                # create a n by n matrix to store all minimum distances
                # set initial values to infinity
                n = self.v
                D = {}
                for i in range(n):
                        D[i] = {v:float('inf') for v in range(n)}
                # let each node is the source
                for i in range(n):
                        D[i][i] = 0  #distance to self is zero
                        visited = []
                        pq = PriorityQueue()
                        pq.put((0,i))
                        for j in range(0,i):
                                # for j<i we've already found the minimum distance from j to i
                                D[i][j] = D[j][i]
                                # add their neighbors to the queue
                                for neighbor in range(i+1,n):
                                        if self.edges[j][neighbor] != -1:
                                                distance = 1
                                                old_cost = D[i][neighbor]
                                                new_cost = D[i][j] + distance
                                                if new_cost < old_cost:
                                                        pq.put((new_cost, neighbor))
                                                        D[i][neighbor] = new_cost
                                visited.append(j)
                        while(not pq.empty() and len(visited) != n ):
                                (dist, current_vertex) = pq.get()
                                visited.append(current_vertex)
                                for neighbor in range(n):
                                                if self.edges[current_vertex][neighbor] != -1:
                                                        if neighbor not in self.visited:
                                                                distance = 1
                                                                old_cost = D[i][neighbor]
                                                                new_cost = D[i][current_vertex] + distance
                                                                if new_cost < old_cost:
                                                                        pq.put((new_cost, neighbor))
                                                                        D[i][neighbor] = new_cost
                return D


	def hops_dist(self, start_vertex):
		#find shortest path for all miners in number of hops
		D = {v:float('inf') for v in range(self.v)}
		D[start_vertex] = 0
		pq = PriorityQueue()
		pq.put((0, start_vertex))
		while (not pq.empty()):
			(dist, current_vertex) = pq.get()
			self.visited.append(current_vertex)
			for neighbor in range(self.v):
				if self.edges[current_vertex][neighbor] != -1:
					distance = 1
					if neighbor not in self.visited:
						old_cost = D[neighbor]
						new_cost = D[current_vertex] + distance
						if new_cost < old_cost:
							pq.put((new_cost, neighbor))
							D[neighbor] = new_cost
		self.visited = []
		return D

	def diameter_hops_dist(self, m, distance = False):
		# calculate the diameter of the network (in hops and maybe distance)
		# as well as the diameter between the miners
		max_hops = 0.0
		max_dist = 0.0
		miner_hops = 0.0
		miner_dist = 0.0
		# for each node, find the distance to the remaining nodes

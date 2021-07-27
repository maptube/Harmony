#!/usr/bin/env python3

import sys
import time
import pandas as pd
import geopandas as gpd
import numpy as np
import networkx as nx
from QUANT.PublicTransportNetwork import PublicTransportNetwork

################################################################################
# Main program
################################################################################

#load luti/zones_data.csv - contains mapping between zonei number and zone code string
dfZonesData = pd.read_csv('athens/Data_LUTI_Athens/zones_data.csv',dtype={'zone':str})
#NOTE: use the zone field as the area key
#BUT with the leading zero suppressed - this matches "NO_" field in the shapefile

print("Building Athens transport network from GTFS files")

print("initialise ptn")
ptn = PublicTransportNetwork()
#Tram = 0, Subway = 1, Rail = 2, Bus = 3, Ferry = 4, CableCar = 5, Gondola = 6, Funicular = 7
#ptn.initialiseGTFSNetworkGraph(["athens/gtfs"],{'*'}) #{1,3,800,900}?
#count = ptn.FixupWalkingNodes(500)

#print("after fixup walking nodes: ", str(ptn.graph.number_of_nodes())," vertices and ", str(ptn.graph.number_of_edges()), " edges in network.")
#ptn.writeGraphML("model-data/graph_OASA.graphml")
#ptn.writeVertexPositions("model-data/vertices_OASA.csv") #note: you also need the vertex position map, otherwise you don't know where the graphml vertices are

#finished network graph creation here
##cheat - load graph for speed
ptn.readGraphML("model-data/graph_OASA.graphml")
###

print("zone codes")
#ZoneLookup is the mapping between the zone code numbers and the MSOA code and lat/lon
#zonecentroids and boundaries has field "NO_" which is the areakey
#NOTE: NO_ is the same as CODE, but with the leading zero missing - must use NO_ to match other data
ZoneLookup = gpd.read_file('athens/Zone_boundaries/zone_centroids_WGS84.shp',dtype={'NO_':str,'CODE':str})
#need to force type on shapefile columns as dtype in above read_file doesn't work!
ZoneLookup['NO_'] = ZoneLookup['NO_'].astype(str)

#use reproject into wgs84 as that's what gtfs is in
#CentroidLookup = ptn.FindCentroids(
#    ZoneLookup,
#    'athens/Zone_boundaries/zones_boundaries_WGS84.shp', 'NO_'
#)
#save it!
#ptn.saveCentroids(CentroidLookup,"model-data/zone_centroid_lookup.csv")
## cheat - load centroids for speed
CentroidLookup = ptn.loadCentroids("model-data/zone_centroid_lookup.csv")
###

################################################################################
## Shortest paths calculations
################################################################################

#make a matrix for the results
N = 1265 #number of zones
Cij = np.zeros(N*N).reshape(N, N) #i,j are the object id
Cij.fill(-1) #use -1 as no data value

#make a lookup here of zonei->zonecode and zonei->vertexi
#dfZonesData contains:
#zonei,zone,pop_tot,floorspace,employment
#0,1001,5580,105999,303
zonei_to_zonecode = {} #mapping between zone i number and string zone code
zonei_to_vertexi = {} #mapping between zone i number and closest vertex on network to centroid
vertex_to_zonei = {} #mapping between graph vertex and closest zone i zone (DESIGNATED CENTROID ONLY)
for index, row in dfZonesData.iterrows():
    zonei = row['zonei']
    zonecode = row['zone']
    if zonecode in CentroidLookup: #otherwise, there is no centroid and no possibility of a shortrst path cost
        gvertex = CentroidLookup[zonecode]
        zonei_to_zonecode[zonei] = zonecode #might need str(zonei)
        zonei_to_vertexi[zonei] = gvertex #might need str(zonei)
        vertex_to_zonei[gvertex] = zonei
        #print(zonei,zonecode,gvertex) #test zonei->(zonecode and graph vertex)
    #endif
#end for

#now onto the shortest paths
print("running shortest paths")

#needs CUDA install to be able to do this...
#shortestPaths = ptn.TestCUDASSSP(ZoneLookup.dt,RailCentroidLookup)
#for k,v in shortestPaths.items():
#    print("SSSP",k,v)

#simpler (slower) code using networkx
start = time.clock()
#result = nx.all_pairs_dijkstra_path_length(graph,weight='weight') #58s
#result = nx.all_pairs_bellman_ford_path_length(graph,weight='weight') #200s
result = nx.all_pairs_shortest_path_length(ptn.graph) #57.5s
#result = nx.floyd_warshall(graph) #takes forever

#NOTE: result: key is the origin vertex, data contains a map of destination vertex with time
#nxgraph is returning a generator, which means we have to cycle through the data in the order
#they specify
#i is origin and j is destination
for keyi, data in result:
    #print("keyi=",keyi)
    #print("keyi=",keyi,"data=",data)
    #so key is EVERY vertex in the data and we only want selected vertices
    if keyi in vertex_to_zonei:
        #key is a vertex designated the closest to a centroid, so fill in this zone's data
        zonei = vertex_to_zonei[keyi]
        count=0
        print("SSSP: ",zonei,keyi,end='')
        for keyj in data:
            if keyj in vertex_to_zonei:
                #we've got a designated centroid vertex for the destination
                zonej = vertex_to_zonei[keyj]
                Cij[zonei,zonej]=data[keyj] #cost
                count+=1
            #endif
        #end for keyj
        print(" count=",count)
    #endif
#end for keyi


secsAPSP = time.clock()-start
print("finish")
print("graphTestNetworkX: all_pairs_dijkstra_path_length test ",secsAPSP," secs")

#and save cij here
with open("model-data/Cij_public.csv","w") as f:
    for i in range(0,N):
        for j in range(0,N):
            f.write(str(Cij[i,j]))
            if j!=N-1:
                f.write(', ')
        #end for j
        f.write('\n')
    #end for i
#end with

#finally, analyse the data
sum=0
min=sys.float_info.max
max=0
missingcount=0
datacount=0
for i in range(0,N):
    for j in range(0,N):
        value = Cij[i,j]
        if value==-1:
            missingcount+=1
        else:
            datacount+=1
            sum+=value
            if value>=max:
                max=value
            if value<=min:
                min=value
        #end else
    #end for j
#end for i
print(
    "Cij stats: mean=",sum/datacount,
    "missingcount=",missingcount,
    " datacount=",datacount,
    "max=",max,
    "min=",min
)


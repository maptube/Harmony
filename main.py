#!/usr/bin/env python3

import sys
import os.path
from os import path
import time
import pandas as pd
import geopandas as gpd
import numpy as np
import networkx as nx
from QUANT.PublicTransportNetwork import PublicTransportNetwork

################################################################################
# Global definitions
################################################################################

#constants definitions

#walking node fixup distance is used to connect network nodes that are walkable e.g.
#where there are two bus routes but none of the stops are shared, then this puts a
#walking link between network segements to make them connected - lower is better
#as it tends to add in lots of additional network links
walkFixupDistMetres = 500

#input files
zonesDataFilename = 'athens/Data_LUTI_Athens/zones_data.csv'
zoneCentroidsShapefilenameWGS84 = 'athens/Zone_boundaries/zone_centroids_WGS84.shp'
zoneBoundariesShapefilenameWGS84 = 'athens/Zone_boundaries/zones_boundaries_WGS84.shp'
athensGTFSDir = 'athens/gtfs'

#output files
modelDataDir = 'model-data' #all results go in here
graphMLFilename = 'model-data/graph_OASA.graphml' #network graph built from GTFS files
graphVertexPositionsFilename = 'model-data/vertices_OASA.csv' #lat/lon of nodes in graphml file above
zoneCentroidLookupFilename = 'model-data/zone_centroid_lookup.csv' #closest graph vertex to zone centroids
cijCostMatrixFilename = 'model-data/Cij_public.csv' #this is the travel cost matrix

################################################################################
# Main program
################################################################################

#make model-data dir if if doesn't already exist
if not path.exists(modelDataDir):
    os.mkdir(modelDataDir)

#load luti/zones_data.csv - contains mapping between zonei number and zone code string
dfZonesData = pd.read_csv(zonesDataFilename,dtype={'zone':str})
#NOTE: use the zone field as the area key
#BUT with the leading zero suppressed - this matches "NO_" field in the shapefile

print("Building Athens transport network from GTFS files")

print("initialise ptn")
ptn = PublicTransportNetwork()

################################################################################
# Make GraphML network file from GTFS data
################################################################################

if not path.exists(graphMLFilename):
    print('File %s does not exist, so creating new' % graphMLFilename)
    #Tram = 0, Subway = 1, Rail = 2, Bus = 3, Ferry = 4, CableCar = 5, Gondola = 6, Funicular = 7
    ptn.initialiseGTFSNetworkGraph([athensGTFSDir],{'*'}) #{1,3,800,900}?
    count = ptn.FixupWalkingNodes(walkFixupDistMetres)

    print('after fixup walking nodes: ', str(ptn.graph.number_of_nodes())," vertices and ", str(ptn.graph.number_of_edges()), " edges in network.")
    ptn.writeGraphML(graphMLFilename)
    ptn.writeVertexPositions(graphVertexPositionsFilename) #note: you also need the vertex position map, otherwise you don't know where the graphml vertices are

    #finished network graph creation here
else:
    #load exsting graph for speed
    print('File %s exists, skipping creation' % graphMLFilename)
    ptn.readGraphML(graphMLFilename)
###

################################################################################
# Centroid lookup
################################################################################

print('Loading zone codes from %s',zoneCentroidsShapefilenameWGS84)
#ZoneLookup is the mapping between the zone code numbers and the MSOA code and lat/lon
#zonecentroids and boundaries has field "NO_" which is the areakey
#NOTE: NO_ is the same as CODE, but with the leading zero missing - must use NO_ to match other data
ZoneLookup = gpd.read_file(zoneCentroidsShapefilenameWGS84,dtype={'NO_':str,'CODE':str})
#need to force type on shapefile columns as dtype in above read_file doesn't work!
ZoneLookup['NO_'] = ZoneLookup['NO_'].astype(str)

if not path.exists(zoneCentroidLookupFilename):
    print('File %s does not exist, creating new' % zoneCentroidLookupFilename)
    #use reproject into wgs84 as that's what gtfs is in
    CentroidLookup = ptn.FindCentroids(
        ZoneLookup,
        zoneBoundariesShapefilenameWGS84, 'NO_'
    )
    #save it!
    ptn.saveCentroids(CentroidLookup,zoneCentroidLookupFilename)
else:
    print('File %s exists, skipping creation' % zoneCentroidLookupFilename)
    ## load centroids for speed
    CentroidLookup = ptn.loadCentroids(zoneCentroidLookupFilename)
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
print("Running shortest paths")

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
    #print('keyi=',keyi)
    #print('keyi=',keyi,'data=',data)
    #so key is EVERY vertex in the data and we only want selected vertices
    if keyi in vertex_to_zonei:
        #key is a vertex designated the closest to a centroid, so fill in this zone's data
        zonei = vertex_to_zonei[keyi]
        count=0
        print('SSSP: ',zonei,keyi,end='')
        for keyj in data:
            if keyj in vertex_to_zonei:
                #we've got a designated centroid vertex for the destination
                zonej = vertex_to_zonei[keyj]
                Cij[zonei,zonej]=data[keyj] #cost
                count+=1
            #endif
        #end for keyj
        print(' count=',count) #this is how many destinations are matched for each origin
    #endif
#end for keyi


secsAPSP = time.clock()-start

print('all_pairs_dijkstra_path_length test ',secsAPSP,' secs')

#and save cij here
with open(cijCostMatrixFilename,'w') as f:
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
    'Cij stats: mean=',sum/datacount,
    'missingcount=',missingcount,
    'datacount=',datacount,
    'max=',max,
    'min=',min
)

#and we're finished...
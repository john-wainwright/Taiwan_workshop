# -*- coding: utf-8 -*-
"""
Replicating connectivity analyses from Tiwari et al. (2025)

Tiwari S, L Turnbull, J Wainwright 2025 ‘Quantification of local- and global-
   scale hydrological and sediment connectivity over grassland and shrubland 
   hillslopes’, Journal of Hydrology 655. doi: 10.1016/j.jhydrol.2025.132896
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import rasterio as rio
import rasterio.plot
import seaborn as sns

#Get DTM data
dataFolder = "../data/"

p1DEMFile = rio.open (dataFolder + "p1dem.asc")
p1DEM = p1DEMFile.read (1)
cellSize = p1DEMFile.transform [0]

fig, ax = plt.subplots ()
rio.plot.show (p1DEM, transform = p1DEMFile.transform)

def flowDirection4 (z, dx):
    '''
    calculate the cardinal flow direction from topography data using the 
    steepest descent method
    

    Parameters
    ----------
    z : numpy array, float
        topography values representing the DEM.
    dx : float
        distance between cells in the DEM (assumed square).

    Returns
    -------
    direction : numpy array, integer
        flow direction code for the corresponding cell in the DEM.
          1 = "S", 2 = "W", 3 = "N", 4 = "E"; 0 = no outflow.

    '''
    m = [1, 0, -1, 0]
    n = [0, -1, 0, 1]

    rows, columns = z.shape
    direction = np.zeros ((rows, columns), dtype = int)
    
    for i in range (1, rows - 1):
        for j in range (1, columns - 1):
            if (np.isnan (z [i, j])):
                direction [i, j] = -9999
            else:
                slpMax = -9999.9
                for k in range (4):
                    inew = i + m [k]
                    jnew = j + n [k]
                    if (inew > 0 and inew < rows and jnew > 0 and jnew < columns):
                        slope = (z [i, j] - z [inew, jnew]) / dx
                        if (slope > slpMax):
                            slpMax = slope
                            direction [i, j] = k + 1
    return direction

p1FlowDir = flowDirection4 (p1DEM, cellSize)

def createDigraphFromDEM (z, zFile, direction):
    '''
    Makes a networkx directed graph from DEM and flow direction data.
    Assumes flow direction codes for the corresponding cell in the DEM as follows:
      1 = "S", 2 = "W", 3 = "N", 4 = "E"; 0 = no outflow.
    Other values are ignored.

    Parameters
    ----------
    z : numpy array, float
        topography values representing the DEM.
    zFile : io.DatasetReader
        file containing the DEM values corresponding to z (needed to extract
                                                           geometry information).
    direction : numpy array, integer
        flow direction code for the corresponding cell in the DEM.

    Returns
    -------
    newFlowGraph : networkx directed graph
        directed graph with the directions determined by the flow direction.
    geometry : dictionary
        x and y values for each node in the directed graph allowing it to be
          plotted geographically.

    '''
    m = [1, 0, -1, 0]
    n = [0, -1, 0, 1]
    newFlowGraph = nx.DiGraph ()
    
    rows, columns = z.shape
    nodes = np.arange (0, rows * columns).tolist ()
    newFlowGraph.add_nodes_from (nodes)
    
    thisNode = 0
    geometry = {}
    for i in range (rows):
        for j in range (columns):
            geometry.update ({thisNode: zFile.xy (i, j)})
            thisNode = thisNode + 1
    
    for i in range (1, rows - 1):
        for j in range (1, columns - 1):
            if (p1FlowDir [i, j] > 0):
                thisNode = i * columns + j
                if (direction [i, j] > 0 and direction [i, j] < 5):
                    thisDir = direction [i, j] - 1
                    inew = i + m [thisDir]
                    jnew = j + n [thisDir]
                    nextNode = inew * columns + jnew
                    newFlowGraph.add_edge (thisNode, nextNode)

    return newFlowGraph, geometry

p1FlowGraph, p1FlowGraphGeom = createDigraphFromDEM (p1DEM, p1DEMFile, p1FlowDir)

fig, ax = plt.subplots ()
nx.draw_networkx_edges (p1FlowGraph, pos = p1FlowGraphGeom, ax = ax, 
                        arrows = True, arrowsize = 2, width = 1, edge_color = "b")
limits = plt.axis ('on') 
plt.axis ('scaled')
ax.autoscale (tight = True)
ax.tick_params (left = True, bottom = True, labelleft = True, labelbottom = True)
plt.savefig ("../figures/p1FlowGraph.png")
plt.show ()        

#Compare for the other plots
p2DEMFile = rio.open (dataFolder + "p2dem.asc")
p2DEM = p2DEMFile.read (1)
p3DEMFile = rio.open (dataFolder + "p3dem.asc")
p3DEM = p3DEMFile.read (1)
p4DEMFile = rio.open (dataFolder + "p4dem.asc")
p4DEM = p4DEMFile.read (1)

fig, axs = plt.subplots (1, 4)
rio.plot.show (p1DEM, transform = p1DEMFile.transform, ax = axs [0], title = "Plot 1")
rio.plot.show (p2DEM, transform = p2DEMFile.transform, ax = axs [1], title = "Plot 2")
rio.plot.show (p3DEM, transform = p3DEMFile.transform, ax = axs [2], title = "Plot 3")
rio.plot.show (p4DEM, transform = p4DEMFile.transform, ax = axs [3], title = "Plot 4")

p2FlowDir = flowDirection4 (p2DEM, cellSize)
p3FlowDir = flowDirection4 (p3DEM, cellSize)
p4FlowDir = flowDirection4 (p4DEM, cellSize)

p2FlowGraph, p2FlowGraphGeom = createDigraphFromDEM (p2DEM, p2DEMFile, p2FlowDir)
p3FlowGraph, p3FlowGraphGeom = createDigraphFromDEM (p3DEM, p3DEMFile, p3FlowDir)
p4FlowGraph, p4FlowGraphGeom = createDigraphFromDEM (p4DEM, p4DEMFile, p4FlowDir)

fig, axs = plt.subplots (1, 4)
nx.draw_networkx_edges (p1FlowGraph, pos = p1FlowGraphGeom, ax = axs [0], 
                        arrows = True, arrowsize = 2, width = 1, edge_color = "b")
axs [0].autoscale (tight = True)
axs [0].tick_params (left = True, bottom = True, labelleft = True, labelbottom = True)
axs [0].title.set_text ("Plot 1")
nx.draw_networkx_edges (p2FlowGraph, pos = p2FlowGraphGeom, ax = axs [1], 
                        arrows = True, arrowsize = 2, width = 1, edge_color = "b")
axs [1].autoscale (tight = True)
axs [1].tick_params (left = True, bottom = True, labelleft = False, labelbottom = True)
axs [1].title.set_text ("Plot 2")
nx.draw_networkx_edges (p3FlowGraph, pos = p3FlowGraphGeom, ax = axs [2], 
                        arrows = True, arrowsize = 2, width = 1, edge_color = "b")
axs [2].autoscale (tight = True)
axs [2].tick_params (left = True, bottom = True, labelleft = False, labelbottom = True)
axs [2].title.set_text ("Plot 3")
nx.draw_networkx_edges (p4FlowGraph, pos = p4FlowGraphGeom, ax = axs [3], 
                        arrows = True, arrowsize = 2, width = 1, edge_color = "b")
axs [3].autoscale (tight = True)
axs [3].tick_params (left = True, bottom = True, labelleft = False, labelbottom = True)
axs [3].title.set_text ("Plot 4")
limits = plt.axis ('on') 
plt.axis ('scaled')
plt.savefig ("../figures/PlotsFlowGraph.png")
plt.show ()        


#Structural Connectivity
#Use vegetation with threshold of 60% to stop flow pathways
p1VegFile = rio.open (dataFolder + "p1vegcover.asc")
p1Veg = p1VegFile.read (1)
p2VegFile = rio.open (dataFolder + "p2vegcover.asc")
p2Veg = p2VegFile.read (1)
p3VegFile = rio.open (dataFolder + "p3vegcover.asc")
p3Veg = p3VegFile.read (1)
p4VegFile = rio.open (dataFolder + "p4vegcover.asc")
p4Veg = p4VegFile.read (1)

fig, axs = plt.subplots (1, 4)
rio.plot.show (p1Veg, transform = p1VegFile.transform, ax = axs [0], title = "Plot 1")
rio.plot.show (p2Veg, transform = p2VegFile.transform, ax = axs [1], title = "Plot 2")
rio.plot.show (p3Veg, transform = p3VegFile.transform, ax = axs [2], title = "Plot 3")
rio.plot.show (p4Veg, transform = p4VegFile.transform, ax = axs [3], title = "Plot 4")

def createSCDigraph (z, zFile, direction, veg, vegThreshold = 60.):
    '''
    Makes a networkx directed graph from DEM and flow direction data.
    Uses vegetation data to disconnect structural flows above a given threshold.
    Assumes flow direction codes for the corresponding cell in the DEM as follows:
      1 = "S", 2 = "W", 3 = "N", 4 = "E"; 0 = no outflow.
    Other values are ignored.

    Parameters
    ----------
    z : numpy array, float
        topography values representing the DEM.
    zFile : io.DatasetReader
        file containing the DEM values corresponding to z (needed to extract
                                                           geometry information).
    direction : numpy array, integer
        flow direction code for the corresponding cell in the DEM.
    veg : numpy array, float
        vegetation cover in percent for the same cells as the topography.
    vegThreshold : float
        threshold value for vegetation cover to break the connectivity in
           flow.  Default = 60 based on percolation theory.

    Returns
    -------
    newFlowGraph : networkx directed graph
        directed graph with the directions determined by the flow direction
        but with connexions broken.
        Contains a weight parameter which is the inverse of the vegetation
          relative to the threshold = (vegThreshold - veg) / vegThreshold
    geometry : dictionary
        x and y values for each node in the directed graph allowing it to be
          plotted geographically.

    '''
    m = [1, 0, -1, 0]
    n = [0, -1, 0, 1]
    newFlowGraph = nx.DiGraph ()
    
    rows, columns = z.shape
    nodes = np.arange (0, rows * columns).tolist ()
    newFlowGraph.add_nodes_from (nodes)
    
    thisNode = 0
    geometry = {}
    for i in range (rows):
        for j in range (columns):
            geometry.update ({thisNode: zFile.xy (i, j)})
            thisNode = thisNode + 1
    
    for i in range (1, rows - 1):
        for j in range (1, columns - 1):
            if (p1FlowDir [i, j] > 0):
                thisNode = i * columns + j
                if (direction [i, j] > 0 and direction [i, j] < 5):
                    thisDir = direction [i, j] - 1
                    inew = i + m [thisDir]
                    jnew = j + n [thisDir]
                    if (veg [inew, jnew] < vegThreshold):
                        nextNode = inew * columns + jnew
                        vegWeight = (vegThreshold - veg [inew, jnew]) / vegThreshold
                        newFlowGraph.add_edge (thisNode, nextNode,
                                               weight = vegWeight)

    return newFlowGraph, geometry

p1SCGraph, p1SCGraphGeom = createSCDigraph (p1DEM, p1DEMFile, p1FlowDir,
                                            p1Veg, vegThreshold = 60.)
edgesP1SC = p1SCGraph.edges ()
weightsP1SC = [p1SCGraph [u][v]['weight'] for u,v in edgesP1SC]

fig, ax = plt.subplots ()
nx.draw_networkx_edges (p1SCGraph, pos = p1SCGraphGeom, ax = ax, 
                        arrows = True, arrowsize = 2, width = weightsP1SC, 
                        edge_color = "b")
limits = plt.axis ('on') 
plt.axis ('scaled')
ax.autoscale (tight = True)
ax.tick_params (left = True, bottom = True, labelleft = True, labelbottom = True)
plt.savefig ("../figures/p1SCGraph.png")
plt.show ()        

#Now repeat for the other plots
p2SCGraph, p2SCGraphGeom = createSCDigraph (p2DEM, p2DEMFile, p2FlowDir,
                                            p2Veg, vegThreshold = 60.)
edgesP2SC = p2SCGraph.edges ()
weightsP2SC = [p2SCGraph [u][v]['weight'] for u,v in edgesP2SC]
p3SCGraph, p3SCGraphGeom = createSCDigraph (p3DEM, p3DEMFile, p3FlowDir,
                                            p3Veg, vegThreshold = 60.)
edgesP3SC = p3SCGraph.edges ()
weightsP3SC = [p3SCGraph [u][v]['weight'] for u,v in edgesP3SC]
p4SCGraph, p4SCGraphGeom = createSCDigraph (p4DEM, p4DEMFile, p4FlowDir,
                                            p4Veg, vegThreshold = 60.)
edgesP4SC = p4SCGraph.edges ()
weightsP4SC = [p4SCGraph [u][v]['weight'] for u,v in edgesP4SC]

fig, axs = plt.subplots (1, 4)
nx.draw_networkx_edges (p1SCGraph, pos = p1SCGraphGeom, ax = axs [0], 
                        arrows = True, arrowsize = 2, width = weightsP1SC, 
                        edge_color = "b")
axs [0].autoscale (tight = True)
axs [0].tick_params (left = True, bottom = True, labelleft = True, labelbottom = True)
axs [0].title.set_text ("Plot 1")
nx.draw_networkx_edges (p2SCGraph, pos = p2SCGraphGeom, ax = axs [1], 
                        arrows = True, arrowsize = 2, width = weightsP2SC, 
                        edge_color = "b")
axs [1].autoscale (tight = True)
axs [1].tick_params (left = True, bottom = True, labelleft = False, labelbottom = True)
axs [1].title.set_text ("Plot 2")
nx.draw_networkx_edges (p3SCGraph, pos = p3SCGraphGeom, ax = axs [2], 
                        arrows = True, arrowsize = 2, width = weightsP3SC, 
                        edge_color = "b")
axs [2].autoscale (tight = True)
axs [2].tick_params (left = True, bottom = True, labelleft = False, labelbottom = True)
axs [2].title.set_text ("Plot 3")
nx.draw_networkx_edges (p4SCGraph, pos = p4SCGraphGeom, ax = axs [3], 
                        arrows = True, arrowsize = 2, width = weightsP4SC, 
                        edge_color = "b")
axs [3].autoscale (tight = True)
axs [3].tick_params (left = True, bottom = True, labelleft = False, labelbottom = True)
axs [3].title.set_text ("Plot 4")
limits = plt.axis ('on') 
plt.axis ('scaled')
plt.savefig ("../figures/PlotsSCGraph.png")
plt.show ()        

#Comparing the weights for the different plots
weightsSC = pd.DataFrame ({"weight": weightsP1SC,
                           "plot": np.repeat ("Plot 1", len (weightsP1SC))})

weightsSC = pd.concat ([weightsSC,
                       pd.DataFrame ({"weight": weightsP2SC,
                                      "plot": np.repeat ("Plot 2", 
                                                         len (weightsP2SC))})])
weightsSC = pd.concat ([weightsSC,
                       pd.DataFrame ({"weight": weightsP3SC,
                                      "plot": np.repeat ("Plot 3", 
                                                         len (weightsP3SC))})])
weightsSC = pd.concat ([weightsSC,
                       pd.DataFrame ({"weight": weightsP4SC,
                                      "plot": np.repeat ("Plot 4", 
                                                         len (weightsP4SC))})])

fig, ax = plt.subplots ()
SCBoxplot = sns.boxplot (data = weightsSC, x = "plot", y = "weight", 
                         hue = "plot", ax = ax)
SCBoxplot.set (xlabel = "", ylabel = "SC weight")
plt.show ()

weightsSC.groupby ("plot").describe ()


#Now for functional connectivity
# 5 rainfall scenarios (A=45mm, B=24mm, C=15mm, D=10mm, E=5mm)
# 3 soil moisture conditions (low=3.8%, med=10.5%, high=21.1%)
# 4 plot types (grassland=1, grass/shrub=2, shrub/grass=3, shrubland=4)

p1AHighQFile = rio.open (dataFolder + "p1_rainA_highsm_dschg.asc")
p1AHighQ = p1AHighQFile.read (1)

fig, ax = plt.subplots ()
#rio.plot.show (p1AHighQ, transform = p1AHighQFile.transform)
plt.imshow (p1AHighQ, extent = rio.plot.plotting_extent (p1DEMFile))
plt.colorbar ()

#Flow connectivity
# Network Construction Method:
# - Uses same D4 steepest descent flow routing as structural connectivity
# - Applies flow threshold (0.8 mm depth (equivalent to 0.249 L for 
#     0.5 × 0.5m cells)) to determine functional connections
# - Normalizes edge weights by maximum simulated flux (3.86 L)



def createFlowFCDigraph (z, zFile, direction, flow, 
                         flowThreshold = 0.249, maxFlow = 3.86):
    '''
    Makes a networkx directed graph from DEM and flow direction data.
    Flow connectivity occurs if flows are above a given threshold
    Assumes flow direction codes for the corresponding cell in the DEM as follows:
      1 = "S", 2 = "W", 3 = "N", 4 = "E"; 0 = no outflow.
    Other values are ignored.

    Parameters
    ----------
    z : numpy array, float
        topography values representing the DEM.
    zFile : io.DatasetReader
        file containing the DEM values corresponding to z (needed to extract
                                                           geometry information).
    direction : numpy array, integer
        flow direction code for the corresponding cell in the DEM.
    flow : numpy array, float
        flow discharge in L for the same cells as the topography.
    flowThreshold : float
        threshold value for flow to produce functional connectivity
        Default = 0.249 L.
    maxFlow : float
        upper limit of flows to be compared, used to scale weightings 
        Default = 3.86 L based on simulations in the Tiwari et al. (2025) paper

    Returns
    -------
    newFlowGraph : networkx directed graph
        directed graph with the directions determined by the flow direction
        but with connexions broken.
        Contains a weight parameter which is the flow scaled by a maximum
           flow rate (based on a range of simulations)
    geometry : dictionary
        x and y values for each node in the directed graph allowing it to be
          plotted geographically.

    '''
    m = [1, 0, -1, 0]
    n = [0, -1, 0, 1]
    newFlowGraph = nx.DiGraph ()
    
    rows, columns = z.shape
    nodes = np.arange (0, rows * columns).tolist ()
    newFlowGraph.add_nodes_from (nodes)
    
    thisNode = 0
    geometry = {}
    for i in range (rows):
        for j in range (columns):
            geometry.update ({thisNode: zFile.xy (i, j)})
            thisNode = thisNode + 1
    
    for i in range (1, rows - 1):
        for j in range (1, columns - 1):
            if (p1FlowDir [i, j] > 0):
                thisNode = i * columns + j
                if (direction [i, j] > 0 and direction [i, j] < 5 and 
                    flow [i, j] >= flowThreshold):
                    thisDir = direction [i, j] - 1
                    inew = i + m [thisDir]
                    jnew = j + n [thisDir]
                    nextNode = inew * columns + jnew
                    flowWeight = (flow [i, j] - flowThreshold) / (maxFlow - flowThreshold)
                    newFlowGraph.add_edge (thisNode, nextNode,
                                           weight = flowWeight)

    return newFlowGraph, geometry

p1AHighQFlowFCGraph, p1AHighQFlowFCGraphGeom = createFlowFCDigraph (p1DEM, p1DEMFile, 
                                                                    p1FlowDir,
                                                                    p1AHighQ, 
                                                                    flowThreshold = 0.249,
                                                                    maxFlow = 3.86)
edgesP1AHighQFlowFC = p1AHighQFlowFCGraph.edges ()
weightsP1AHighQFlowFC = [p1AHighQFlowFCGraph [u][v]['weight'] for u,v in edgesP1AHighQFlowFC]

fig, ax = plt.subplots ()
nx.draw_networkx_edges (p1AHighQFlowFCGraph, pos = p1AHighQFlowFCGraphGeom, ax = ax, 
                        arrows = True, arrowsize = 2, width = weightsP1AHighQFlowFC, 
                        edge_color = "b")
limits = plt.axis ('on') 
plt.axis ('scaled')
ax.autoscale (tight = True)
ax.tick_params (left = True, bottom = True, labelleft = True, labelbottom = True)
plt.savefig ("../figures/p1AHighQFlowFCGraph.png")
plt.show ()        

#Sediment-transport connectivity
p1AHighSedFile = rio.open (dataFolder + "p1_rainA_highsm_sedtr.asc")
p1AHighSed = p1AHighSedFile.read (1)

fig, ax = plt.subplots ()
#rio.plot.show (p1AHighSed, transform = p1AHighSedFile.transform)
plt.imshow (p1AHighSed, extent = rio.plot.plotting_extent (p1DEMFile))
plt.colorbar ()


def createSedFCDigraph (z, zFile, direction, sedTrans, 
                        sedThreshold = 1.0, maxSed = 14.18):
    '''
    Makes a networkx directed graph from DEM and flow direction data.
    Sediment connectivity occurs if transport is above a given threshold
    Assumes flow direction codes for the corresponding cell in the DEM as follows:
      1 = "S", 2 = "W", 3 = "N", 4 = "E"; 0 = no outflow.
    Other values are ignored.

    Parameters
    ----------
    z : numpy array, float
        topography values representing the DEM.
    zFile : io.DatasetReader
        file containing the DEM values corresponding to z (needed to extract
                                                           geometry information).
    direction : numpy array, integer
        flow direction code for the corresponding cell in the DEM.
    sedTrans : numpy array, float
        sedimentTransport in g for the same cells as the topography.
    flowThreshold : float
        threshold value for flow to produce functional connectivity
        Default = 1 g.
    maxSed : float
        upper limit of sediment transport to be compared, used to scale weightings 
        Default = 14.18 g based on simulations in the Tiwari et al. (2025) paper

    Returns
    -------
    newFlowGraph : networkx directed graph
        directed graph with the directions determined by the flow direction
        but with connexions broken.
        Contains a weight parameter which is the sediment transport scaled by a 
          maximum transport rate (based on a range of simulations)
    geometry : dictionary
        x and y values for each node in the directed graph allowing it to be
          plotted geographically.

    '''
    m = [1, 0, -1, 0]
    n = [0, -1, 0, 1]
    newFlowGraph = nx.DiGraph ()
    
    rows, columns = z.shape
    nodes = np.arange (0, rows * columns).tolist ()
    newFlowGraph.add_nodes_from (nodes)
    
    thisNode = 0
    geometry = {}
    for i in range (rows):
        for j in range (columns):
            geometry.update ({thisNode: zFile.xy (i, j)})
            thisNode = thisNode + 1
    
    for i in range (1, rows - 1):
        for j in range (1, columns - 1):
            if (p1FlowDir [i, j] > 0):
                thisNode = i * columns + j
                if (direction [i, j] > 0 and direction [i, j] < 5 and 
                    sedTrans [i, j] >= sedThreshold):
                    thisDir = direction [i, j] - 1
                    inew = i + m [thisDir]
                    jnew = j + n [thisDir]
                    nextNode = inew * columns + jnew
                    sedWeight = (sedTrans [i, j] - sedThreshold) / (maxSed - sedThreshold)
                    newFlowGraph.add_edge (thisNode, nextNode,
                                           weight = sedWeight)

    return newFlowGraph, geometry

p1AHighQSedFCGraph, p1AHighQSedFCGraphGeom = createSedFCDigraph (p1DEM, p1DEMFile, 
                                                                 p1FlowDir,
                                                                 p1AHighSed, 
                                                                 sedThreshold = 1.0, 
                                                                 maxSed = 14.18)
edgesP1AHighQSedFC = p1AHighQSedFCGraph.edges ()
weightsP1AHighQSedFC = [p1AHighQSedFCGraph [u][v]['weight'] for u,v in edgesP1AHighQSedFC]

fig, ax = plt.subplots ()
nx.draw_networkx_edges (p1AHighQSedFCGraph, pos = p1AHighQSedFCGraphGeom, ax = ax, 
                        arrows = True, arrowsize = 2, width = weightsP1AHighQSedFC, 
                        edge_color = "b")
limits = plt.axis ('on') 
plt.axis ('scaled')
ax.autoscale (tight = True)
ax.tick_params (left = True, bottom = True, labelleft = True, labelbottom = True)
plt.savefig ("../figures/p1AHighQSedFCGraph.png")
plt.show ()        


fig, axs = plt.subplots (1, 3)
nx.draw_networkx_edges (p1SCGraph, pos = p1SCGraphGeom, ax = axs [0], 
                        arrows = True, arrowsize = 2, width = weightsP1SC, 
                        edge_color = "b")
nx.draw_networkx_edges (p1AHighQFlowFCGraph, pos = p1AHighQFlowFCGraphGeom, ax = axs [1], 
                        arrows = True, arrowsize = 2, width = weightsP1AHighQFlowFC, 
                        edge_color = "b")
nx.draw_networkx_edges (p1AHighQSedFCGraph, pos = p1AHighQSedFCGraphGeom, ax = axs [2], 
                        arrows = True, arrowsize = 2, width = weightsP1AHighQSedFC, 
                        edge_color = "b")
limits = plt.axis ('on') 
axs [0].axis ('scaled')
axs [1].axis ('scaled')
axs [2].axis ('scaled')
axs [0].autoscale (tight = True)
axs [1].autoscale (tight = True)
axs [2].autoscale (tight = True)
axs [0].set_xlim ([0, 11])
axs [0].set_ylim ([0, 31])
axs [1].set_xlim ([0, 11])
axs [1].set_ylim ([0, 31])
axs [2].set_xlim ([0, 11])
axs [2].set_ylim ([0, 31])
axs [0].tick_params (left = True, bottom = True, labelleft = True, labelbottom = True)
axs [1].tick_params (left = True, bottom = True, labelleft = False, labelbottom = True)
axs [2].tick_params (left = True, bottom = True, labelleft = False, labelbottom = True)
axs [0].title.set_text ("SC")
axs [1].title.set_text ("FC - Flow")
axs [2].title.set_text ("FC - Sediment")
plt.savefig ("../figures/p1AHighQSCFCGraph.png")
plt.show ()        

#Comparison of same event for shrubland
p4AHighQFile = rio.open (dataFolder + "p4_rainA_highsm_dschg.asc")
p4AHighQ = p4AHighQFile.read (1)

p4AHighQFlowFCGraph, p4AHighQFlowFCGraphGeom = createFlowFCDigraph (p4DEM, p4DEMFile, 
                                                                    p4FlowDir,
                                                                    p4AHighQ, 
                                                                    flowThreshold = 0.249,
                                                                    maxFlow = 3.86)
edgesP4AHighQFlowFC = p4AHighQFlowFCGraph.edges ()
weightsP4AHighQFlowFC = [p4AHighQFlowFCGraph [u][v]['weight'] for u,v in edgesP4AHighQFlowFC]

p4AHighSedFile = rio.open (dataFolder + "p4_rainA_highsm_sedtr.asc")
p4AHighSed = p4AHighSedFile.read (1)

p4AHighQSedFCGraph, p4AHighQSedFCGraphGeom = createSedFCDigraph (p4DEM, p4DEMFile, 
                                                                 p4FlowDir,
                                                                 p4AHighSed, 
                                                                 sedThreshold = 1.0, 
                                                                 maxSed = 14.18)
edgesP4AHighQSedFC = p4AHighQSedFCGraph.edges ()
weightsP4AHighQSedFC = [p4AHighQSedFCGraph [u][v]['weight'] for u,v in edgesP4AHighQSedFC]

fig, axs = plt.subplots (2, 3)
nx.draw_networkx_edges (p1SCGraph, pos = p1SCGraphGeom, ax = axs [0, 0], 
                        arrows = True, arrowsize = 2, width = weightsP1SC, 
                        edge_color = "b")
nx.draw_networkx_edges (p1AHighQFlowFCGraph, pos = p1AHighQFlowFCGraphGeom, ax = axs [0, 1], 
                        arrows = True, arrowsize = 2, width = weightsP1AHighQFlowFC, 
                        edge_color = "b")
nx.draw_networkx_edges (p1AHighQSedFCGraph, pos = p1AHighQSedFCGraphGeom, ax = axs [0, 2], 
                        arrows = True, arrowsize = 2, width = weightsP1AHighQSedFC, 
                        edge_color = "b")
nx.draw_networkx_edges (p4SCGraph, pos = p4SCGraphGeom, ax = axs [1, 0], 
                        arrows = True, arrowsize = 2, width = weightsP4SC, 
                        edge_color = "b")
nx.draw_networkx_edges (p4AHighQFlowFCGraph, pos = p4AHighQFlowFCGraphGeom, ax = axs [1, 1], 
                        arrows = True, arrowsize = 2, width = weightsP4AHighQFlowFC, 
                        edge_color = "b")
nx.draw_networkx_edges (p4AHighQSedFCGraph, pos = p4AHighQSedFCGraphGeom, ax = axs [1, 2], 
                        arrows = True, arrowsize = 2, width = weightsP4AHighQSedFC, 
                        edge_color = "b")
limits = plt.axis ('on') 
axs [0, 0].axis ('scaled')
axs [0, 1].axis ('scaled')
axs [0, 2].axis ('scaled')
axs [1, 0].axis ('scaled')
axs [1, 1].axis ('scaled')
axs [1, 2].axis ('scaled')
axs [0, 0].autoscale (tight = True)
axs [0, 1].autoscale (tight = True)
axs [0, 2].autoscale (tight = True)
axs [1, 0].autoscale (tight = True)
axs [1, 1].autoscale (tight = True)
axs [1, 2].autoscale (tight = True)
axs [0, 0].set_xlim ([0, 11])
axs [0, 0].set_ylim ([0, 31])
axs [0, 1].set_xlim ([0, 11])
axs [0, 1].set_ylim ([0, 31])
axs [0, 2].set_xlim ([0, 11])
axs [0, 2].set_ylim ([0, 31])
axs [1, 0].set_xlim ([0, 11])
axs [1, 0].set_ylim ([0, 31])
axs [1, 1].set_xlim ([0, 11])
axs [1, 1].set_ylim ([0, 31])
axs [1, 2].set_xlim ([0, 11])
axs [1, 2].set_ylim ([0, 31])
axs [0, 0].tick_params (left = True, bottom = True, labelleft = True, labelbottom = True)
axs [0, 1].tick_params (left = True, bottom = True, labelleft = False, labelbottom = True)
axs [0, 2].tick_params (left = True, bottom = True, labelleft = False, labelbottom = True)
axs [1, 0].tick_params (left = True, bottom = True, labelleft = True, labelbottom = True)
axs [1, 1].tick_params (left = True, bottom = True, labelleft = False, labelbottom = True)
axs [1, 2].tick_params (left = True, bottom = True, labelleft = False, labelbottom = True)
axs [0, 0].title.set_text ("SC")
axs [0, 1].title.set_text ("FC - Flow")
axs [0, 2].title.set_text ("FC - Sediment")
axs [0, 0].set_ylabel ("Grassland")
axs [1, 0].set_ylabel ("Shrubland")
plt.savefig ("../figures/p1p4AHighQSCFCGraph.png")
plt.show ()        

weightsSCFCp1 = pd.DataFrame ({"weight": weightsP1SC,
                               "connectivity": np.repeat ("SC", len (weightsP1SC))})

weightsSCFCp1 = pd.concat ([weightsSCFCp1,
                            pd.DataFrame ({"weight": weightsP1AHighQFlowFC,
                                           "connectivity": np.repeat ("FC - Flow", 
                                                         len (weightsP1AHighQFlowFC))})])
weightsSCFCp1 = pd.concat ([weightsSCFCp1,
                            pd.DataFrame ({"weight": weightsP1AHighQSedFC,
                                           "connectivity": np.repeat ("FC - Sed", 
                                                         len (weightsP1AHighQSedFC))})])
weightsSCFCp4 = pd.DataFrame ({"weight": weightsP4SC,
                               "connectivity": np.repeat ("SC", len (weightsP4SC))})

weightsSCFCp4 = pd.concat ([weightsSCFCp4,
                            pd.DataFrame ({"weight": weightsP4AHighQFlowFC,
                                           "connectivity": np.repeat ("FC - Flow", 
                                                         len (weightsP4AHighQFlowFC))})])
weightsSCFCp4 = pd.concat ([weightsSCFCp4,
                            pd.DataFrame ({"weight": weightsP4AHighQSedFC,
                                           "connectivity": np.repeat ("FC - Sed", 
                                                         len (weightsP4AHighQSedFC))})])

fig, axs = plt.subplots (2, 1)
SCBoxplot1 = sns.boxplot (data = weightsSCFCp1, x = "connectivity", y = "weight", 
                         hue = "connectivity", ax = axs [0])
SCBoxplot1.set (xlabel = "", ylabel = "grassland weights")
SCBoxplot1 = sns.boxplot (data = weightsSCFCp4, x = "connectivity", y = "weight", 
                         hue = "connectivity", ax = axs [1])
SCBoxplot1.set (xlabel = "", ylabel = "shrubland weights")
plt.show ()

weightsSCFCp1.groupby ("connectivity").describe ()
weightsSCFCp4.groupby ("connectivity").describe ()


#Global Metrics
##Centralization Degree
p1SCDegCent = nx.degree_centrality (p1SCGraph)
p1SCDegCentNonZero = pd.Series (p1SCDegCent)[pd.Series (p1SCDegCent) > 0]
n = p1SCDegCentNonZero.count ()
p1SCCentDeg = sum (max (p1SCDegCentNonZero) - p1SCDegCentNonZero) / (n - 2)

p2SCDegCent = nx.degree_centrality (p2SCGraph)
p2SCDegCentNonZero = pd.Series (p2SCDegCent)[pd.Series (p2SCDegCent) > 0]
n = p2SCDegCentNonZero.count ()
p2SCCentDeg = sum (max (p2SCDegCentNonZero) - p2SCDegCentNonZero) / (n - 2)

p3SCDegCent = nx.degree_centrality (p3SCGraph)
p3SCDegCentNonZero = pd.Series (p3SCDegCent)[pd.Series (p3SCDegCent) > 0]
n = p3SCDegCentNonZero.count ()
p3SCCentDeg = sum (max (p3SCDegCentNonZero) - p3SCDegCentNonZero) / (n - 2)

p4SCDegCent = nx.degree_centrality (p4SCGraph)
p4SCDegCentNonZero = pd.Series (p4SCDegCent)[pd.Series (p4SCDegCent) > 0]
n = p4SCDegCentNonZero.count ()
p4SCCentDeg = sum (max (p4SCDegCentNonZero) - p4SCDegCentNonZero) / (n - 2)
print (p1SCCentDeg, p2SCCentDeg, p3SCCentDeg, p4SCCentDeg)

p1SCGE = pd.Series (nx.closeness_centrality (p1SCGraph)).sum ()
p2SCGE = pd.Series (nx.closeness_centrality (p2SCGraph)).sum ()
p3SCGE = pd.Series (nx.closeness_centrality (p3SCGraph)).sum ()
p4SCGE = pd.Series (nx.closeness_centrality (p4SCGraph)).sum ()
print (p1SCGE, p2SCGE, p3SCGE, p4SCGE)

p1SCAC = nx.degree_pearson_correlation_coefficient (p1SCGraph, weight = "weight")
p2SCAC = nx.degree_pearson_correlation_coefficient (p2SCGraph, weight = "weight")
p3SCAC = nx.degree_pearson_correlation_coefficient (p3SCGraph, weight = "weight")
p4SCAC = nx.degree_pearson_correlation_coefficient (p4SCGraph, weight = "weight")
print (p1SCAC, p2SCAC, p3SCAC, p4SCAC)



#Local Metrics
##Weighted Mean Length of Connected Pathway
p1SCLOCOP = [nx.shortest_path_length (p1SCGraph, target = i, weight = "weight") 
             for i in p1SCGraph]
p1SCMeanLOCOP = [np.array (list (p1SCLOCOP [i].values ())).mean () 
                 for i in range (len (p1SCLOCOP))]
p1SCWLOCOP = np.array (p1SCMeanLOCOP).mean ()
#nx.set_node_attributes (p1SCGraph, dict (enumerate (np.array (p1SCMeanLOCOP).flatten (), 0)), "WLOCOP")


p4SCLOCOP = [nx.shortest_path_length (p4SCGraph, target = i, weight = "weight") 
             for i in p4SCGraph]
p4SCMeanLOCOP = [np.array (list (p4SCLOCOP [i].values ())).mean () 
                 for i in range (len (p4SCLOCOP))]
p4SCWLOCOP = np.array (p4SCMeanLOCOP).mean ()
#nx.set_node_attributes (p1SCGraph, dict (enumerate (np.array (p4SCMeanLOCOP).flatten (), 0)), "WLOCOP")

print (p1SCWLOCOP, p4SCWLOCOP)

fig, axs = plt.subplots (1, 2)
nx.draw_networkx_edges (p1SCGraph, pos = p1SCGraphGeom, ax = axs [0], 
                        arrows = True, arrowsize = 2, width = 1, 
                        edge_color = "b")
nx.draw_networkx_nodes (p1SCGraph, pos = p1SCGraphGeom, ax = axs [0], node_size = p1SCMeanLOCOP, node_color = "k")
nx.draw_networkx_edges (p4SCGraph, pos = p4SCGraphGeom, ax = axs [1], 
                        arrows = True, arrowsize = 2, width = 1, 
                        edge_color = "b")
nx.draw_networkx_nodes (p4SCGraph, pos = p4SCGraphGeom, ax = axs [1], node_size = p4SCMeanLOCOP, node_color = "k")
limits = plt.axis ('on') 
axs [0].axis ('scaled')
axs [1].axis ('scaled')
axs [0].autoscale (tight = True)
axs [1].autoscale (tight = True)
axs [0].set_xlim ([0, 11])
axs [0].set_ylim ([0, 31])
axs [1].set_xlim ([0, 11])
axs [1].set_ylim ([0, 31])
axs [0].tick_params (left = True, bottom = True, labelleft = True, labelbottom = True)
axs [1].tick_params (left = True, bottom = True, labelleft = False, labelbottom = True)
axs [0].title.set_text ("Plot 1")
axs [1].title.set_text ("Plot 4")
plt.savefig ("../figures/p1p4SCMeanLocopGraph.png")
plt.show ()     

##Betweenness centrality
p1SCBC = nx.betweenness_centrality (p1SCGraph, weight = "weight")
p1SCMeanBC = pd.Series (p1SCBC).mean ()
p4SCBC = nx.betweenness_centrality (p4SCGraph, weight = "weight")
p4SCMeanBC = pd.Series (p4SCBC).mean ()
print (p1SCMeanBC, p4SCMeanBC)

sizeScale = 20. / max (pd.Series (p1SCBC).max (), pd.Series (p4SCBC).max ())
fig, axs = plt.subplots (1, 2)
nx.draw_networkx_edges (p1SCGraph, pos = p1SCGraphGeom, ax = axs [0], 
                        arrows = True, arrowsize = 2, width = 1, 
                        edge_color = "b")
nx.draw_networkx_nodes (p1SCGraph, pos = p1SCGraphGeom, ax = axs [0], node_size = pd.Series (p1SCBC) * sizeScale, node_color = "k")
nx.draw_networkx_edges (p4SCGraph, pos = p4SCGraphGeom, ax = axs [1], 
                        arrows = True, arrowsize = 2, width = 1, 
                        edge_color = "b")
nx.draw_networkx_nodes (p4SCGraph, pos = p4SCGraphGeom, ax = axs [1], node_size = pd.Series (p4SCBC) * sizeScale, node_color = "k")
limits = plt.axis ('on') 
axs [0].axis ('scaled')
axs [1].axis ('scaled')
axs [0].autoscale (tight = True)
axs [1].autoscale (tight = True)
axs [0].set_xlim ([0, 11])
axs [0].set_ylim ([0, 31])
axs [1].set_xlim ([0, 11])
axs [1].set_ylim ([0, 31])
axs [0].tick_params (left = True, bottom = True, labelleft = True, labelbottom = True)
axs [1].tick_params (left = True, bottom = True, labelleft = False, labelbottom = True)
axs [0].title.set_text ("Plot 1")
axs [1].title.set_text ("Plot 4")
plt.savefig ("../figures/p1p4SCBCGraph.png")
plt.show ()     

##Relative node efficiency
p1SCGE = nx.global_efficiency (p1SCGraph.to_undirected ())
p1SCRNE = []
for i in range (len (p1SCGraph)):
    tempGraph = p1SCGraph.to_undirected ()
    tempGraph.remove_node (i)
    thisRNE = 100. * (p1SCGE - nx.global_efficiency (tempGraph)) / p1SCGE
    p1SCRNE.append (thisRNE)

p4SCGE = nx.global_efficiency (p4SCGraph.to_undirected ())
p4SCRNE = []
for i in range (len (p4SCGraph)):
    tempGraph = p4SCGraph.to_undirected ()
    tempGraph.remove_node (i)
    thisRNE = 100. * (p4SCGE - nx.global_efficiency (tempGraph)) / p4SCGE
    p4SCRNE.append (thisRNE)

sizeScale1 = -min (pd.Series (p1SCRNE).min (), pd.Series (p4SCRNE).min ())
sizeScale2 = 20. / (sizeScale1 + max (pd.Series (p1SCRNE).max (), pd.Series (p4SCRNE).max ()))
fig, axs = plt.subplots (1, 2)
nx.draw_networkx_edges (p1SCGraph, pos = p1SCGraphGeom, ax = axs [0], 
                        arrows = True, arrowsize = 2, width = 1, 
                        edge_color = "b")
nx.draw_networkx_nodes (p1SCGraph, pos = p1SCGraphGeom, ax = axs [0], node_size = sizeScale1 + pd.Series (p1SCRNE) * sizeScale2, node_color = "k")
nx.draw_networkx_edges (p4SCGraph, pos = p4SCGraphGeom, ax = axs [1], 
                        arrows = True, arrowsize = 2, width = 1, 
                        edge_color = "b")
nx.draw_networkx_nodes (p4SCGraph, pos = p4SCGraphGeom, ax = axs [1], node_size = sizeScale1 + pd.Series (p4SCRNE) * sizeScale2, node_color = "k")
limits = plt.axis ('on') 
axs [0].axis ('scaled')
axs [1].axis ('scaled')
axs [0].autoscale (tight = True)
axs [1].autoscale (tight = True)
axs [0].set_xlim ([0, 11])
axs [0].set_ylim ([0, 31])
axs [1].set_xlim ([0, 11])
axs [1].set_ylim ([0, 31])
axs [0].tick_params (left = True, bottom = True, labelleft = True, labelbottom = True)
axs [1].tick_params (left = True, bottom = True, labelleft = False, labelbottom = True)
axs [0].title.set_text ("Plot 1")
axs [1].title.set_text ("Plot 4")
plt.savefig ("../figures/p1p4SCRNEGraph.png")
plt.show ()     


##Reading in all the Functional Connectivity Files

plotLabel = ["1", "2", "3", "4"]
rainfallLabel = ["A", "B", "C", "D", "E"]
soilMoistureLabel = ["_lowsm", "_medsm", "_highsm"]
demData = [p1DEM, p2DEM, p3DEM, p4DEM]
demDataFile = [p1DEMFile, p2DEMFile, p3DEMFile, p4DEMFile]
flowDirData = [p1FlowDir, p2FlowDir, p3FlowDir, p4FlowDir]
plotOutput = []
rainOutput = []
smOutput = []
flowOutput = []
centDegOutput = []

for plot in range (len (plotLabel)):
    for rain in range (len (rainfallLabel)):
        for soilMoisture in range (len (soilMoistureLabel)):
            filename = dataFolder + "p" + plotLabel [plot] + "_rain" + rainfallLabel [rain] + soilMoistureLabel [soilMoisture] + "_dschg.asc"
            with (rio.open (filename)) as dataFile:
                data = dataFile.read (1)
            thisFCGraph, thisFCGraphGeom = createFlowFCDigraph (demData [plot], 
                                                                demDataFile [plot],
                                                                flowDirData [plot],
                                                                data,
                                                                flowThreshold = 0.249,
                                                                maxFlow = 3.86)
            theseEdges = thisFCGraph.edges ()
            theseWeights = [thisFCGraph [u][v]['weight'] for u,v in theseEdges]
            thisDegCent = nx.degree_centrality (thisFCGraph)
            thisDegCentNonZero = pd.Series (thisDegCent)[pd.Series (thisDegCent) > 0]
            n = thisDegCentNonZero.count ()
            if (n > 2):
                thisCentDeg = sum (max (thisDegCentNonZero) - thisDegCentNonZero) / (n - 2)
            else:
                thisCentDeg = 0.
            
            plotOutput.append ("Plot " + plotLabel [plot]) 
            rainOutput.append ("Rain " + rainfallLabel [rain]) 
            smOutput.append (soilMoistureLabel [soilMoisture] [1:])
            flowOutput.append (data.sum ())
            centDegOutput.append (thisCentDeg)

flowFCOutput = pd.DataFrame ({"Plot": plotOutput, 
                              "Rain": rainOutput, 
                              "SM": smOutput, 
                              "Total flow": flowOutput,
                              "Centralization Degree": centDegOutput})

fig, ax = plt.subplots ()
gridPlot = sns.FacetGrid (data = flowFCOutput,
                          col = "Plot", col_wrap = 4)
gridPlot.map_dataframe (sns.scatterplot, x = "Total flow",
                        y = "Centralization Degree", hue = "Rain", style = "SM")
gridPlot.add_legend ()
plt.savefig ("../figures/FlowFCCDGraph.png")
plt.show ()     

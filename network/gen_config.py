#!/usr/bin/env python2.7
# gen_config.py

"""
    gen_config.py pattern.pat


        Generate a network based on a given configuration. This may be specified in a .pat file, or directed as arguments of the script.


    Options
    -------
        topology : string, default : hexagonal
            Topology of the network can be chosen amongst : triangular, hexagonal, chain (1d)

        outdir : string, default : ~/git/cavitation/network/network
            Directory path where to write the configuration files

        subdir : string, default : config
            Subdirectory

        nsim : int, default : 1
            Number of simulations

        tpl : str, default : test.ini
            Template used for the generation of files.

        nlayers : int, default : 1
            Number of layers for the network. Minimum is 1.

        noisy : boolean, default : False
            True : noisy lattice from regular lattive
            False : regular lattice

        seed : int, default : None
            If a integer is given, it corresponds to the input seed for random generation.

        vol_avg : float, default : 1.
            Average of the lumen volume distribution

        vol_std : float, default : 0.1
            Standard deviation of the lumen volume distribution

        lumen_pos_avg : float, default : 0
            Average of the lumen position distribution

        lumen_pos_std : float, default : 0.05
            Standard deviation of the lumen position distribution

        center : tuple, default : 20. 20.

        gamma_border : float, optional, default : 1.0
    
        gamma_c_border : float, optional, default : 1.0
    
        chimeras : boolean

        gamma_chimera : float, optional, default : 1.0

        nbicellular : int, optional
            If specified, nbicellular < 1 is the number of lumens at bicellular contacts. 
            By default, the bicellular lumen are evenly spaced. Bi- and multicellular lumens
            have the same initial area distribution.

        pbc : boolean, optional, default : False
            Include periodic boundary conditions. Implemented only for 1d chains.

        show : boolean, default : False
            Show the result. Not available on cluster.
    
    Examples
    --------

    1. Direct arguments
        ./gen_config.py outdir=network2 tpl=test.ini noisy=True seed=1 vol_avg=2.1 vol_std=0.4 lumen_pos_avg=1. lumen_pos_std=0.1

        Generate one configuration in the folder network2, with noisy lattice and given parameters. 

    2.
        ./gen_config config.pat
        
        Generate the configuration from the config.pat file.

    Functions
    ---------
    read_config
        ## Chain networks
    chain_network
    chain_config
        ## Triangular networks
    reset
    giveN
    find_vertex_index
    find_edge
    vertexhere
    sort_edge_list
    calc_connectivity
    calc_border_list
    close_borders
    neighbours
    triangular_config
        ## Hexagonal networks
    give_radius_hexagonal
    hexa_config
    gen_bicellular_graph

    main

    Imports
    -------
    Libraries : numpy (np), networkx (nx), sys, os, shutil, config_lib (*), 
                matplotlib.pyplot (plt), matplotlib.patches (ptc)
    Functions : scipy.spatial.Voronoi, scipy.spatial.voronoi_plot_2d, matplotlib.cm
    
    
    Created September 2018 by M. Le Verge--Serandour
    Last update : Feb. 2019
"""

import numpy as np
import networkx as nx
import sys, os
import shutil
from scipy.spatial import Voronoi, voronoi_plot_2d

from config_lib import *

matplotlib_ON = False

try :
    import matplotlib.pyplot as plt
    import matplotlib.patches as ptc
    from matplotlib import cm
    matplotlib_ON = True
    
except :
    pass

tpl= 'test.ini'

# =================================================================================================
# ==================================  READ CONFIGURATION ==========================================
# =================================================================================================

def read_config(conf_name) :
    config = configparser.ConfigParser()
    config.read(conf_name)
    return config

# =================================================================================================
# ==============================  CHAIN (1D) CONFIGURATION ========================================
# =================================================================================================

def chain_network(N, a=1., pbc=False) :
    """
    chain_network(N, a=1., pbc=False)
    
        Generate a chain network
    
        Inputs
        ------
        N : int
            Number of lumens in the chain
        a : float, optional, default : 1.0
            Distance between lumens
        pbc : boolean, optional, default : False
            Include periodic boundary conditions (N+1 = 1) if True
    
        Returns
        -------
        c : array
            Coordinates of the graph
        v : array
            Vertices of the graph
        e : array
            Edges of the graph
    """
    v = np.arange(N, dtype=int)
    c = np.zeros((N, 2))
    c[:, 0] = a*np.arange(0, N, dtype=float)
    
    e = np.zeros((N-1, 2), dtype=int)
    
    for i in range(N-1) :
        e[i] = np.array([[i, i+1]])
    if pbc == True :
        e = np.insert(e, obj =1, values=np.array([[0, N-1]]), axis=0)
    return c, v, e

def chain_config(nlayers, noisy, lumen_pos_avg, lumen_pos_std, show, pbc, center) :
    """
    chain_config(nlayers, noisy, seed, lumen_pos_avg, lumen_pos_std, show, pbc, center)
    
        Generate a chain configuration
    
       Inputs
       ------
       nlayers : int
           Number of layers of the graph. Minimum is one, no maximum.
   
       noisy : float
           Makes the output graph a noisy graph if True, a regular if False.
   
       lumen_pos_avg : float, optional
           If specified, correspond to the mean deviation of the noisy lattice with respect to the regular lattice
       lumen_pos_std : float, optional
           If specified, correspond to the standard deviation of the noisy lattice with respect to the regular lattice
   
       show : boolean, optional
           If specified True, shows the graph.
   
       pbc : boolean, optional, not implemented
           If specified, makes the graph periodic
       
       center : float array, optional, default : np.array([0., 0.])
           If specified, sets the center of the graph to the given position.
    
        Returns
        -------
        coord : array
            Array of coordinates of the graph.
    
        vertices : array
            Array of vertices of the graph.
    
        edge_list : array
            Array of vertices of the graph.
    
        resistances : array
            Array of resistance of the graph's edges.
    
        borders : set
            Set of vertices belonging to the broders of the graph.
    """
    
    if nlayers < 1 :
        print 'Not enough vertices !'
        return;
    
    coordinates, vertices, edges = chain_network(nlayers, pbc=pbc)
    
    resistances = resistance_edge(coord=coordinates, edge_list=edges)
    
    if pbc :
        borders = []
        resistances[1] = resistances[0] # the second link is the periodic one, so its length is set to be the same as the others
    else :
        borders = [np.min(vertices), np.max(vertices)]
        
    if show :
        plt.figure(figsize = (10, 10))

        for v in vertices :
            plt.scatter(coordinates[v, 0], coordinates[v, 1])
            #plt.text(coord[v, 0]+0.05, coord[v, 1]+0.02, s = str(v))
        for e in edges :
            plt.plot( (coordinates[e[0]][0], coordinates[e[1]][0]), (coordinates[e[0]][1], coordinates[e[1]][1]), color = 'k', alpha = 0.4)
        
    
    return coordinates, vertices, edges, resistances, borders

# =================================================================================================
# ==============================  TRIANGULAR CONFIGURATION ========================================
# =================================================================================================
def reset(center) :
    """
    reset(center)
    
        Reset or initalize the graph
    
        Inputs
        ------
        center : array of floats
            Coordinates of the center of the graph
        Returns
        -------
        coord : array
            Array of coordinates
        vertex : array
            Array of vertices
        edge_list : array
            Array of edges
    """
    coord = center
    vertex = np.array([0])
    edge_list = np.array([[0, 0]], dtype=int)
    return coord, vertex, edge_list

def giveN(nb_layers) :
    """
    giveN(nb_layers)
    
        Give a number to generate a graph with nb_layers layers
    """
    k=0
    for i in range(nb_layers) :
        k+=i
    return k
    
def find_vertex_index(x, y, coord, eps=0.1) :
    """
    find_vertex_index(x, y, coord, eps=0.1)
    
        Find the index of the vertex at given (x,y) coordinates
        
        Inputs
        ------
        x, y : floats
            coordinates of position
        coord : array
            Array of coordinates of vertices
        eps : float, optional, default : 0.1
            Tolerance to find the given vertex.
        Returns
        -------
        v : index of the vertex
    """
    try  :
        v = np.where(( np.abs(coord - np.array([[x, y]])) < eps ).all(axis=1))[0][0]
        return v
    except :
        pass

def find_edge(index1, index2, edge_list) :
    """
    find_edge(index1, index2, edge_list)
    
        Search in the edge_list if there is already an edge [index1, index2]
    
        Input
        -----
        index1 : int
            index of the first vertex
        index2 : int
            index of the second vertex
    
        Returns
        -------
        e : array
            If exists, index of the edge [index1, index2]            
    """
    if index1 > index2 :
        index1, index2 = index2, index1
    try :
        e = np.where((edge_list == (index1, index2)).all(axis=1))[0][0]
        return e
    except :
        return None

def vertexhere(x, y, coord, eps) :
    """
    vertexhere(x, y, coord, eps)
    
        Search if there is a vertex at x, y coordinates
    
        Inputs
        ------
        x : float
    
        y : float
    
        coord : array
            Array of coordinates of vertices of the graph
    
        eps : float
            Tolerance of the search
    
        Returns
        -------
        True if there is a vertex
        False otherwise
    """
    u = np.array([x, y])
    if len(np.shape(coord)) == 1 :
        if np.sum(coord - u) == 0. :
            return True
        else : return False
    
    else :
        if np.sum(np.linalg.norm(coord - u, axis = 1) < eps) :
            # there is a vertex
            return True
        else :
            # there is no vertex
            return False

def sort_edge_list(edge_list) :
    """
    sort_edge_list(edge_list)
    
        Sort the edge list from top to bottom, left to right.
    
        Inputs
        ------
        edge_list : array
            Array of edges
    
        Returns
        -------
        e : array
            Sorted edge array
    """
    e = np.zeros(np.shape(edge_list), dtype = int)
    for i in range(len(edge_list)) :
        e[i] = np.array([np.min(edge_list[i]), np.max(edge_list[i])])
    return e

def calc_connectivity(edge_list, coord) :
    """
    calc_connectivity(edge_list, coord)
    
        Calculates the connectivity of each vertex in the graph
    
        Inputs
        ------
        edge_list : array
            Array of edges of the graph
        coord : array
            Array of coordinates of the graph
    
        Returns
        -------
        connectivity_list : array
            Array of connectivities of each vertex
    """
    connectivity_list = np.zeros(len(coord), dtype=int)
    for v in range(len(coord)) :
        connectivity_list[v] = np.sum(edge_list == v)
    return connectivity_list

def calc_border_list(edge_list, coord) :
    """
    calc_border_list(edge_list, coord)
    
        Find the borders of the triangular graph by finding the vertices with connectivity lower than 6
        
        Input
        -----
        edge_list : array
            Array of edges of the graph
        coord : array
            Array of coordinates of the graph
        
        Returns
        -------
        borders : array
            Array of vertices belonging to the borders of the graph.
    """
    # find the borders : the lumens with less than 6 neighbors for triangular lattice
    c = calc_connectivity(edge_list, coord)
    return np.where(c < 6)[0]

def close_borders(edge_list, coord) :
    """
    close_borders(edge_list, coord)
    
        For a triangular graph, closes the borders, ie : makes the edges at the border in case they miss.
    
        Inputs
        ------
        edge_list : array
            Array of edges
        coord : array
            Array of coordinates
        Returns
        -------
        edge_list : array
            Updated edge_list.
    """
    border = calc_border_list(edge_list, coord)
    for b in border :
        # find the neighbors (vertices at less than 1.1 in distance)
        ngb = np.where(np.linalg.norm(coord[b] - coord, axis = 1) <= 1.1)[0]    # NB : inferior to 1.1 for numerical errors
        for n in ngb : # for each neighbor n of the border vertex b
            if find_edge(b, n, edge_list) == None and n != b : # if edge (n, b) not in graph, then add it
                edge_list = np.append(edge_list, np.array([[min(n, b), max(n, b)]], dtype = int), axis = 0)
    return edge_list
    
def neighbours(index, xp, yp, coord, vertex, edge_list) :
    """
    neighbours(index, xp, yp, coord, vertex, edge_list)
        
        For a given vertex with given coordinates, generates new neighbours in order to make a triangular graph
        
        Inputs
        ------
        index : int
            Index of the current vert
        xp : float
            x-axis position of the current vertice where to generate new neighbours
        yp : float
            y-axis position of the current vertice where to generate new neighbours
        coord : array
            Current list of coordinates
        vertex : array
            Current list of vertices
        edge_list :
            Current list of edges
        Returns
        -------
        coord : array
            Updated coordinates array
        vertex  : array
            Updated vertices array
        edge_list : array
            Updated edges array
    """
    eps = 0.1
    p1 = np.array([xp + 1., yp])
    p2 = np.array([xp + .5, yp + np.sqrt(3.)/2.])
    p3 = np.array([xp - .5, yp + np.sqrt(3.)/2.])
    p4 = np.array([xp - 1., yp])
    p5 = np.array([xp - .5, yp - np.sqrt(3.)/2.])
    p6 = np.array([xp + .5, yp - np.sqrt(3.)/2.])
    p_list = [p1, p2, p3, p4, p5, p6]
    
    
    for i in range(len(p_list)) :
        # if there is no vertex at the position : generate and connect it to the lumen        
        if not vertexhere(p_list[i][0], p_list[i][1], coord, eps) :
            new_index = np.max(vertex)+1
            vertex = np.append(vertex, new_index)
            coord = np.append(coord, np.array([[p_list[i][0], p_list[i][1]]]), axis = 0 )
            edge_list = np.append(edge_list, np.array([[index, new_index]], dtype=int), axis = 0)            
            
        # else : there is already a lumen, so check if connected
        else :
            v = find_vertex_index(p_list[i][0], p_list[i][1], coord, eps)
            e = find_edge(index, v, edge_list)
            if e == None :
                if index > v :
                    i1, i2 = v, index
                else :
                    i1, i2 = index, v
                edge_list = np.append(edge_list, np.array([[i1, i2]], dtype=int), axis = 0)
    return coord, vertex, edge_list
    
def triangular_config(nlayers, noisy, lumen_pos_avg, lumen_pos_std, show, pbc, center) :
    """
    triangular_config(nlayers, noisy, seed, lumen_pos_avg, lumen_pos_std, show, pbc, center)
        Generate a triangular graph
    
        Inputs
        ------
        nlayers : int
            Number of layers of the graph. Minimum is one, no maximum.
    
        noisy : float
            Makes the output graph a noisy graph if True, a regular if False.
    
        lumen_pos_avg : float, optional
            If specified, correspond to the mean deviation of the noisy lattice with respect to the regular lattice
        lumen_pos_std : float, optional
            If specified, correspond to the standard deviation of the noisy lattice with respect to the regular lattice
    
        show : boolean, optional
            If specified True, shows the graph.
    
        pbc : boolean, optional, not implemented
            If specified, makes the graph periodic
        
        center : float array, optional
            If specified, sets the center of the graph to the given position.
            
        Returns
        -------
        coord : array
            Array of coordinates of the graph.
    
        vertices : array
            Array of vertices of the graph.
    
        edge_list : array
            Array of vertices of the graph.
    
        resistances : array
            Array of resistance of the graph's edges.
    
        borders : set
            Set of vertices belonging to the broders of the graph.
    """
    if nlayers > 0 :    
        N = giveN(nlayers)
    else :
        print 'Error : you have less than 1 layer !'
        return;
    
    if pbc :
        print 'Periodic boundary conditions are not implemented yet !'
    
    # initialize the graph
    coord, vertices, edge_list = reset(center)
    
    coord, vertices, edge_list = neighbours(0, 0., 0., coord, vertices, edge_list)
    
    edge_list = np.delete(edge_list, 0, 0)
    
    # run the graph
    for i in range(N*6 + 1) :
        coord, vertices, edge_list = neighbours(i, coord[i, 0], coord[i, 1], coord, vertices, edge_list)

    # close the borders
    edge_list = close_borders(edge_list, coord)

    # sort the list
    edge_list = sort_edge_list(edge_list)
    
    # borders
    
    borders = calc_border_list(edge_list, coord)
    v = []
    
    for i in range(len(vertices)) :
        if i in borders :
            v += [[i, 1]]
        else :
            v += [[i, 0]]
    
    vertices = np.array(v)
    
    # if noisy
    if noisy :
        coord = coord + np.random.normal(loc=lumen_pos_avg, scale=lumen_pos_std, size=(len(coord), 2))
        
    # resistances
    
    resistances = resistance_edge(coord, edge_list)
    
    if show :
        plt.figure(figsize = (10, 10))

        for v in vertices :
            plt.scatter(coord[v, 0], coord[v, 1])
            #plt.text(coord[v, 0]+0.05, coord[v, 1]+0.02, s = str(v))
        for e in edge_list :
            plt.plot( (coord[e[0]][0], coord[e[1]][0]), (coord[e[0]][1], coord[e[1]][1]), color = 'k', alpha = 0.4)

    return coord, vertices, edge_list, resistances, borders

# =================================================================================================
# ==============================  HEXAGONAL CONFIGURATION =========================================
# =================================================================================================
def give_radius_hexagonal(nlayers) :
    """
    give_radius(nlayers)
        
        Returns the radius for selection of vertices of the future graph in hexagonal graphs.
    
        Inputs
        ------
        nlayers : int
            Number of layers, min being 1, max is 10.
        Returns
        -------
        rad : float 
            Radius of selection.
    """
    rad = [1.5, 2., 3., 3.7, 4.7, 5.7, 6.6, 7.5, 8.3, 9.2]
    if nlayers < len(rad)+1 :
        return rad[nlayers-1]
    else :
        print 'Error : too large number of layers ! Maximum layers is ' + str(len(rad)+1)
        return;

def hexa_config(nlayers, noisy, lumen_pos_avg, lumen_pos_std, show, pbc, cen1 = np.array([20., 17.32050808])) :
    """
    hexa_config(nlayers, noisy, seed, lumen_pos_avg, lumen_pos_std, show, pbc, cen1 = np.array([20., 17.32050808]))
        Generate a hexagonal graph
    
        Inputs
        ------
        nlayers : int
            Number of layers of the graph. Minimum is one, max is set to be 10 layers.
    
        noisy : float
            Makes the output graph a noisy graph if True, a regular if False.
    
        lumen_pos_avg : float, optional
            If specified, correspond to the mean deviation of the noisy lattice with respect to the regular lattice
        lumen_pos_std : float, optional
            If specified, correspond to the standard deviation of the noisy lattice with respect to the regular lattice
    
        show : boolean, optional
            If specified True, shows the graph.
    
        pbc : boolean, optional, not implemented
            If specified, makes the graph periodic
        
        cen1 : float array, optional, default : np.array([20., 17.32050808])
            If specified, sets the center of the graph to the given position.
            
        Returns
        -------
        coord : array
            Array of coordinates of the graph.
    
        vertices : array
            Array of vertices of the graph.
    
        edge_list : array
            Array of vertices of the graph.
    
        resistances : array
            Array of resistance of the graph's edges.
    
        borders : set
            Set of vertices belonging to the broders of the graph.
    """
    if nlayers > 0 :    
        radius = give_radius_hexagonal(nlayers)
    else :
        print 'Error : you have less than 1 layer !'
        return;
    
    if pbc :
        print 'Periodic boundary conditions are not implemented yet !'
    N=40
    M=N*2
    a=1
    rad = 10
    center = np.array([20., 20.])
    mu, var = lumen_pos_avg, lumen_pos_std
    points = np.zeros((N, M, 2))
    
    # Generate a set of points on a regular lattice
    for i in range(N) :
        if i%2 == 0 :
            points[i] = np.array([0, i*a*np.sqrt(3)/2])*np.ones((M, 2)) + np.column_stack((a*np.arange(M), np.zeros(M)))
        else :
            points[i] = np.array([0.5*a, i*a*np.sqrt(3)/2])*np.ones((M, 2)) + np.column_stack((a*np.arange(M), np.zeros(M)))
    
    c = center*np.ones((2*N*N,2))
    cen1 = cen1
    radius1 = radius
    
    # Noisy hexagonal graph
    hexa = np.reshape(points, newshape=(N*M, 2))
    if noisy :
        hexa_noisy = hexa + np.random.normal(mu, var, size=(N*M, 2))
        vor_noisy2 = Voronoi(hexa_noisy)
        hexa_noisy_sel = np.column_stack((np.linalg.norm(hexa_noisy-c, axis=1) <= rad, np.linalg.norm(hexa_noisy-c, axis=1) <= rad)) * hexa_noisy
        hexa_noisy_sel = hexa_noisy_sel[np.nonzero(hexa_noisy_sel)[0]]
        vor_noisy = Voronoi(hexa_noisy_sel)
        vor_noisy_vert = np.column_stack((np.linalg.norm(vor_noisy.vertices-c[:len(vor_noisy.vertices)], axis=1) <= rad, np.linalg.norm(vor_noisy.vertices-c[:len(vor_noisy.vertices)], axis=1) <= rad)) * vor_noisy.vertices
        vor_noisy_vert = vor_noisy_vert[np.nonzero(vor_noisy_vert)[0]]
        
        vor, net = vor_noisy2, hexa_noisy
    # Regular hexagonal graph
    else :
        vor2 = Voronoi(hexa)
        hexa_sel = np.column_stack((np.linalg.norm(hexa-c, axis=1) <= rad, np.linalg.norm(hexa-c, axis=1) <= rad)) * hexa
        hexa_sel = hexa_sel[np.nonzero(hexa_sel)[0]]
        vor = Voronoi(hexa_sel)
        vor_vert = np.column_stack((np.linalg.norm(vor.vertices-c[:len(vor.vertices)], axis=1) <= rad, np.linalg.norm(vor.vertices-c[:len(vor.vertices)], axis=1) <= rad)) * vor.vertices
        vor_vert = vor_vert[np.nonzero(vor_vert)[0]]

        vor, net = vor2, hexa
    ############

    if show :
        f1, a1 = plt.subplots(figsize=(20, 20))
        v, s = plot_inside(vor, net, cen1, radius1, a1, show=show)
    else :
        v, s = plot_inside(vor, net, cen1, radius1)


    #### SORT ARRAY

    S = np.zeros((len(list(s)), 2), dtype=int)  
    V = np.zeros((len(list(v)), 2), dtype=int)

    for i in range(len(list(s))) :
        S[i] =  np.array([list(s)[i], i])
    
    for j in range(len(list(v))) :
        j1 = min(np.argwhere(S[:,0] == v[j,0])[0,0], np.argwhere(S[:,0] == v[j,1])[0,0] )
        j2 = max(np.argwhere(S[:,0] == v[j,0])[0,0] , np.argwhere(S[:,0] == v[j,1])[0,0] )
        V[j] = np.array([j1, j2 ])
    
    # list of the edges of the network
    edge_list = np.array(sorted_edge_list(V, S[:,1]))
    
    # borders of the graph/embryo
    borders = find_border(edge_list)
    
    # vertices and type of each vertex (TE, ICM, bicellular, ...)
    vertices = []
    for i in range(len(S[:,1])) :
        if i in borders :
            vertices += [[i, 1]]
        else :
            vertices += [[i, 0]]
    vertices = np.array(vertices)
    
    # coordinates of vertices
    coord = np.array([[vor.vertices[S[i,0]][0], vor.vertices[S[i,0]][1]] for i in range(len(S))])
    
    #resistances
    resistances = resistance_edge(coord=coord, edge_list=edge_list)
                
    return coord, vertices, edge_list, resistances, borders

# =================================================================================================
# =============================   BICELLULAR CONFIGURATION    =====================================
# =================================================================================================   
    
def gen_bicellular_graph(nbicellular, c, v, e, borders) :
    """
    gen_bicellular_graph(nbicellular, c, v, e, borders)
    
        Generate a graph with bicellular lumens
    
        Input
        -----
        nbicellular : int
            Number of bicellular lumens per edge. If 0, returns the graph as it is.
    
        Returns
        -------
        c_bicell : array
            Array of coordinates of each vertex
    
        v_bicell : array
            Array of vertices of the graph, with first column being the index of the vertex, 
            the second indicates the type of the vertex (0=ICM multi, 1=TE multi, 2=ICM bi, 3=TE bi)
            
            Example : v_bicell = np.array([[0, 1], [1, 3], [2, 0], [3, 0]])
            
        e_bicell : array
            Array of the edges of the graph.
            
        resistances : array
            Array of resistances associated to each connection between lumens.
    
        borders : set
            Set of vertices that belong to the borders of the graph.
    """
    borders = set(borders)
    if nbicellular == 0 :
        resistances = resistance_edge(c, e)
        return c, v, e, resistances, borders
    
    c_bicell, e_bicell, v_bicell = list(c), [], [[i, 0] for i in range(len(v))]
    for i in range(len(v_bicell)) :
        if i in borders :
            v_bicell[i][1] = 1
    
    v_start = len(v)-1
    for edge in e :
        nseg = nbicellular+1
        a = c[edge[0]]
        b = c[edge[1]]
        
        if edge[0] in borders and edge[1] in borders :
            # the new bicellular lumen is on the border
            flag = 3
        else :
            # the new bicellular lumen is in the inside
            flag = 2
            
        for i in range(nbicellular) :
            v_start += 1
            if i == 0 :
                e_bicell += [[edge[0], v_start]]
            else :
                e_bicell += [[v_start-1, v_start]]
            
            coord = a + (i+1) * (b-a) / (nseg)
            c_bicell += [coord]
            v_bicell += [[v_start, flag]]
            if flag == 3 :
                borders.add(v_start)
            
        e_bicell += [[v_start, edge[1]]]
        
    c_bicell = np.array(c_bicell)
    v_bicell = np.array(v_bicell)
    
    e_bicell = np.array(e_bicell)
    e_bicell = sorted_edge_list(edge_list=e_bicell, sorted_vertex=v_bicell[:, 0])
    
    resistances = resistance_edge(c_bicell, e_bicell)
    
    return c_bicell, v_bicell, e_bicell, resistances, borders

# =================================================================================================
# ======================================    MAIN    ===============================================
# =================================================================================================

def main(topology='hexagonal', outdir = 'cavitation/outputs/', subdir = 'config', nsim=1, nlayers=1, noisy=False, seed=None, vol_avg = 1., vol_std = 0.1, lumen_pos_avg = 0, lumen_pos_std = 0.05, show=False, pbc=False, args=[]) :
    
    global matplotlib_ON, src, configfile
    # default values
    chimeras = False
    gamma_border = 1.
    gamma_c_border = 1.
    pbc=False
    nbicellular = 0
    center = np.array([0., 0.])
    
    # Import arguments
    if args[0].endswith('.pat') :       # if pattern file is specified
        print 'Pattern file found'
        
        config = read_config(args[0])
        
        outdir, tpl, subdir, nsim = str(config['files']['outdir']), str(config['files']['tpl']), str(config['files']['subdir']), int(config['files']['nsim'])
        topology, nlayers, pbc, nbicellular, seed = str(config['network']['topology']), int(config['network']['nlayers']), eval(config['network']['pbc']), int(config['network']['nbicellular']), eval(config['network']['seed'])
        gamma_border, gamma_c_border = float(config['tensions']['gamma_border']), float(config['tensions']['gamma_c_border'])
        noisy, lumen_pos_avg, lumen_pos_std = eval(config['noisy']['noisy']), float(config['noisy']['lumen_pos_avg']), float(config['noisy']['lumen_pos_std'])
        vol_avg, vol_std = float(config['volume']['vol_avg']), float(config['volume']['vol_std'])
        chimeras, gamma_chimera = eval(config['chimera']['chimeras']), float(config['chimera']['gamma_chimera'])
        show = eval(config['display']['show'])
        
        if seed is not None :
            np.random.seed(int(seed))

        parameters = [[outdir, tpl, subdir, nsim], [topology, nlayers, pbc, nbicellular, seed], [gamma_border, gamma_c_border], [noisy, lumen_pos_avg, lumen_pos_std], [vol_avg, vol_std], [chimeras, gamma_chimera], [show]]

    elif len(args) > 0 :                # if use of command lines
        for arg in args :
            if arg.startswith('topology=') :
                topology = arg[len('topology='):]
                if not topology in ['hexagonal', 'triangular', 'chain'] :
                    print 'Topology not recognized. Stop generation of network.'
                    return;
            elif arg.startswith('outdir=') :
                outdir = arg[len('outdir='):]
                if not os.path.isdir(outdir) :
                    print 'Outdir does not exists... creation'
                    os.mkdir(outdir)

            elif arg.startswith('subdir=') :
                subdir = arg[len('subdir='):]
                
            elif arg.startswith('nsim=') :
                nsim = int(arg[len('nsim='):])
            
            elif arg.startswith('tpl=') :
                tpl = arg[len('tpl='):]
            
            elif arg.startswith('noisy=') :
                noisy = eval(arg[len('noisy='):])
            elif arg.startswith('nlayers=') :
                nlayers = int(arg[len('nlayers='):])
            elif arg.startswith('pbc=') :
                pbc = eval(arg[len('pbc='):])
            elif arg.startswith('seed=') :
                seed = int(arg[len('seed='):])
                np.random.seed(seed)                # initialize the seed
                
            elif arg.startswith('vol_avg=') :
                vol_avg = float(arg[len('vol_avg='):])
            elif arg.startswith('vol_std=') :
                vol_std = float(arg[len('vol_std='):])
            elif arg.startswith('lumen_pos_avg=') :
                lumen_pos_avg = float(arg[len('lumen_pos_avg='):])
            elif arg.startswith('lumen_pos_std=') :
                lumen_pos_std = float(arg[len('lumen_pos_std='):])
            elif arg.startswith('gamma_border=') :
                gamma_border = float(arg[len('gamma_border='):])
            elif arg.startswith('gamma_chimera=') :
                gamma_chimera = float(arg[len('gamma_chimera='):])
            elif arg.startswith('gamma_c_border=') :
                gamma_c_border = float(arg[len('gamma_c_border='):])
            elif arg.startswith('chimeras=') :
                chimeras = eval(arg[len('chimeras='):])
            elif arg.startswith('nbicellular=') :
                nbicellular = int(arg[len('nbicellular=')])
            elif arg.startswith('show=') :
                if matplotlib_ON :
                    show=eval(arg[len('show='):])
                else :
                    show=False
            elif arg.startswith('center=') :
                center = np.array([float(elem) for elem in arg[len('center='):].split(' ')])

    if topology == 'hexagonal' :
        center = np.array([20., 17.32050808]) # for hexagonal lattices
        coordinates, vertices, edges, resistances, borders = hexa_config(nlayers, noisy, lumen_pos_avg, lumen_pos_std, show, pbc, center)
        if nbicellular :
            coordinates, vertices, edges, resistances, borders = gen_bicellular_graph(nbicellular, coordinates, vertices, edges, borders)
    elif topology == 'triangular' :
        center = np.array([[0., 0.]])
        coordinates, vertices, edges, resistances, borders = triangular_config(nlayers, noisy, lumen_pos_avg, lumen_pos_std, show, pbc, center)
        if nbicellular :
            coordinates, vertices, edges, resistances, borders = gen_bicellular_graph(nbicellular, coordinates, vertices, edges, borders)
    elif topology == 'chain' :
        coordinates, vertices, edges, resistances, borders = chain_config(nlayers, noisy, lumen_pos_avg, lumen_pos_std, show, pbc, center)
        if nbicellular :
            print 'No bicelullar lumens for a chain...'
    
    if show :
        plt.axis('equal')
        plt.savefig(os.path.join(outdir, 'config.png'))
        plt.show()

    # Surface tensions
    gamma_list   = {elem : 1. for elem in range(len(vertices))}
    gamma_c_list = {elem : 1. for elem in range(len(vertices))}
    bridge_list  = {elem : 0 for elem in range(len(vertices))}
    
    
    # Chimeras
    if chimeras  and not topology == 'chain':
        print 'CHIMERAS'
        gamma_w_ch = 0.5 * (gamma_border + gamma_chimera)

        n = list(borders)[11]
        wild, ch, wild_ch = generate_wild_chimeras(n, edges, borders)
        
        for b in borders :
            if b in wild :
                gamma_list[b] = gamma_border
            elif b in ch :
                gamma_list[b] = gamma_chimera
            elif b in wild_ch :
                gamma_list[b] = gamma_w_ch
        l_types = lumen_type(vertices, borders, wild_list=wild, mutants_list=ch, wild_mutants_list=wild_ch, topology=topology)

    else :
        for b in borders :
            gamma_list[b] = gamma_border
            gamma_c_list[b] = gamma_c_border
        l_types = lumen_type(vertices, borders, topology=topology)
    
    # Generate multiple configurations
    if nsim > 1 :
        list_dir = gen_folders(nfold=nsim, name=subdir, abs_path=outdir)
        for n in range(nsim) :
            path = os.path.join(list_dir[n], os.path.split(tpl)[-1])
            shutil.copyfile(tpl, path)
            
            write_init(path, os.path.split(path)[0])
            
            path2 = os.path.join(list_dir[n], 'network')
            area_list = area_distribution(range(len(vertices)), vol_avg, vol_std)
            write_files(path2, range(len(vertices)), coordinates, edges, bridge_list, resistances, gamma_list, gamma_c_list, borders, area_list, l_types)
    else :
        area_list = area_distribution(range(len(vertices)), vol_avg, vol_std)
        dir_path = os.path.join(outdir, subdir)
        path = os.path.join(dir_path, 'network')
        try :
            os.mkdir(dir_path)
            os.mkdir(path)
        except :
            pass
        shutil.copyfile(tpl, os.path.join(dir_path, os.path.split(tpl)[1]))
        write_files(path, range(len(vertices)), coordinates, edges, bridge_list, resistances, gamma_list, gamma_c_list, borders, area_list, l_types)
    
    return;

if __name__ == "__main__" :

    if sys.argv[1] == 'help' :
        print(__doc__)
    else :
        if len(sys.argv) > 1 :
            args = sys.argv[1:]
            main(args=args)
        else :
            main()
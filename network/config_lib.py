"""
config_lib.py is a file containing useful functions for generating configurations for the graph of lumens

    Contents
    --------
    remove_single_edge
    find_border
X    plot_inside
    sorted_edge_list
    plot_edges
    plot_lumens
    resistance_edge
    get_neighbours
    generate_wild_chimeras
    lumen_type
    gen_folders
    area_distribution
    write_init
    write_lumen_coord
    write_lumen
    write_lumen_lumen
    write_bridge_lumen
    write_bridge_bridge
    write_files

    Libraries
    ---------
    numpy (np), networkx (nx), sys, os, configparser
    
    Imported Functions
    ------------------
    scipy.spatial : Voronoi, voronoi_plot_2d

    Notes
    -----
    Default used version is Python2.7

    Created by M. Le Verge--Serandour, September 2018
    Update : Feb, 2019
"""

import numpy as np
import networkx as nx
import sys
import os
from scipy.spatial import Voronoi, voronoi_plot_2d
import configparser
	
def remove_single_edge(edge_list, set_bridges) :
    """
    remove_single_edge(edge_list)
    
        Remove single edges from a list of edges
        
        Parameters
        ----------
        edge_list : array
            list of edge in the shape 
            np.array(['a', 'b'], ['b', 'c'], ['b', 'd'], ['d', 'c'])
            
            where ['a', 'b'] is a single edge (one occurence of vertex 'a')
            
        set_bridges : set
            set of bridges.
        
        Returns
        -------
        
        
    """
    s = set(set_bridges) # copy for safety
    e = np.array(edge_list) # copy for safety
    
    # make an array of the occurence of each vertex of the graph
    occurence_list = np.unique(np.reshape(e, 2*len(e)), return_counts=True)
    
    # select all the vertex appearing only once ("single vertex")
    single_edge_list = occurence_list[0][np.argwhere(occurence_list[:][1] == 1)]
    
    index_list = []
    
    for elem in single_edge_list :
        # remove "single vertex"
        s.remove(elem[0])
        
        # the "single vertex" can be either i or j in edge (i, j), so check where it is
        w1, w2 = np.argwhere(e[:,0] == elem), np.argwhere(e[:,1] == elem)

        if len(w1) != 0 :
            index = np.argwhere(e[:,0] == elem[0])[0,0]
            #print 'index = ' + str(index)
        
        else :
            index = np.argwhere(e[:,1] == elem[0])[0,0]
            #print 'index = ' + str(index)
        
        index_list += [index]
    
    return np.delete(e, index_list, 0), s

def find_border(edge_list) :
    """
    find_border(edge_list)
    
        Find the borders of a hexagonal graph
    
        Input
        -----
        edge_list : array
            List of edges of the graph
    
        Returns
        -------
        border_set : set
            Set of vertices of the border
    """
    G = nx.Graph([(edge_list[i,0], edge_list[i,1]) for i in range(len(edge_list))])
    occurence_list = np.unique(np.reshape(edge_list, 2*len(edge_list)), return_counts=True)
    
    # list of vertex of degree 2
    sec_edge_list = occurence_list[0][np.argwhere(occurence_list[:][1] == 2)]
    # list of vertex of degree 3
    three_edge_list = occurence_list[0][np.argwhere(occurence_list[:][1] == 3)]


    sec = np.reshape(sec_edge_list, newshape=(len(sec_edge_list)))
    border_set = set(sec)
    inner_set = set()

    for elem in three_edge_list :
        for neigh in G[elem[0]].keys() :
            if len(G[neigh]) == 2 :
                border_set.add(elem[0])
    return border_set

def plot_inside(vr, net, center, radius, axes=None, show=False) :
    """
    plot_inside
    
        Returns the graph
    """
    list_edges = []
    set_bridges = set()
    for vpair in vr.ridge_vertices :
        
        if vpair[0] >= 0 and vpair[1] >= 0 : # ie : if the edge is full            
                
            v0 = vr.vertices[vpair[0]]
            v1 = vr.vertices[vpair[1]]
            
            # if the edge is inside the defined circle (and not only the bridge...)
            if np.sum((v0-center)*(v0-center)) <= radius*radius and np.sum((v1-center)*(v1-center)) <= radius*radius :
                list_edges += [vpair]
                set_bridges.add(vpair[0])
                set_bridges.add(vpair[1])
                #axes.plot([v0[0], v1[0]], [v0[1], v1[1]], color='blue', linewidth=4, alpha = 0.2)
    
    list_edges, set_bridges = remove_single_edge(list_edges, set_bridges)
    
    if show :
        for vpair in list_edges :
            v0 = vr.vertices[vpair[0]]
            v1 = vr.vertices[vpair[1]]
            axes.plot([v0[0], v1[0]], [v0[1], v1[1]], color='black', linewidth=4, alpha = 0.4)
    
    return list_edges, set_bridges

def sorted_edge_list(edge_list, sorted_vertex) :
    """
    sorted_edge_list(edge_list, sorted_vertex)
    
        Sort the edge list : if [[i1, i2], [j1, j2]] then i1 < i2 AND i1 < j1
    
        Inputs
        ------
        edge_list : array
            Array of edges of the graph
        sorted_vertex : array
            Sorted array of the vertices of the graph
        Returns
        -------
        sorted_list : array
            Sorted array of the edges of the graph.
    """
    G = nx.Graph([(edge_list[i, 0], edge_list[i,1]) for i in range(len(edge_list))])
    sorted_list = []
    for elem in sorted_vertex :
        #print elem, np.sort(G[elem].keys())
        for neigh in np.sort(G[elem].keys()) :
            #print elem, neigh
            if elem < neigh : 
                sorted_list += [[elem, neigh]]
    return sorted_list

def plot_edges(coord, edge_list) :
    """
    plot_edges(coord, edge_list)
    
        Plot the edges of the graph
    
        Inputs
        ------
        coord : array
            Coordinates of the vertices of the graph
        edge_list : array
            Array of the edges of the graph
    """
    
    
    for i in range(len(edge_list)) :
        plt.plot((coord[edge_list[i][0]][0], coord[edge_list[i][1]][0]), (coord[edge_list[i][0]][1], coord[edge_list[i][1]][1]), color = 'black', alpha=0.5)
          
def plot_lumens(step, coord, edge_list, area, save = False, show = False, savename='pic.png', folder = 'pic/') :
    """
    plot_lumens(step, coord, edge_list, area, save = False, show = False, savename='pic.png', folder = 'pic/')
    
        Plot the graph of lumens at a given step.
    
        Inputs
        ------
        step : int
            Time step when to plot the graph
        coord : array
            Coordinates of the vertices of the graph
        edge_list : array
            Array of the edges of the graph
        area : array
            Array of the areas at step
        save : boolean, optional, default : False
            True if saving the picture, False otherwise
        show : boolean, optional, default : False
            True if showing the picture, False otherwise
        savename : str, optional, default : pic.png
            Name of the picture if saved
        folder : str, optional, default : pic
            Name of the folder where the picture is saved.
    
    """
    
    plt.figure(figsize=(8, 8))
    
    plot_edges(coord, edge_list)
        
    for i in range(len(coord)) :
        plt.scatter(coord[i, 0], coord[i, 1], s = 500*area[i])
        
    plt.axis('equal')
    if save :
        plt.savefig(folder + savename)
    if show :
        plt.show()
    plt.close()

    
def resistance_edge(coord, edge_list) :
    """
    resistance_edge(coord, edge_list)
    
        Calculate the resistances of each edge
        
        Inputs
        ------
        coord : array
            Coordinates of the vertices of the graph
        edge_list : array
            Array of the edges of the graph
    
        Returns
        -------
        R : array
            Resistances of the edges
    """
    R = np.zeros(len(edge_list))
    for i in range(len(edge_list)) :
        R[i] = np.linalg.norm(coord[edge_list[i][0]] - coord[edge_list[i][1]])
    return R

def get_neighbours(n, edge_list) :
    """
    get_neighbours(n, edge_list)
    
        Get the neighbours of a vertex n
    
        Inputs
        ------
        n : int
            Index if the vertex
        edge_list : array
            Array of the edges of the graph
    
        Returns
        -------
        neigh : list
            List of the neighbours of n
    """

    neigh = []
    for i in edge_list :
        if n in i :
            if i[0] == n :
                neigh += [i[1]]
            else :
                neigh += [i[0]]
    return neigh

def generate_wild_chimeras(n, edge_list, border_set) :
    """
    generate_wild_chimeras(n, edge_list, border_set)
    
        Given a graph, generate the list of WT, mutants and WT-mutants lumens
        Inputs
        ------
        n : int
            Initial lumen
        Returns
        -------
        wild : list
            List of WT lumens
        mutants : list
            List of mutant lumens
        wild_mutants : list
            List of WT-mutant lumens
    """
    if n not in border_set :
        print 'The initial chimera lumen is not on the border !'
        return;
    wild = list(border_set)
    
    wild_mutants = []
    mutants = []

    mutants += [n]
    wild.remove(n)

    while len(mutants) < 8 :
        next_ch = mutants[-1]
        for neighbor in get_neighbours(next_ch, edge_list) :
            if neighbor in border_set and neighbor not in mutants :
                mutants += [neighbor]
                wild.remove(neighbor)
                

    # find the lumen at interface between wild and mutants
    for elem in mutants :
        for neighbor in get_neighbours(elem, edge_list) :
            if neighbor not in mutants and neighbor in border_set :
                wild_mutants += [neighbor]
                wild.remove(neighbor)
    
    return wild, mutants, wild_mutants

def lumen_type(lumen_list, border_set, wild_list=[], mutants_list = [], wild_mutants_list=[], topology='hexagonal') :
    """
    lumen_type(lumen_list, border_set, wild_list=[], mutants_list = [], wild_mutants_list=[])
    
        Returns the type of each lumen
    
        Inputs
        ------
        lumen_list : list
            List of the vertices of the graph
        border_set : set
            Set of vertices of the borders
    
        wild_list : list, optional, default : []
            
        mutants_list : list, optional, default : []
    
        wild_mutants_list : list, optional, default : []
            
        Returns
        -------
        l : list
            List of types of lumens
    """
    l = {}
    if topology == 'chain' :
        for n in range(len(lumen_list)) :
            l[n] = 'ICMmi'
        return l
    # if no mutants :
    if len(mutants_list) == 0 :
        for elem in range(len(lumen_list)) :#[:, 0] :
            if elem in border_set :
                s = 'TE'
                if lumen_list[elem][1] == 1 :
                    s += 'mi'
                elif lumen_list[elem][1] == 3 :
                    s += 'bi'
                l[elem] = s
            else :
                s = 'ICM'
                if lumen_list[elem][1] == 0 :
                    s += 'mi'
                elif lumen_list[elem][1] == 2 :
                    s += 'bi'
                l[elem] = s
    
    # if mutants :
    else :
        for elem in lumen_list[:, 0] :
            if elem in border_set :
                if elem in wild_list :
                    l[elem] = 'TEmi'
                elif elem in mutants_list :
                    l[elem] = 'mutantsmi'
                else :
                    l[elem] = 'wild_mutantsmi'
            else :
                l[elem] = 'ICMmi'
    return l

def gen_folders(nfold, name='config', abs_path = '~/cavitation/network/') :
    """
    gen_folders(nfold, name='config', abs_path = '~/cavitation/outputs/')
    
        Generate folders for simulations
        
        Inputs
        ------
        nfold : int
            Number of folders to generate
        name : str, optional, default : config
            Name of the folders
        abs_path : str, path, default : ~/cavitation/network/
            Absolute path for the folders.

        Returns
        -------
        list_dir : list
            List of the directories created
    """
    list_dir = []
    for n in range(nfold) :
        s = os.path.join(abs_path, name+ str(n).zfill(4))
        list_dir += [s]
        if not os.path.isdir(s) :
            os.mkdir(s)
            os.mkdir(os.path.join(s, 'network'))
    return list_dir
        
def area_distribution(vertices, vol_avg, vol_std, threshold = 0.1) :
    """
    volume_distribution(conversion_list, vol_avg, vol_std, threshold = 0.1)
    
    Parameters
    ----------
    vertices : list
        List of the vertex of the network
    vol_avg : float
        Average volume
    vol_std : float
        Average volume
    threshold : float, optional, default : 0.1
        Minimal area allowed for the distribution
    
    Returns
    -------
    area_list : dict
        Dictionnary of areas, indexed by corresponding vertex index.
    """
    # Area distribution
    threshold = 0.1

    area_list = {}
    for elem in vertices :
        p = np.random.normal(loc=vol_avg, scale=vol_std)
        while p <= threshold :
            p = np.random.normal(loc=vol_avg, scale=vol_std)
        area_list[elem] = p
    return area_list
    
def write_init(init_config, pathname) :
    """
    write_init(init_config, pathname)
    
        Write into the init_config file the path of the folder where to run the simulation
    
        Inputs
        ------
        init_config : str, path
            Init configuration file, usually filename.ini
    
        pathname : str, path
            Path of the folder where to run the simulation.
    
    
    """
    config = configparser.ConfigParser()
    config.read(init_config)

    config.set('network', 'path', pathname)

    with open(init_config, 'wb') as configfile :
        config.write(configfile)

def write_lumen_coord(folder, lumen_list, lumen_positions) :
    """
    write_lumen_coord(folder, lumen_list, lumen_positions)
    
        Write the coordinates of the lumens in file 'lumen_coord.dat'. The line index in the file is the vertex index.
    
        Inputs
        ------
        folder : str, path
            Folder where to store the file.
        lumen_list : array or list
            Array of vertices of the graph, namely the lumens.
        lumen_positions : array or list
            Array of the coordinates of the lumens.
        
        Returns
        -------
        0 if there is no problem.  
    
    
        Structure
        ---------
        pos_x   pos_y
        0.      0.
        1.4     3.5
        ...      ...
    
    """
    filename = 'lumen_coord.dat'

    pos_x, pos_y = lumen_positions[lumen_list[:]][:,0], lumen_positions[lumen_list[:]][:,1]
    
    fi = open(os.path.join(folder, filename), 'w')
    
    fi.write('#  # coordinates of each lumen, row represents the ID of the lumen \n')
    fi.write('#  # coordinates in this model do not change over time, empty lumen are included\n')
    fi.write('#  # x   y\n')
    
    for i in range(len(pos_x)) :
        fi.write(str(pos_x[i]) + '\t' + str(pos_y[i]) + '\n' )
        
    fi.close()
    return 0

def write_lumen(folder, lumen_list, gamma_list, gamma_c_list, border_list, area_list, lum_type) :
    """
    write_lumen(folder, lumen_list, gamma_list, gamma_c_list, border_list, area_list, lum_type)
        
        Write the lumen.dat file
    
        Parameters
        ----------
        folder : str, path
            Folder where to store the file.
        lumen_list : array or list
            Array of vertices of the graph, namely the lumens.
        gamma_list : list
            List of the tension gamma associated with each vertex
        gamma_c_list : list
            List of the (adhesion) tension gamma_c associated with each vertex
        border_list : list
            List of the vertex belonging to the borders of the graph.
        area_list : list    
            List of the initial area of each vertex.
        lum_type : list
            List of the type of each vertex (0 : ICM-multi, 1 : TE-multi, 2 : ICM-bi, 3 : TE-multi, ...)
    
        Returns
        -------
        1 if no problem.
    
        Structure
        ---------
        gamma1   gamma2     gamma_c    area     boundary   
        1         1.        0.4       3          0
        1         1.        0.4       2.3        1
    """
    def flag(type_lum) :
        if type_lum == 'ICMmi' :
            return 0
        elif type_lum == 'TEmi' :
            return 1
        elif type_lum == 'ICMbi' :
            return 2
        elif type_lum == 'TEbi' :
            return 3
        elif type_lum == 'mutantsmi' :
            return 4
        elif type_lum == 'wild_mutantsmi' :
            return 5
        elif type_lum == 'mutantsbi' :
            return 6
        
    filename = 'lumen.dat'
    fi = open(os.path.join(folder, filename), 'w')
    
    e = np.sort(list(lumen_list))
    
    for i in range(len(e)) :        
        fi.write(str(gamma_list[e[i]]) + '\t' + str(gamma_list[e[i]]) + '\t' + str(gamma_c_list[e[i]]) + '\t' + str(area_list[e[i]]) + '\t' + str(flag(lum_type[i])) + '\n')

    fi.close()
    
    return 1
    
def write_lumen_lumen(folder, edge_list, R_list) :
    """
    write_lumen_lumen(folder, edge_list, R_list)
        
        Write the lumen_lumen.dat file
    
        Parameters
        ----------
        folder : str, path
            Folder where to store the file.
        edge_list : array or list
            Array of the edges of the graph.
        R_list : list
            List of the hydraulic resistance of each edge (no friction for the moment)
      
        Returns
        -------
        2 if no problem.
        NB : the file is such that for edge = [ID1, ID2], ID1 < ID2. Moreover, from one line to the next,
                the ID1 is in growing order.
    
        Structure
        ---------
        ID1     ID2     distance
        0       1       1.0
        0       4       4.
        1       2       1.0
        1       4       3.
        2       3       4.5
        3       4       3.1
    """
    
    filename = 'lumen_lumen.dat'
    
    fi = open(os.path.join(folder, filename), 'w')
    fi.write('# # triangle with lumen in the center\n')
    fi.write('#  # containes the connections between lumen, in the structure of a graph \n')
    fi.write('#  # lumen have an ID starting from 0, sorted such that, smallest on top and ID1 < ID2 \n')
    fi.write('#  # ID1 ID2 distance\n')
    for i in range(len(edge_list)) :
        fi.write(str(edge_list[i][0]) + '\t' + str(edge_list[i][1]) + '\t' + str(R_list[i]) + '\n' )
        
    fi.close()
    return 2

def write_bridge_lumen(folder, edge_list, bridge_list, R_list) :
    """
    write_bridge_lumen(folder, edge_list, bridge_list, R_list)
        
        Write the bridge_lumen.dat file
    
        Parameters
        ----------
        folder : str, path
            Folder where to store the file.
        edge_list : array or list
            Array of the edges of the graph.
        bridge_list : array or list
            Array of the bridges of the graph.
        R_list : list
            List of the hydraulic resistance of each edge (no friction for the moment)
      
        Returns
        -------
        3 if no problem.
    
    
        NB : the file is such that for edge = [ID1, ID2], ID1 < ID2. Moreover, from one line to the next,
                the ID1 is in growing order.
        NB : the ID1 is ALWAYS a bridge
        NB : the indices of the bridges start from 0 when created. 
    
        Structure
        ---------
        ID1     ID2     distance
        0       1       1.0
        0       4       4.
        1       2       1.0
        1       4       3.
        2       3       4.5
        3       4       3.1
    """
    filename = 'bridge_lumen.dat'
    
    fi = open(os.path.join(folder, filename), 'w')
    fi.write('# # triangle with lumen in the center\n')
    fi.write('#  # containes the connections between lumen and bridges (empty lumen), in the structure of a graph\n')
    fi.write('#  # lumen and bridges have an ID starting from 0; there are lumen and bridges with the same ID, sorted by the ID of the bridge\n')
    fi.write('#  # bridge lumen distance\n')

    fi.close()
    return 3

def write_bridge_bridge(folder, bridge_list, R_list) :
    """
    write_bridge_bridge(folder, bridge_list, R_list)
        
        Write the bridge_bridge.dat file
    
        Parameters
        ----------
        folder : str, path
            Folder where to store the file.
        bridge_list : array or list
            Array of the bridges of the graph.
        bridge_list : array or list
            Array of the bridges of the graph.
        R_list : list
            List of the hydraulic resistance of each edge (no friction for the moment)
      
        Returns
        -------
        4 if no problem.
    
        NB : the file is such that for edge = [ID1, ID2], ID1 < ID2. Moreover, from one line to the next,
                the ID1 is in growing order.
        NB : the indices of the bridges start from 0 when created. 
                See the correspondance with preexisting vertices in conversion.dat
    
        Structure
        ---------
        ID1     ID2     distance
        0       1       1.0
        0       4       4.
        1       2       1.0
        1       4       3.
        2       3       4.5
        3       4       3.1
    """
    filename = 'bridge_bridge.dat'
    
    fi = open(os.path.join(folder, filename), 'w')
    fi.write('# # triangle with lumen in the center  \n')
    fi.write('#  # containes the connections between bridges, in the structure of a graph \n')
    fi.write('#  # bridges have an ID starting from 0, sorted such that, smallest on top and ID1 < ID2 \n')
    fi.write('#  # ID1 ID2 distance\n')
    fi.close()
    return 4

def write_files(folder, vertices, coordinates, edges, bridge_list, resistances, gamma_list, gamma_c_list, borders, area_list, lum_type) :
    """
    write_files(folder, vertices, coordinates, edges, bridge_list, resistances, gamma_list, gamma_c_list, borders, area_list, lum_type
        
        Write the files for the simulation.
    
        Inputs
        ------
        folder : str, path
            Folder where to store the file.
        vertices : array or list
            Array of vertices of the graph, namely the lumens.
        coordinates : array or list
            Array of the coordinates of the lumens.
        edges : array or list
            Array of the edges of the graph.
        bridge_list : array or list
            Array of the bridges of the graph.
        resistances : list
            List of the hydraulic resistance of each edge (no friction for the moment)
        gamma_list : list
            List of the tension gamma associated with each vertex
        gamma_c_list : list
            List of the (adhesion) tension gamma_c associated with each vertex
        borders : list
            List of the vertex belonging to the borders of the graph.
        area_list : list    
            List of the initial area of each vertex.
        lum_type : list
            List of the type of each vertex (0 : ICM-multi, 1 : TE-multi, 2 : ICM-bi, 3 : TE-multi, ...)
    
        Returns
        -------
        0 if no problem
    """
    wlc = write_lumen_coord(folder, vertices, coordinates)
    wl = write_lumen(folder, vertices, gamma_list, gamma_c_list, borders, area_list, lum_type)
    wll = write_lumen_lumen(folder, edges, resistances)
    wlb = write_bridge_lumen(folder, edges, bridge_list, resistances)
    wbb = write_bridge_bridge(folder, bridge_list, resistances)
    return 0
    

# Library for the dynamics of a lumen network
# The lumen are 2 dimensional and symmetric and connected with 1 dimensional tubes
#
# Created by A. Mielke, 2018
# Modified by M. Le Verge--Serandour on 8/04/2019

"""
    network.py conf.init
    
    Defines the class network and associated functions
    
    Imports
    -------
    Libraries : numpy, os, math
    
    Created by A. Mielke
    Modified by H. Turlier on 8/06/2018
    Modified by M. Le Verge--Serandour on 8/04/2019
"""

import numpy as np
import math
import os

class network:

    def __init__(self, network_folder, out_path, t_step, tube_radius = 0.01, friction = 1, swelling = False, swelling_rate=0., save_area_dat=False):
        """
        Initialization of the object network
        All properties needed for the simulation are read and initialized
            
            Input
            -----
            network_folder : str
        
            out_path : str, path-like
                
            t_step : float
                Time step of the simulation. Note that if the simulation is adaptative, this time step will change.
            tube_radius : float, optional, default = 0.01
                Radius of the tube connecting lumens. Define the condition for empty lumens.
            friction : float, optional, default = 1
                Friction constant for the fluid circulating through pipes.
            swelling : bool, optional, default = False
                Swelling option for the simulation. True if swelling is included, False otherwise.
            swelling_rate : float, optional, default = 0.
                Swelling rate value in case the swelling is considered. Make sure the rate is not to big to avoid non-converging simulations.
            save_area_dat : bool, optional, default = False
                Save area option. True if areas are saved in area.dat, False otherwise.
            """
        self.network_folder = network_folder
        
        # Reading  properties of the lumen
        self.gamma_lumen, self.gamma_contact, self.area = np.loadtxt(os.path.join(network_folder, 'lumen.dat'), dtype = float, usecols = [0,2,3], unpack = True)
        

        # Reading links between two lumen
        self.lumen_lumen = self.read_lumen_lumen(os.path.join(network_folder, 'lumen_lumen.dat'))
        
        # Reading links between bridge and lumen
        self.bridge_lumen, self.num_bridges = self.read_bridge_lumen(os.path.join(network_folder, 'bridge_lumen.dat'))

        # Reading links between two bridges
        self.bridge_bridge, self.num_bridges = self.read_bridge_bridge(os.path.join(network_folder, 'bridge_bridge.dat'), self.num_bridges)

        # Surface tension ratio
        self.alpha = self.gamma_contact/(2*self.gamma_lumen)
        self.delta = np.full(len(self.alpha), 1) # Possibility of asymmetric lumen is not included
        
        # Resistances
        self.tube_radius = tube_radius # Radius of the tube connecting the lumen and the bridges
        self.friction = friction # Friction coefficient; friction * length = resistance
        
        # Opening angle of the lumen (angle between curvature and tube)
        self.theta = self.set_theta()
        # Area factor for expressing the pressure in terms of the area instead of the radius
        self.area_factor = self.set_area_factor()


        # Ending time: time at which only one lumen is remaining
        self.end_time = 0
        # Time step for the output of the area evolution
        self.time_step = t_step

        # Creating output file for the area evolution, events, error messages
        self.save_area(start = True, out_path = out_path)
        self.save_event('', start = True, out_path = out_path)
        self.save_error('', start = True, out_path = out_path)
        
        # Area distribution after only one lumen is remaining
        self.final_area = []
        
        # Current time step of the simulation
        self.current_time = 0
    
        # List of empty lumen (area < tube_radius **2)
        self.empty_list = np.zeros(len(self.alpha))

        # Swelling
        self.swelling_bool = swelling
        self.swelling_rate = swelling_rate
        
        # Save area
        self.save_area_dat = save_area_dat

############################################################################################################################
########################################################## Dynamics ########################################################
############################################################################################################################

    def flux(self, t, state):
        """
        Determines the flux/ area change for each lumen of the network, main function of network.py
            
            Input
            -----
            self : network object
                Needs to be called by a class object
            t : float
                Actual time step (not needed for the calculation of the flux, but required for the used integration method in network_simulation.py
            state : float array
                The current area of the lumens
            
            Returns
            -------
            flux : float array
                Contains the area change for each lumen in dt
        
                """
        # Initialization of the array containing the area change (index == lumen ID)
        flux = []
        self.current_time = t
        for i in range(len(self.alpha)):
            flux.append(0)
        
        # If only one lumen remains -> End of simulation, flux is zero (needed as for the integration method used, no dynamic stop is possible)
        if(np.sum(self.empty_list) >= len(self.alpha) - 1):
            if(self.end_time == 0):
                # Setting the end time for the output file area.log
                self.end_time = t
        
        # more than one lumen remaining: calculation of the flux
        else:
            
            # Adapting network to new state: Empty lumen are removed and graph is reconnected
            self.area = state
            self.remove_empty_lumen()
            
            # Area change between directly connected lumen
            flux = self.flux_lumen(flux)
            
            # Calculating artificial pressure at each bridge; linear system of equations, with flux(bridge) = 0, the bridge does not gain or loose area
            pressure_bridges = self.pressure_bridges()
            
            # Area change between lumen-bridges
            flux = self.flux_bridges(flux, pressure_bridges)
        
            # Area change due to swelling
            if self.swelling_bool :
                flux = self.flux_swelling(flux)
            
        # Saving area for the time step given in the configuration file
        if self.save_area_dat :
            self.save_area()
            
        self.t_old = t
        
        
        if(np.abs(np.sum(flux)) > self.tube_radius ** 2):
            error = 'total flux is non-zero: total flux = %f' % (np.sum(flux))
            self.save_error(error)

        return flux

    def flux_lumen(self,flux):
        """
        Determines the flux/ area change for each lumen due to the connection between lumen and lumen
            
            Input
            -----
            self    network object
                    needs to be called by a class object
            flux    float array
                    vector containing the area change for each lumen; index = lumen ID
            
            Returns
            -------
            flux    float array
                    area changes due to lumen-lumen connection added to the vector passed
                """
        
        # for each connection between two lumen
        for line in range(len(self.lumen_lumen)):
            lumen_1 = int (self.lumen_lumen[line][0]) # first lumen
            lumen_2 = int (self.lumen_lumen[line][1]) # second lumen
            # flux from lumen 2 to lumen 1
            fl = (self.pressure(lumen_2) - self.pressure(lumen_1))*self.friction/self.lumen_lumen[line][2]
            flux[lumen_1] += fl
            flux[lumen_2] -= fl
        return flux
    
    def pressure_bridges(self):
        """
        Determines the pressure at each bridge
        for each bridge the total flux is 0, meaning that the bridge does not gain or loose area
        this gives a linear equation system, which can be solved
        The connections are taken from the files bridge_lumen.dat and bridge_bridge.dat
        For Information about the equations see the documentation to the code
            
            Input
            -----
            self : network object
                Needs to be called by a class object
            
            Returns
            -------
            pressure_bridges : float array
                Pressure at each bridge
            """
            
        R_sum = np.zeros(self.num_bridges, dtype = float) # sum of the resistences around one bridge
        P_over_R_sum = np.zeros(self.num_bridges, dtype = float) # sum of pressure over resistance between one bridge and all directly connected lumen
        matrix_bridges = np.zeros([self.num_bridges, self.num_bridges], dtype= float) # matrix to calculate the pressure at each bridge
        
        # For each connection between bridge and lumen
        for line in self.bridge_lumen:
            bridge = int(line[0])
            lumen = int(line[1])
            R_sum[bridge] += 1./line[2]*self.friction
            P_over_R_sum[bridge] += self.pressure(lumen)/line[2]*self.friction
        
        # For each connection between bridge and bridge
        for line in self.bridge_bridge:
            bridge1 = int(line[0])
            bridge2 = int(line[1])
            matrix_bridges[bridge1][bridge2] = 1./line[2]*self.friction
            matrix_bridges[bridge2][bridge1] = 1./line[2]*self.friction
            R_sum[bridge1] += 1./line[2]*self.friction
            R_sum[bridge2] += 1./line[2]*self.friction
        for line in range(self.num_bridges):
            matrix_bridges[line][line] = -R_sum[line]
        
        # Solving linear problem with the pressure at each bridge as solution
        pressure_bridges = np.linalg.solve(matrix_bridges, -P_over_R_sum)
        
        return pressure_bridges;

    def flux_bridges(self, flux, pressure_bridges):
        """
        Determines the flux/ area change for each lumen due to the connection between lumen and bridge
            
            Input
            -----
            self : network object
                Needs to be called by a class object
            
            Returns
            -------
            flux : float array
                Area changes due to bridge-lumen connection added to the vector passed
            """
        
        # Area change in one bridge; should be 0; calculated as control value
        flux_bridge = np.zeros(self.num_bridges, dtype = float)
        
        # For each connection between bridge and bridge
        for line in self.bridge_bridge:
            bridge1 = int(line[0])
            bridge2 = int(line[1])
            fb = (pressure_bridges[bridge2] - pressure_bridges[bridge1])*self.friction/line[2]
            flux_bridge[bridge1] += fb
            flux_bridge[bridge2] -= fb
            
        # For each connection between bridge and lumen
        for line in self.bridge_lumen:
            bridge = int(line[0])
            lumen = int(line[1])
            fl = (pressure_bridges[bridge] - self.pressure(lumen))*self.friction/line[2]
            flux[lumen] += fl
            flux_bridge[bridge] -= fl

        for i in range(len(flux_bridge)):
            if (np.abs(flux_bridge[i]) > self.tube_radius ** 2):
                error = 'total flux of bridge %d is non-zero: total flux = %f' % (i,flux_bridge[i])
                self.save_error(error)
        return flux


    def flux_swelling(self, flux) :
        """
        Determines the flux/ area change for each lumen due to sewlling

            Input
            -----
            self : network object
                Needs to be called by a class object
        
            Returns
            -------
            flux : float array
                Area changes due to bridge-lumen connection added to the vector passed 
        """
        
        # for each lumen     (lumen is the index of the lumen's area)
        for lumen in range(len(self.area)) :
            # if not empty
            if not self.area[lumen] < 2*self.tube_radius ** 2 :
                # then add the swelling contribution
                flux[lumen] += self.swelling(lumen)
            
        return flux

############################################################################################################################
###################################################### Removing Functions #####################################################
############################################################################################################################

    def remove_empty_lumen(self):
        """
        Determines and removes empty lumen
        Calls a function to obtain a list of empty lumen and passes the list to a function to remove them and reconnect the network
            
            Input
            -----
            self : network object
                Needs to be called by a class object
            
            Returns
            -------
            no return
            """
            
        empty_lumen_list = []
        # Creating a list of empty lumen
        empty_lumen_list = self.get_empty_lumen()
        # Removing empty lumen and reconnecting the network
        if (len(empty_lumen_list) > 0 ):
            event = 'empty lumen: ' + ' '.join(map(str, empty_lumen_list))
            #print event
            self.save_event(event)
            self.remove_lumen(empty_lumen_list)
        return;
    

    def remove_lumen(self, lumen_to_remove):
        """
        Removes the lumen that are passed and connects the neighbors of these lumen
            
            Input
            -----
            self : network object
                Needs to be called by a class object
            lumen_to_remove : int list
                List of lumen to be removed
            
            Returns
            -------
            no return
            """
        
        # For each lumen that has to be removed
        for lumen in lumen_to_remove:
            neighbours = self.get_neighbours(lumen) # List of connected lumen
            bridges = self.get_bridges(lumen)       # List of connected bridges
            self.save_event('lumen ' + str(lumen) + ' neighbours ' + str(neighbours))
            self.save_event('lumen ' + str(lumen) + ' bridges ' + str(bridges))
            # Lumen had two connections, this means that it disappears and the two connected parts get directly connected, the resistance for the new link is the sum of the resistance of the two previous connections
            test=True
            if(len(neighbours) + len(bridges) == 2):
                # Lumen was connected to two lumen -> new connection between lumen and lumen
                if(len(neighbours) == 2):
                    self.create_link([neighbours[0][0], neighbours[1][0], neighbours[0][1] + neighbours[1][1]])
                    #print 'lumen_lumen connexion (' + str(neighbours[0][0]) + ', ' + str(neighbours[1][0]) + ')'
                # Lumen was connected to a lumen and a bridge -> new connection between lumen and bridge
                if(len(neighbours) == 1 and len(bridges)==1):
                    self.create_bridge_lumen([bridges[0][0], neighbours[0][0], bridges[0][1] + neighbours[0][1]])
                    #print 'lumen_bridge connexion (' + str(bridges[0][0]) + ', ' + str(neighbours[0][0]) + ')'
                # Lumen was connected to two bridges -> new connection between bridge and bridge
                if(len(bridges)==2):
                    self.create_bridge_bridge([bridges[0][0], bridges[1][0], bridges[0][1] + bridges[1][1]])
                    #print 'bridge_bridge connexion (' + str(bridges[0][0]) + ', ' + str(bridges[1][0]) + ')'
                
                self.create_bridge(neighbours, bridges, lumid=lumen)
                
            # Lumen had more than two connections -> becomes a bridge, the resistances remain the same but the connections are changed to connections to a bridge
            if(len(neighbours) + len(bridges) > 2):
                self.create_bridge(neighbours, bridges, lumid=lumen)
        return;


    def remove_link(self, lumen_1, lumen_2):
        """
        Removes a connection between two lumen
            
            Input
            -----
            self : network object
                Needs to be called by a class object
            lumen_1 : int
                First lumen of the connection
            lumen_2 : 
                Second lumen of the connection
            
            Returns
            -------
            no return
            """
        # Due to data structure first lumen must be smaller than second lumen
        if(lumen_1 > lumen_2):
            n = lumen_1
            lumen_1 = lumen_2
            lumen_2 = n

        # Find connection in lumen_lumen file and remove it
        line = 0
        # For each line in lumen_lumen until connection is found
        while (line < len(self.lumen_lumen)):
            # If connection is found removing it
            if(self.lumen_lumen[line][0] == lumen_1 and self.lumen_lumen[line][1] == lumen_2):
                
                event = 'link lumen %d to lumen %d removed' % (lumen_1, lumen_2)
                #print event
                self.save_event(event)
                
                link = [lumen_1, lumen_2, self.lumen_lumen[line][2]]
                self.lumen_lumen.remove(link)

                break;
            # Look at next line
            else: line += 1



############################################################################################################################
###################################################### Get Functions #####################################################
############################################################################################################################

    def get_empty_lumen(self):
        """
        Gets the IDs of the empty lumen
        Empty means that the area is smaller than the tube_radius^2
            
            Input
            -----
            self : network object
                Needs to be called by a class object
            
            Returns
            -------
            empty_lumen_list : int list
                Contains the IDs of the empty lumens
            """

        empty_lumen_list = []
        # For each lumen ID
        for i in range(len(self.area)):
            # If area is smaller than the treshhold
            if(self.area[i] < self.tube_radius ** 2 and self.empty_list[i] == 0):
                self.empty_list[i] = 1
                self.area[i] = 0
                empty_lumen_list.append(i)
        return empty_lumen_list


    def get_neighbours(self, lumen):
        """
        Gets the lumen that are directly connected to the lumen passed on and deletes the connections
            
            Input
            -----
            self : network object
                Needs to be called by a class object
            lumen : int
                ID of a lumen
            
            Returns
            -------
            neighbour_list : int list
                ID of all lumen that are directly connected to the lumen passed on
            """
        neighbour_list = []
        line = 0
        # Going through links in lumen_lumen.dat
        while line < len(self.lumen_lumen) and self.lumen_lumen[line][0] < lumen :
            if self.lumen_lumen[line][1] == lumen :
                neighbour_list.append([self.lumen_lumen[line][0], self.lumen_lumen[line][2]])

                event = 'link lumen %d to lumen %d removed' % (self.lumen_lumen[line][0], lumen)
                self.save_event(event)
        
                link = [self.lumen_lumen[line][0], self.lumen_lumen[line][1], self.lumen_lumen[line][2]]
                self.lumen_lumen.remove(link)
            
            else : line += 1
        while line < len(self.lumen_lumen) and self.lumen_lumen[line][0] < lumen : 
            line += 1
        while(line < len(self.lumen_lumen) and self.lumen_lumen[line][0] == lumen):
            neighbour_list.append([self.lumen_lumen[line][1], self.lumen_lumen[line][2]])

            event = 'link lumen %d to lumen %d removed' % (lumen, self.lumen_lumen[line][1])
            self.save_event(event)
            
            link = [self.lumen_lumen[line][0], self.lumen_lumen[line][1], self.lumen_lumen[line][2]]
            self.lumen_lumen.remove(link)

        return neighbour_list
        
    def get_bridges(self, lumen):
        """
        Gets the bridges that are directly connected to the lumen passed on
            
            Input
            -----
            self : network object
                Needs to be called by a class object
            lumen : int
                ID of a lumen
            
            Returns
            -------
            neighbour_list : int list
                ID of all lumen that are directly connected to the lumen passed on
            """
        bridge_list = []
        line = 0
        # Going through the links in bridge_lumen.dat
        while(line < len(self.bridge_lumen)):
            if (self.bridge_lumen[line][1] == lumen):
                bridge_list.append([self.bridge_lumen[line][0], self.bridge_lumen[line][2]])
                event = 'link bridge %d to lumen %d removed' % (self.bridge_lumen[line][0], lumen)
                self.save_event(event)
               
                self.bridge_lumen.remove(self.bridge_lumen[line])
            else: line += 1
        return bridge_list

############################################################################################################################
#################################################### Creating Functions ###################################################
############################################################################################################################



    def create_link(self, link):
        """
        Creates a link between two lumen in lumen_lumen.dat
            
            Input
            -----
            self : network object
                Needs to be called by a class object
            link : float array
                [ID lumen1, ID lumen2, length]
            
            Returns
            -------
            no return
            """
        # no self-loops allowed
        if(len(link) == 4 and link[0] != link[1]):
            # Ensuring: lumen_1 < lumen_2
            if(link[0] < link[2]):
                lumen_1 = link[0]
                lumen_2 = link[1]
            else:
                lumen_1 = link[1]
                lumen_2 = link[0]
            length = link[2]
            line = 0
            # Finding line in lumen_lumen.dat, to keep the sorting
            while(line < len(self.lumen_lumen) and lumen_1 > self.lumen_lumen[line][0]): line += 1
            if(line < len(self.lumen_lumen) - 1):
                while(line < len(self.lumen_lumen) and lumen_2 > self.lumen_lumen[line][1] and lumen_1 == self.lumen_lumen[line][0]): line += 1

            # Creating the link in lumen_lumen.dat
            
            self.lumen_lumen.append([lumen_1,lumen_2, length])
            self.lumen_lumen.sort()
            
            event = 'link lumen %d to lumen %d created' % (lumen_1,lumen_2)
            self.save_event(event)
        return;

    def create_bridge_lumen(self, link):
        """
        Creates a link between a lumen and a bridge in bridge_lumen.dat
            
            Input
            -----
            self : network object
                Needs to be called by a class object
            link : float array
                [ID bridge, ID lumen, length]
            
            Returns
            -------
            no return
            """
        bridge = link[0]
        lumen = link[1]
        length = link[2]
        line = 0

        # Creating the link in bridge_lumen.dat
        self.bridge_lumen.append(link)
        
        self.bridge_lumen.sort()
        
        event = 'link bridge %d to lumen %d created' % (bridge,lumen)
        self.save_event(event)
        return;

    def create_bridge_bridge(self, link):
        """
        Creates a link between two bridges in bridge_bridge.dat
            
            Input
            -----
            self : network object
                Needs to be called by a class object
            link : float array
                [ID bridge1, ID bridge2, length]
            
            Returns
            -------
            no return
            """
        if(link[0] == link[1]): return;
        if(link[0] < link[1]):
            bridge_1 = link[0]
            bridge_2 = link[1]
        else:
            bridge_1 = link[1]
            bridge_2 = link[0]
        length = link[2]
        line = 0

        # Creating the link in bridge_bridge.dat
        self.bridge_bridge.append([bridge_1,bridge_2, length])
        
        self.bridge_bridge.sort()
        
        event = 'link bridge %d to bridge %d created' % (bridge_1,bridge_2)
        self.save_event(event)
        return;


    def create_bridge(self, lumen, bridge, lumid):
        """
        Creates a new bridge connected with the lumen and bridges passed on
            
            Input
            -----
            self : network object
                Needs to be called by a class object
            lumen : int list
                [[lumen ID, length], [lumen ID, length],.....]
                lumen IDs to which the new bridge should be connected to
            bridge : int list
                [[bridge ID, length], [bridge ID, length],.....]
                bridge IDs to which the new bridge should be connected to
            
            Returns
            -------
            no return
            """
        #####
        bridge_conversionfile = os.path.join(self.network_folder,'bridgesconversion.txt')
        
        # ID of the new bridge
        bridge_number = self.num_bridges
        # Bridge ID counter, contains the ID of the next new bridge
        self.num_bridges += 1
        event = 'new bridge %d' % (bridge_number) + ' (' + str(lumid) + ')'
        self.save_event(event)
        line = 0
        lumen.sort()
        bridge.sort()
        
        # For each lumen that should be connected to the new bridge
        for i in range(len(lumen)):
            new_link = [bridge_number, lumen[i][0], lumen[i][1]]
            # Create link in bridge_lumen.dat
            self.create_bridge_lumen(new_link)
        
        # For each lumen that should be connected to the new bridge
        for i in range(len(bridge)):
            new_link = [bridge[i][0], bridge_number, bridge[i][1]]
            # Create link in bridge_bridge.dat
            self.create_bridge_bridge(new_link)
    
        open(bridge_conversionfile, 'a').write(str(bridge_number) + ' ' + str(lumid)+ '\n')
        return;


############################################################################################################################
################################ Geometric Functions for area and Pressure ###############################################
############################################################################################################################

    def set_theta(self):
        """
        Sets the angle theta
        Calculates the angle theta, angle between the lumen and the tube
            
            Input
            -----
            self : network object
                Needs to be called by a class object
            
            Returns
            -------
            theta : float list  
                Theta value for each lumen
            """
        theta = []
        
        for i in range(len(self.alpha)):
            
            #cos = (2*self.alpha[i]-(4*self.alpha[i]**2-self.delta[i]**2+1)/(4*self.alpha[i]))/self.delta[i] ## Old version, for assymmetric lumen
            #theta.append(math.acos(cos))
            
            theta.append(np.arccos(self.alpha[i]))
            
        return theta;
    
    def set_area_factor(self):
        """
        Sets the area factor, needed to express the pressure in terms of the area instead of the curvature radius
            
            Input
            -----
            self : network object
                Needs to be called by a class object
            
            Returns
            -------
            area_factor : float list
                Area factor for each lumen
            """
        area_factor = []
        for i in range(len(self.alpha)):
            area_factor.append(np.sqrt((2*self.theta[i]-np.sin(2*self.theta[i]))))
        return area_factor;
    
    def opening_radius(self, lumen):
        """
        Calculates the length/2 parallel to the 'tube' where the membrane is not attached for a given lumen
            
            Input
            -----
            lumen : int
                ID of the lumen
            
            Returns
            -------
            radius : float
                Length/2 of the opening radius
            """
        return np.sqrt(2*self.area[lumen]/(2*self.theta[lumen]-np.sin(2*self.theta[lumen])))*np.sin(self.theta[lumen])

    def get_area(self, lumen):
        """
        Calculates the area in one half of the lumen (for symmetric lumen)
            
            Input
            -----
            lumen : int
                ID of the lumen
            
            Returns
            -------
            area : float
                Area/2 of the lumen
            """
        area = self.area[lumen]
        return area

    def pressure(self,lumen):
        """
        Calculates the pressure inside the lumen (for symmetric lumen)
            
            Input
            -----
            lumen : int
                ID of the lumen
            
            Returns
            -------
            pressure : float
                Pressure of the lumen
            """
        
        area = self.get_area(lumen)
        # Avoid dividing by zero
        if(area < 0.1 * self.tube_radius**2 ):
            error = 'division by zero in pressure: lumen ID: %d' % (lumen)
            self.save_error(error)
        pressure = self.gamma_lumen[lumen]*self.area_factor[lumen]/np.sqrt(area)
        return pressure


############################################################################################################################
################################################# Reading Functions ########################################################
############################################################################################################################

    def read_lumen_lumen(self, lumen_lumen_file):
        """
        Reading the file with links between two lumens
            
            Input
            -----
            lumen_lumen_file : str
                File path to file with the links between two lumens
            
            Returns
            -------
            lumen_lumen : float list  [lumen1, lumen2, length]
                Information about the links between two lumens
            """
        
        if (os.path.getsize(lumen_lumen_file)>0): # If the file is not empty
            lumen_1, lumen_2 = np.loadtxt(lumen_lumen_file, dtype = int, usecols = [0,1], unpack = True)
            length = np.loadtxt(lumen_lumen_file, dtype = float, usecols = [2])
            lumen_lumen = np.column_stack([lumen_1, lumen_2, length]).tolist()
        
        else:
            lumen_lumen = []
        return lumen_lumen


    def read_bridge_lumen(self, bridge_lumen_file):
        """
        Reading the file with links between bridge and lumen
            
            Input
            -----
            bridge_lumen_file : str
                File path to file with the links between bridge and lumen
            
            Returns
            -------
            bridge_lumen : float list [bridge, lumen, length]
                Information about the links between bridge and lumen
            num_bridges : int
                Number of bridge_lumen links
            """
        
        with open(bridge_lumen_file, 'r') as f:
            lines = f.read().splitlines()
            last_line = lines[-1]
        if ('#' in last_line): # If the file is empty
            bridge_lumen = []
            num_bridges = 0 # number of existing bridges
        else:
            bridge, lumen = np.loadtxt(bridge_lumen_file, dtype = int, usecols = [0,1], unpack = True)
            length = np.loadtxt(bridge_lumen_file, dtype = float, usecols = [2])
            bridge_lumen = np.column_stack([bridge, lumen, length]).tolist()
            num_bridges = max(bridge)+1 # number of existing bridges
        return bridge_lumen, num_bridges

    def read_bridge_bridge(self, bridge_bridge_file, num_bridges):
        """
        Reading the file with links between two bridge
            
            Input
            -----
            bridge_bridge_file : str
                File path to file with the links between two bridge
            
            Returns
            -------
            bridge_bridge : float list  [bridge1, bridge2, length]
                Information about the links between two bridge
            num : int
                Number of bridge_bridge links
            """
        with open(bridge_bridge_file, 'r') as f:
            lines = f.read().splitlines()
            last_line = lines[-1]
        if ('#' in last_line>0): # If the file is empty
            bridge_bridge = []
            num = num_bridges
        else:
            bridge1, bridge2 = np.loadtxt(bridge_bridge_file, dtype = int, usecols = [0,1], unpack = True)
            length = np.loadtxt(bridge_bridge_file, dtype = float, usecols = [2])
            bridge_bridge = np.column_stack([bridge1, bridge2, length]).tolist()
            if (max(bridge2)+1 > num_bridges): num = max(bridge2)+1

        return bridge_bridge, num

############################################################################################################################
################################################# Output functions #########################################################
############################################################################################################################
    def save_event(self, event, start = False, out_path = ''):
        """
        Saves each event in the output folder in the file event.dat
        Events like a lumen disappearing, reconnections in the graph
            
            Input
            -----
            event : str
                Message of the event
            start : boolean
                True: File is created
                False: the message is stored in the file
            Returns
            ------
            no return
            """
        if(start):
            header_event = '# Saves each event during the simulation; event is a disappearing lumen, graph reconnection \n'
            self.file_event = os.path.join(out_path, 'event.dat')
            fevent = open(self.file_event, 'w')
            fevent.write(header_event)
            fevent.close()
        else:
            fevent = open(self.file_event, 'a')
            fevent.write('%.5f' % self.current_time)
            fevent.write(' ')
            fevent.write(event)
            fevent.write('\n')
            fevent.close()
        return;

    def save_error(self, error, start = False, out_path = ''):
        """
        Saves errors in the output folder in the file error.dat
        Errors like volume loss
            
            Input
            -----
            error : string
                Message of the event
            start : boolean
                True: File is created
                False: the message is stored in the file
            Returns
            ------
            no return
            """
        if(start):
            header_error = '# Saves each warning like volume loss \n'
            self.file_error = os.path.join(out_path, 'error.dat')
            ferror = open(self.file_error, 'w')
            ferror.write(header_error)
            ferror.close()
        else:
            ferror = open(self.file_error, 'a')
            ferror.write('%.5f' % self.current_time)
            ferror.write(' ')
            ferror.write(error)
            ferror.write('\n')
            ferror.close()
        return;


    def save_area(self, start = False, out_path = ''):
        """
        Saves the volume evolution in the output folder in the file area.dat
            
            Input
            -----
            start : boolean
                True: File is created
                False: the message is stored in the file
            Returns
            ------
            no return
            """
        if(start):
            header_volume = '# Saves the volume evolution of each lumen for the time step %f \n' %(self.time_step)
            self.file_area = os.path.join(out_path, 'area.dat')
            farea = open(self.file_area, 'w')
            farea.write(header_volume)
            farea.close()
            self.t_old = 0
        else:
            farea = open(self.file_area, 'a')
            farea.write('%.5f' % self.current_time)
            farea.write(' ')
            farea.write(' '.join(map(str, self.area)))
            farea.write('\n')
            farea.close()
        return;


############################################################################################################################
################################################# Swelling functions #######################################################
############################################################################################################################


    def swelling(self, lumen) :
        """
        self.swelling(lumen)
        
            Calculates the input flux for the area fo a given lumen, due to swelling.
        
            Input
            -----
            lumen : int
                Index of the lumen
        
        """
        area = self.get_area(lumen)
        theta = self.theta[lumen]
        
        flux_swelling = self.swelling_rate * 4 * theta * np.sqrt(area)/ self.area_factor[lumen]
        #print flux_swelling
        return flux_swelling

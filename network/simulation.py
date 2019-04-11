#!/usr/bin/env python2.7
# simulation.py

"""
    simulation.py conf.init
    
    Python script for simulating a network of lumens etc.
    
    Arguments
    ---------
    conf.init: the initial network configuration file
    
    Options
    -------
    swelling : boolean, optional, default : False
        Include swelling or not to the simulation
    
    Examples:
    
    1. network_simulation.py conf.init
    Submit one job to run the config file provided
    
    ...

    
    Functions
    ---------
    initialize_area
    initialize_tension
    end_simulation_odeint
    end_simulation_ivp
    give_times
    non_empty_lumens
    save_final_state
    run_integration
    simulation

    main
    
    
    Needed
    ------
    Libraries : numpy, os, sys, scipy.integrate, configparser
    Scripts : pattern, network
    
    Created by A. Mielke
    Modified by H. Turlier on 8/06/2018
    Modified by M. Le Verge--Serandour on 08/04/2019
    """

import numpy as np
import os
import network as net
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
import configparser
import sys


def non_empty_lumens(system) :
    """
    non_empty_lumens(system)
    
    Gives the number of non empty lumens (fulflling the conditionthat A > d*d)
    
        Input
        -----
        system : network object
    
        Returns
        -------
        n : int
            Number of non empty lumens
    """
    n = len(system.area)-np.count_nonzero(system.empty_list)
    return n

def give_times(system, nint = 20, cst = 1e-2) :
    """
    give_times(system, nint = 20, cst = 1e-2)
    
    Give the time step and ending time for the integration
    
        Input
        -----
        system : network object
    
        nint : int, optional, default : 20
            Number of intervals between 
        cst : float, optional, default : 1e-2
            An arbitrary constant for the integration
    
        Returns
        -------
        time_step : float
            cst / (n*n)
        time_end : float
            nint * time_step
    """
    n = len(system.area)-np.count_nonzero(system.empty_list)
    #print len(system.area), system.empty_list
    a_tot = np.sum(system.area)
    a_mf = a_tot / n

    kappa, alpha = 1e-2, 2.
    time_step = kappa*a_mf**alpha
    
    t_end = nint * time_step
    
    return time_step, t_end

def initialize_area(number_lumen, config):
    """
    Initializes the area with a given distribution
        
    Possible distributions: gaussian

        Input
        -----
        number_lumen : int
            Number of lumen in the system
        config : dictionary
            Contains all the information of the configuration file

        Returns
        -------
        init_area : list of floats
            List with the area distribution.
        """
    area_distribution = config['volume']['distribution']
    if(area_distribution == 'gaussian'):
        vol_mean = float(config['volume']['mean'])
        vol_sigma = float(config['volume']['sigma'])
        init_area = pattern.random_area(number_lumen, vol_mean, vol_sigma)
    else :
        print 'Non recognized distribution !'
        return;
    return init_area

def initialize_tension(init_area, network_folder, config):
    """
    Initializes the tension and saves the properties of the lumen in lumen.dat
        
        Input
        -----
        network_folder : string
            Path to the folder containing the files about the network
        config : dictionary
            Contains all the information of the configuration file
        init_area : int list
            Initial area distribution
        Returns
        -------
        no return
        """
    tension_type = config['tension']['type']
    gamma = float(config['tension']['gamma_cell'])
    gamma_contact = float(config['tension']['gamma_contact'])
    # create lumen file containing information about lumen
    if(tension_type == 'TM'):
        delta = float(config['tension']['delta'])
        pattern.gamma_TM(network_folder, init_area, gamma, delta, gamma_contact)
    if(tension_type == 'TM_contact'):
        delta = float(config['tension']['delta'])
        pattern.gamma_contact_TM(network_folder, init_area, gamma, gamma_contact, delta)
    elif(tension_type == 'constant'):
        pattern.lumen_equal_gamma(network_folder, init_area, gamma, gamma_contact)
    return;

def end_simulation_odeint(t,vol, system, time_step):
    """
    Stopping condition for integration with odeint
        
        Input
        -----
        t : float
            Current time step
        vol : float list
            Current area distribution
        time_step : float
            Time step of the simulation (absolute units)
    
        Returns
        -------
        b : boolean
            True if the simulation continues ; False if the Simulation is stopped
        """
    b = True
    c = 0
    for v in vol:
        if(v > system.tube_radius ** 2): c+=1           # counts the number of non-empty lumens
    if (c < 2): return False
    if(system.end_time > 0 and system.end_time < t - 2 * time_step):
        b = False
        system.final_area = vol
    return b

def end_simulation_ivp(t, vol):
    """
    Stopping condition for integration with solve_ivp
        
        Input
        -----
        t : float
            Current time step
        vol : float list
            Current area distribution
    
        Returns
        -------
        b : boolean
            True if the simulation continues ; False if the Simulation is stopped
        """
    b = True
    global system, time_step

    if(t - system.end_time > 2*time_step  and system.end_time > 0):
        b = False
        system.final_volume = vol
    return b

def save_final_state(init_area, system, out_folder):
    """
    Saves initial and final state
        
        Input
        -----
        init_area : float list
            Initial area distribution of the lumens
        system : network object
    
        out_folder : str
            Path of the output folder
    
        Returns
        -------
        no return
        """
    
    log_file = configparser.ConfigParser()
    log_file['volume'] = {}
    log_file['volume']['total'] = str(np.sum(init_area))
    log_file['volume']['final'] = str(np.sum(system.final_area))
    log_file['volume']['initial'] = ' '.join(str(e) for e in init_area)
    ID = []
    log_file['lumen'] = {}
    log_file['lumen']['ID'] = str(np.argmax(system.final_area))
    log_file['time'] = {}
    log_file['time']['end'] = str(system.end_time)
    
    with open( os.path.join(out_folder,'area.log'), 'w') as configfile:
        log_file.write(configfile)
    return;

def run_integration(init_area, config, system):
    """"
    Runs the integration
        
        Input
        -----
        init_area : int list
            Initial area distribution of the lumens
        config : dictionary
            Contains all the information of the configuration file
        system : network object
        
        Returns
        -------
        no return
        """


    def flux_odeint(vol, t):
        """
            Changes Input order for the integrator odeint
            
            Input
            -----
            vol : float list
                Current area distribution
            t : float
                Current time step
        
            Returns
            -------
            float array     
                Area change from network.py
            """
        return system.flux(t, vol)


    # initializing the time steps for the simulation
    global time_step
    time_max = float(config['time']['time_max'])
    time_step = float(config['time']['time_step'])
    
    a_tot = np.sum(system.area)
    
    # adaptative time step
    if(config.has_option('integration', 'adaptative')) : 
        adaptative = eval(config['integration']['adaptative'])
    else : 
        adaptative = False
    
    
    if(config.has_option('integration', 'min_step')) :
        min_step = float(config['integration']['min_step'])
    else : 
        min_step = 0.
    
    if(config.has_option('integration', 'max_step')) :
        step_max = float(config['integration']['max_step'])
    else :
        step_max = 0.1
    
    # define the integrator ; default is odeint
    if(config.has_option('integration', 'integrator')) :
        integrator = str(config['integration']['integrator'])
    else : integrator = 'odeint'
    
    # ===================== run the simulation
    
    if(integrator == 'ivp') :
        if(config.has_option('integration', 'method')) :
            method = str(config['integration']['method'])
        else :
            method = 'LSODA'

        # start Simulation with solve_ivp
        volume = solve_ivp(system.flux, time_span, init_area, method = method, events = end_simulation_ivp, min_step = min_step, max_step = step_max)
    
    else:
        # start Simulation using odeint
        sim = True
        nint = 1
        cst = 0.1
        
        
        time_start = 0

        vol_start = init_area
        
        
        if adaptative :
            time_step, t_end = give_times(system, nint = nint, cst = cst)
            time_end = t_end
            
            n = non_empty_lumens(system)
            
            while(sim and time_start < time_max):
                n = non_empty_lumens(system)                                                        # counts the number of non-empty lumens
                
                area = odeint(flux_odeint, vol_start, [time_start, time_end], hmax=step_max)        # Calculate the area changes
                sim = end_simulation_odeint(time_end, area[-1,:], system, time_step)                # Test the ending conditions
                time_step, t_end = give_times(system, nint = nint, cst = cst)
                
                vol_start = system.area
                time_start = time_end
                time_end += t_end
            
        else :
            while(sim and time_start < time_max):
                time_end = time_start + time_step
                
                area = odeint(flux_odeint, vol_start, [time_start, time_end], hmax=step_max)        # Calculate the area changes
                sim = end_simulation_odeint(time_end, area[-1,:], system, time_step)                # Test the ending conditions
                vol_start = system.area
                time_start = time_end
                time_end += nint * time_step
    return;


def simulation(config):
    global system
    
    print "Simulation running..."

    network_folder  = config['network']['path']

    try :
        out_folder = os.path.join(network_folder, 'out')
        os.mkdir(out_folder)
    except :
        pass
    os.chmod(out_folder, 504)
    out_folder = os.path.join(network_folder, 'out')
    try :
        os.mkdir('network')
    except :
        pass
    os.chmod('network', 504)
    
    swelling        = eval(config['swelling']['swelling_bool'])
    swelling_rate   = float(config['swelling']['swelling_rate'])
    
    save_area       = eval(config['integration']['save_area'])

    if len(network_folder) > 0 :
        network_folder = os.path.join(network_folder, 'network')
        # Initializing the area
        init_area = np.loadtxt(os.path.join(network_folder, 'lumen.dat'), unpack=True, usecols = (3))
        # Initializing the number of lumens
        number_lumen = len(init_area)
        # Initializing the number of edges
        number_edges = len(np.loadtxt(os.path.join(network_folder, 'lumen_lumen.dat'), unpack=True, usecols = (0)))
    else :
        print 'Configuration files not found, please give the folder path.'
        return;
        
    # initializing the network
    system = net.network(network_folder, out_folder, float(config['time']['time_step']), tube_radius=float(config['network']['tube_radius']), friction=float(config['network']['friction']), swelling = swelling, swelling_rate=swelling_rate, save_area_dat=save_area)
    
    # start Simulation
    run_integration(init_area, config, system)

    # Saving the final state (winning lumen, area)
    save_final_state(init_area, system, out_folder)

    return;
#------------------------------------------------------------------------

def main(conf_name):
    """return the parsed simulation configuration"""
    config = configparser.ConfigParser()
    config.read(conf_name)
    simulation(config)

#------------------------------------------------------------------------

if __name__ == "__main__" :
    if len(sys.argv) < 2:
        print('[network_simulation.py] Error: missing config file, type help for instructions.')
    elif sys.argv[1]=='help':
        print(__doc__)
    # first argument should be a readable config file:
    elif os.access(sys.argv[1], os.R_OK):
        conf_name = sys.argv[1]
        main(conf_name)
    else :
        print 'Error : no readable config file found'
        print sys.argv[1]
    sys.exit()

#------------------------------------------------------------------------
#------------------------------------------------------------------------

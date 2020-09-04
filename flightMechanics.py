# -*- coding: utf-8 -*-
"""
   FLIGHT MECHANICS
"""
import os
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, Poly
from scipy import signal
from scipy.optimize import fsolve
from scipy.interpolate import interp1d

# Import utilities
sys.path.insert(1, '/Users/Public/Python/Utilities')
from atmos import Atmos
from speed import Speed

def read_input_file(file):
    """Function to read input file into a dictionary

    Args:
        file (str): Path to the input file.

    Returns:
        input_data (dict): Dictionary with the input file variables.
        adf_data (dict): Dictionary with the aircraft data for the selected condition.
        aircraft_config (dict): Dictionary with the aircraft configuration.

    """
    with open(file,'r') as f:
        lines = [line.strip().split(':') for line in f.readlines() if (line.strip() and not line.startswith("#"))]
        
    variable_name = [line[0].strip() for line in lines]
    variable_data = [line[1].replace(" ","").split(',') for line in lines]
    
    input_data = {}   
    for index, title in enumerate(variable_name):    
        if title == 'aircraft':
            file_path = os.path.split(file)
            input_data[title] = file_path[0] + '/' + variable_data[index][0]
        elif title in ['condition', 'type', 'plot']:
            input_data[title] = np.array([int(var) for var in variable_data[index]])[0]
        elif title in ['CG']:
            input_data[title] = np.asarray(variable_data[index], dtype=np.float64, order='C')[0]
        else:
            input_data[title] = np.array([float(var) for var in variable_data[index]])

    if 'plot' not in input_data:
        input_data['plot'] = 1

    # Read adf data            
    adf_data = read_adf_file(input_data)
    
    # Read aircraft config
    aircraft_config = read_aircraft_config(adf_data)
    
    return input_data, adf_data, aircraft_config

def read_adf_file(input_data):
    """Function to read adf file into a dictionary

    Args:
        input_data (dict): Dictionary with the input file variables.

    Returns:
        adf_data (dict): Dictionary with the aircraft data for the selected condition.

    """
    file = input_data['aircraft']
    condition = input_data['condition']

    ind = condition - 1

    with open(file,'r') as f:
        lines = [line.strip().split(':') for line in f.readlines() if (line.strip() and not line.startswith("#"))]
        
    variable_name = [line[0].strip() for line in lines]
    variable_data = [line[1].replace(" ","").split(',') for line in lines]

    adf_data = {}
    common_str    = ['A', 'B', 'C', 'DDD', 'E', 'F', 'V']
    common_single = ['n_conditions','TO_condition', 'LD_condition',
                     'TOW_min', 'TOW_max', 'LW_min', 'LW_max', 'dW', 'delta_e_max', 'delta_a_max', 'delta_r_max',
                     'S', 'c', 'b', 'n_eng', 'mu_gr', 'CA_x', 'CG_y', 'CG_z']
    common_array  = ['TOGA', 'epsilon', 'ni', 'eng_x', 'eng_y', 'eng_z', 'lg_x', 'lg_y', 'lg_z']
    cond_filter   = ['ISA', 'altitude', 'speed', 'CG_position', 'm', 'Ix', 'Iy', 'Iz', 'Jxz']                 
    T_array       = ['T' + str(i) for i in range(20)]
    T_condition   = 'T' + str(condition)
    for index, title in enumerate(variable_name):
        if title in common_str: # Common data in str format
            adf_data[title] = variable_data[index][0]
        elif title in common_single: # Common data as single value
            adf_data[title] = np.asarray(variable_data[index], dtype=np.float64, order='C')[0]
        elif title in common_array: # Common data in an array
            adf_data[title] = np.asarray(variable_data[index], dtype=np.float64, order='C')
        elif title in cond_filter: # Select data for the desired condition
            adf_data[title] = np.asarray(variable_data[index], dtype=np.float64, order='C')[ind]
        elif title in T_array: # Select thrust por the desired condition
            if title == T_condition:
                adf_data['T'] = np.asarray(variable_data[index], dtype=np.float64, order='C')
        else: # Get coeffs for desired condition, TO and LD (if exist)
            adf_data[title] = np.asarray(variable_data[index], dtype=np.float64, order='C')[ind]
            if adf_data['TO_condition'] != 0:
                ind_TO = int(adf_data['TO_condition'] - 1)
                adf_data[title + '_TO'] = np.asarray(variable_data[index], dtype=np.float64, order='C')[ind_TO]
            if adf_data['LD_condition'] != 0:
                ind_LD = int(adf_data['LD_condition'] - 1)
                adf_data[title + '_LD'] = np.asarray(variable_data[index], dtype=np.float64, order='C')[ind_LD]
    
    adf_data['n_conditions'] = int(adf_data['n_conditions'])
    adf_data['TO_condition'] = int(adf_data['TO_condition'])
    adf_data['LD_condition'] = int(adf_data['LD_condition'])
    adf_data['n_eng']        = int(adf_data['n_eng'])

    # Get speed in Mach, KCAS and KTAS
    if adf_data['speed'] < 5.0 :
        adf_data['Mach'] = adf_data['speed']
        adf_data['KCAS'] = Speed.mach_to_cas(adf_data['Mach'], adf_data['altitude']*feet_to_meters) / kt_to_ms
        adf_data['KTAS'] = Speed.mach_to_tas(adf_data['Mach'], adf_data['altitude']*feet_to_meters, adf_data['ISA']) / kt_to_ms
    else:
        adf_data['KCAS'] = adf_data['speed']
        adf_data['KTAS'] = Speed.cas_to_tas(adf_data['KCAS']*kt_to_ms, adf_data['altitude']*feet_to_meters, adf_data['ISA']) / kt_to_ms
        adf_data['Mach'] = Speed.cas_to_mach(adf_data['KCAS']*kt_to_ms, adf_data['altitude']*feet_to_meters)
        
    # Get Cz_s for the stability condition
    rho = Atmos.get_density(feet_to_meters * adf_data['altitude'], adf_data['ISA'])
    p_d = get_dynamic_pressure(rho, kt_to_ms * adf_data['KTAS'])
    adf_data['Cz_s'] = -(adf_data['m'] * g) / (p_d * adf_data['S'])

    # Get parameters affected by CG position
    if 'CG' in input_data:
        adf_data['CG'] = input_data['CG']
    else:
        adf_data['CG'] = adf_data['CG_position']
    
    get_cg_parameters(adf_data)
    
    return adf_data

def read_aircraft_config(adf_data):
    """Function to read the aircraf configuration into a dictionary

    Args:
        adf_data (dict): Dictionary with the aircraft data for the selected condition.

    Returns:
        aircraft_config (dict): Dictionary with the aircraft configuration.

    """

    aircraft_config = {}

    component = ['Fuselage', 'Wing', 'Tail', 'Prop_dist', 'Prop_tip', 'Prop_tail','Air_inlet', 'Canard']
    option = [['Streamlined', 'Conventional'],                 # Fuselage
              ['Low', 'High'],                                 # Wing
              ['T', 'V', 'U', 'Conventional'],                 # Tail
              ['Lower', 'Upper', 'Zigzag', 'Embedded', 'No'],  # Prop_dist
              ['Puller', 'Pusher', 'No'],                      # Prop_tip
              ['BLI', 'No'],                                   # Prop_tail
              ['Upper', 'Lower', 'Lateral', 'No'],             # Air_inlet
              ['Yes', 'No']]                                   # Canard
    code   = [['1', '0'],                # Fuselage
              ['1', '0'],                # Wing
              ['1', '2', '3', '0'],      # Tail
              ['1', '2', '3', '4', '0'], # Prop_dist
              ['1', '2', '0'],           # Prop_tip
              ['1', '0'],                # Prop_tail
              ['1', '2', '3', '0'],      # Air_inlet
              ['1', '0']]                # Canard

    letter = ['A', 'B', 'C', 'DDD', 'E', 'F']
    
    i = 0
    for let in letter:
        ac_code = adf_data[let]
        if let == 'DDD': # Prop
            ac_code = list(ac_code)
            for c in ac_code:
                code_ind  =  code[i].index(c)
                aircraft_config[component[i]] = option[i][code_ind]
                i += 1
        else:
            code_ind  =  code[i].index(ac_code)
            aircraft_config[component[i]] = option[i][code_ind]
            i += 1

    if 'V' in adf_data:
        aircraft_config['Version'] = adf_data['V']
    
    return aircraft_config

def get_cg_parameters(adf_data):
    """Function to obtain all the parameters afected by CG position

    Args:
        adf_data (dict): Dictionary with the aircraft data for the selected condition.

    Returns:
        adf_data (dict): Dictionary with the aircraft data for the selected condition,
                         adding parameters afected by CG position.

    """

    # Center of gravity
    adf_data['CG_x'] = adf_data['CA_x'] + (adf_data['CG'] - 0.25) * adf_data['c']

    # Engine distances in body axes
    adf_data['CG_eng_x'] = -(adf_data['eng_x'] - adf_data['CG_x'])
    adf_data['CG_eng_y'] = adf_data['eng_y'] - adf_data['CG_y']
    adf_data['CG_eng_z'] = -(adf_data['eng_z'] - adf_data['CG_z'])

    # Landing gear distances in body axes
    adf_data['CG_lg_x'] = -(adf_data['lg_x'] - adf_data['CG_x'])
    adf_data['CG_lg_y'] = adf_data['lg_y'] - adf_data['CG_y']
    adf_data['CG_lg_z'] = -(adf_data['lg_z'] - adf_data['CG_z'])

    # Aerodynamic coefficients (Only the ones affected)
    Cm_0 = adf_data['Cm_0'] - (adf_data['CG'] - adf_data['CG_position'])*adf_data['Cz_0']
    Cm_alpha = adf_data['Cm_alpha'] - (adf_data['CG'] - adf_data['CG_position'])*adf_data['Cz_alpha']
    Cz_q = adf_data['Cz_q'] - (adf_data['CG'] - adf_data['CG_position'])*adf_data['Cz_q']
    Cx_q = adf_data['Cx_q'] - (adf_data['CG'] - adf_data['CG_position'])*adf_data['Cx_q']
    Cm_q = adf_data['Cm_q'] - (adf_data['CG'] - adf_data['CG_position'])*adf_data['Cz_q'] + (adf_data['CG'] - adf_data['CG_position'])**2 * adf_data['Cz_q']
    Cm_u = adf_data['Cm_u'] - (adf_data['CG'] - adf_data['CG_position'])*adf_data['Cz_u']
    Cm_delta_e = adf_data['Cm_delta_e'] - (adf_data['CG'] - adf_data['CG_position'])*adf_data['Cz_delta_e']
    
    Cn_beta = adf_data['Cn_beta'] + (adf_data['CG'] - adf_data['CG_position'])*adf_data['c']/adf_data['b'] * adf_data['Cy_beta']
    Cn_p = adf_data['Cn_p'] + (adf_data['CG'] - adf_data['CG_position'])*adf_data['c']/adf_data['b'] * adf_data['Cy_p']
    Cy_r = adf_data['Cy_r'] + (adf_data['CG'] - adf_data['CG_position'])*adf_data['c']/adf_data['b'] * adf_data['Cy_r']**2 / adf_data['Cn_r']
    Cl_r = adf_data['Cl_r'] + (adf_data['CG'] - adf_data['CG_position'])*adf_data['c']/adf_data['b'] * adf_data['Cy_r'] +\
                       (adf_data['CG'] - adf_data['CG_position'])*adf_data['c']/adf_data['b'] * adf_data['Cy_r']*adf_data['Cl_r']/adf_data['Cn_r'] +\
                       ((adf_data['CG'] - adf_data['CG_position'])*adf_data['c']/adf_data['b'])**2 * adf_data['Cy_r']**2/adf_data['Cn_r']
    Cn_r = adf_data['Cn_r'] + (adf_data['CG'] - adf_data['CG_position'])*adf_data['c']/adf_data['b'] * adf_data['Cy_r'] +\
                       (adf_data['CG'] - adf_data['CG_position'])*adf_data['c']/adf_data['b'] * adf_data['Cy_r']**2 +\
                       ((adf_data['CG'] - adf_data['CG_position'])*adf_data['c']/adf_data['b'])**2 * adf_data['Cy_r']**3/adf_data['Cn_r']
    Cn_delta_r = adf_data['Cn_delta_r'] + (adf_data['CG'] - adf_data['CG_position'])*adf_data['c']/adf_data['b'] * adf_data['Cy_delta_r']

    adf_data['Cm_0'] = Cm_0
    adf_data['Cm_alpha'] = Cm_alpha
    adf_data['Cz_q'] = Cz_q
    adf_data['Cx_q'] = Cx_q
    adf_data['Cm_q'] = Cm_q
    adf_data['Cm_u'] = Cm_u
    adf_data['Cm_delta_e'] = Cm_delta_e

    adf_data['Cn_beta'] = Cn_beta
    adf_data['Cn_p']= Cn_p
    adf_data['Cy_r']= Cy_r
    adf_data['Cl_r']= Cl_r
    adf_data['Cn_r']= Cn_r
    adf_data['Cn_delta_r'] = Cn_delta_r 
    

    # Moment of inertia ??

def get_dynamic_pressure(density, speed):
    """Function to calculate the dynamic pressure for a given density and speed

    Args:
        density (float): Air density (kg/m3).
        speed (float): True air speed (m/s).

    Returns:
        p_d (float): Dynamic pressure (N/m2).

    """
    p_d = 0.5*density*speed ** 2
    
    return p_d

def get_stall_speed(m, density, S, CL_max):
    """Function to calculate the reference stall speed (Vsr)

    Args:
        m (float): Weight (kg)
        density (float): Air density (kg/m3).
        S (float): Wing surface (m2)
        CL_max (float): Maximum lift coefficient (-)

    Returns:
        Vsr (float): Reference stall speed (m/s).

    """

    Vsr = ((m*g)/(0.5*density*S*CL_max))**0.5

    return Vsr

def get_critical_engine(adf_data):
    """Function to calculate the most critical engine in case of failure

    Args:
        adf_data (dict): Dictionary with the aircraft data for the selected condition.

    Returns:
       critical_eng (int): Index of the most critical engine.

    """

    Nt = [adf_data['TOGA'][i] * (- math.cos(np.radians(adf_data['epsilon'][i])) * math.cos(np.radians(adf_data['ni'][i])) * adf_data['CG_eng_y'][i] +
          math.cos(np.radians(adf_data['epsilon'][i])) * math.sin(np.radians(adf_data['ni'][i])) * adf_data['CG_eng_x'][i]) for i in range(adf_data['n_eng'])]
    
    critical_eng = Nt.index(min(Nt)) # Failure in the right side

    return critical_eng

def calculate_static_long(adf_data, input_data):
    """Function to calculate the static longitudinal stability and control

    Args:
        adf_data (dict): Dictionary with the aircraft data for the selected condition.
        input_data (dict): Dictionary with the input file variables.

    Returns:
        static_margin (float) : Static margin.
        neutral_point (float) : Neutral point.
        alpha_trim (float): Trim angle of attack (deg).
        delta_e_trim (float): Elevator deflection for trim (deg).
        delta_e (float): Elevator deflection for alpha (deg).

    """

    # Intermediate calculations   
    rho = atmos.density
    p_d = get_dynamic_pressure(rho, kt_to_ms * adf_data['KTAS'])
    
    """
    # Propulsion
    Ft_z = sum(-adf_data['T'][i] * math.sin(np.radians(adf_data['epsilon'][i]))  for i in range(adf_data['n_eng']))

    Mt_y = sum(adf_data['T'][i] * (math.cos(np.radians(adf_data['epsilon'][i])) * math.cos(np.radians(adf_data['ni'][i])) * adf_data['CG_eng_z'][i] +
               math.sin(np.radians(adf_data['epsilon'][i])) * adf_data['CG_eng_x'][i]) for i in range(adf_data['n_eng']))
    """
    
    # NEUTRAL POINT & STATIC MARGIN
    static_margin = adf_data['Cm_alpha']/adf_data['Cz_alpha']
    neutral_point = adf_data['CG'] + static_margin
    
    # TRIM
    CL_alpha   = - adf_data['Cz_alpha']
    CL_delta_e = - adf_data['Cz_delta_e']
    CL_0       = - adf_data['Cz_0']
    
    A = np.array([[            CL_alpha,             CL_delta_e],
                  [adf_data['Cm_alpha'], adf_data['Cm_delta_e']]])
    
    """
    B = np.array([[(adf_data['m'] * g + Ft_z) / (p_d * adf_data['S']) - CL_0],
                  [-Mt_y / (p_d * adf_data['S'] * adf_data['c']) -adf_data['Cm_0']]])
    """ 

    B = np.array([[(adf_data['m'] * g) / (p_d * adf_data['S']) - CL_0],
                  [-adf_data['Cm_0']]]) 
    
    X = np.linalg.solve(A,B)
    
    alpha_trim   = np.degrees(X[0])[0]
    delta_e_trim = np.degrees(X[1])[0]
    
    # CONTROL
    if 'alpha' in input_data:
        delta_e = np.array([])
        for alpha in input_data['alpha']:
            #d_e = ((- Mt_y) / (p_d * adf_data['S'] * adf_data['c']) - adf_data['Cm_0'] - adf_data['Cm_alpha'] * (np.radians(alpha))) / adf_data['Cm_delta_e']
            d_e = (- adf_data['Cm_0'] - adf_data['Cm_alpha'] * (np.radians(alpha))) / adf_data['Cm_delta_e']
    
            delta_e = np.append(delta_e, np.degrees(d_e))
    else:
        delta_e = np.nan
        
    return static_margin, neutral_point, alpha_trim, delta_e_trim, delta_e

def calculate_static_latdir(adf_data, input_data):
    """Function to calculate the static lateral-directionsl stability 

    Args:
        adf_data (dict): Dictionary with the aircraft data for the selected condition.
        input_data (dict): Dictionary with the input file variables.

    Returns:
        phi_b (float): Roll angle for trimming at given beta (deg).
        delta_a_b (float): Aileron deflection for trimming at given beta (deg).
        delta_r_b (float): Rudder deflection for trimming at given beta (deg).
        beta_p (float): Sidesleep angle for stationary turn at given phi (deg).
        delta_a_p (float): Aileron deflection for stationary turn at given phi (deg).
        delta_r_p (float): Rudder deflection for stationary turn at given phi (deg).
        n_p (float): Load factor in a stationary turn at given phi (deg).
    """


    # Intermediate calculations   
    rho = atmos.density
    p_d = get_dynamic_pressure(rho, kt_to_ms * adf_data['KTAS'])
    
    # Propulsion   
    Ft_y = sum(input_data['T_rate'][i] * adf_data['T'][i] * math.cos(np.radians(adf_data['epsilon'][i])) * math.sin(np.radians(adf_data['ni'][i])) for i in range(adf_data['n_eng']))
    
    Lt   = sum(-input_data['T_rate'][i] * adf_data['T'][i] * (math.cos(np.radians(adf_data['epsilon'][i])) * math.sin(np.radians(adf_data['ni'][i])) * adf_data['CG_eng_z'][i] +
             math.sin(np.radians(adf_data['epsilon'][i])) * adf_data['CG_eng_y'][i]) for i in range(adf_data['n_eng']))
    
    Nt   = sum(input_data['T_rate'][i] * adf_data['T'][i] * (- math.cos(np.radians(adf_data['epsilon'][i])) * math.cos(np.radians(adf_data['ni'][i])) * adf_data['CG_eng_y'][i] +
             math.cos(np.radians(adf_data['epsilon'][i])) * math.sin(np.radians(adf_data['ni'][i])) * adf_data['CG_eng_x'][i]) for i in range(adf_data['n_eng']))
    
    # TRIM 
    phi_b     = np.array([])
    delta_a_b = np.array([])
    delta_r_b = np.array([])
    
    beta_p    = np.array([])
    delta_a_p = np.array([])
    delta_r_p = np.array([])
    n_p       = np.array([])  # Load factor
    
    if 'beta' in input_data: # Fixed beta 
        A = np.array([[adf_data['m'] * g / (p_d * adf_data['S']), adf_data['Cy_delta_a'], adf_data['Cy_delta_r']],
                      [                                      0.0, adf_data['Cl_delta_a'], adf_data['Cl_delta_r']],
                      [                                      0.0, adf_data['Cn_delta_a'], adf_data['Cn_delta_r']]])

        for beta_ in input_data['beta']:          
            B = np.array([[-adf_data['Cy_beta'] * np.radians(beta_) - Ft_y / (p_d * adf_data['S'])],
                          [-adf_data['Cl_beta'] * np.radians(beta_) - Lt / (p_d * adf_data['S'] * adf_data['b'])],
                          [-adf_data['Cn_beta'] * np.radians(beta_) - Nt / (p_d * adf_data['S'] * adf_data['b'])]])
            
            X = np.linalg.solve(A,B)
        
            phi_b     = np.append(phi_b, np.degrees(X[0]))
            delta_a_b = np.append(delta_a_b, np.degrees(X[1]))
            delta_r_b = np.append(delta_r_b, np.degrees(X[2]))
        
    if 'phi' in input_data: # Stationary turn
        
        A = np.array([[adf_data['Cy_beta'], adf_data['Cy_delta_a'], adf_data['Cy_delta_r']],
                      [adf_data['Cl_beta'], adf_data['Cl_delta_a'], adf_data['Cl_delta_r']],
                      [adf_data['Cn_beta'], adf_data['Cn_delta_a'], adf_data['Cn_delta_r']]])
        
        u  = kt_to_ms * adf_data['KTAS']
        
        for phi_ in input_data['phi']:
            k1 = (adf_data['b'] * g * math.sin(np.radians(phi_))) / (2.0 * u ** 2)
            k2 = (g ** 2 * math.sin(np.radians(phi_)) ** 3) / (p_d * adf_data['S'] * adf_data['b'] * u ** 2 * math.cos(np.radians(phi_)))
            
            B = np.array([[-adf_data['Cy_r'] * k1 - Ft_y / (p_d * adf_data['S'])],
                          [(adf_data['Iz'] - adf_data['Iy']) * k2 - adf_data['Cl_r'] * k1 - Lt / (p_d * adf_data['S'] * adf_data['b'])],
                          [adf_data['Jxz'] * k2 - adf_data['Cn_r'] * k1 - Nt / (p_d * adf_data['S'] * adf_data['b'])]])
            
            X = np.linalg.solve(A,B)
            
            beta_p    = np.append(beta_p, np.degrees(X[0]))
            delta_a_p = np.append(delta_a_p, np.degrees(X[1]))
            delta_r_p = np.append(delta_r_p, np.degrees(X[2]))
            n_p       = np.append(n_p, 1.0 / math.cos(np.radians(phi_)))
        
    return phi_b, delta_a_b, delta_r_b, beta_p, delta_a_p, delta_r_p, n_p

def calculate_dynamic_long(adf_data, input_data):
    """Function to calculate the static lateral-directionsl stability 

    Args:
        adf_data (dict): Dictionary with the aircraft data for the selected condition.
        input_data (dict): Dictionary with the input file variables.

    Returns:
        lambda_ph (float): Eigenvalues in phugoid mode (-).
        wn_ph (float): Natural frequency in phugoid mode (1/s).
        zeta_ph (float): Dumping coefficient in phugoid mode (-).
        t_12_ph (float): Time constant in phugoid mode (s).
        T_ph (float): Period of the phugoid mode (s).
        lambda_sp (float): Natural frequency in short period mode (1/s).
        wn_sp (float): Natural frequency in short period mode (1/s).
        zeta_sp (float): Dumping coefficient in short period mode (-).
        t_12_sp (float): Time constant in short period mode (s).
        T_sp (float): Period of the short period mode (s).
        Control response charts
    """


    rho = atmos.density
    
    mu = adf_data['m'] / (0.5 * rho * adf_data['S'] *  adf_data['c'])
    
    Iy_ad = adf_data['Iy'] / (rho * adf_data['S'] * (adf_data['c'] / 2) ** 3)
    
    s= symbols('s')

    # STABILITY QUARTIC
    As = np.array([[               2*mu*s - adf_data['Cx_u'],                                      -adf_data['Cx_alpha'],               -adf_data['Cz_s']],
                   [-(adf_data['Cz_u'] + 2*adf_data['Cz_s']), (2*mu - adf_data['Cz_alpha_dot'])*s - adf_data['Cz_alpha'],    -(2*mu + adf_data['Cz_q'])*s],
                   [                       -adf_data['Cm_u'],       -(adf_data['Cm_alpha_dot']*s + adf_data['Cm_alpha']), Iy_ad*s**2 - adf_data['Cm_q']*s]])

    stability_quartic = Poly(As[0,0]*As[1,1]*As[2,2] + As[0,1]*As[1,2]*As[2,0] + As[0,2]*As[1,0]*As[2,1] - (As[0,2]*As[1,1]*As[2,0] + As[0,1]*As[1,0]*As[2,2] + As[0,0]*As[1,2]*As[2,1]),s)
    
    [A,B,C,D,E] = [float(stability_quartic.all_coeffs()[k]) for k in range(len(stability_quartic.all_coeffs()))]
    
    L_ad = np.roots([A, B, C, D, E])
    L = L_ad / (adf_data['c'] / (2 * kt_to_ms * adf_data['KTAS']))

    ind_ph = np.argmin(np.absolute(L.real))
    ind_sp = np.argmax(np.absolute(L.real)) 
    
    # PHUGOID
    lambda_ph = np.array([[L.real[ind_ph],  np.absolute(L.imag[ind_ph])],
                          [L.real[ind_ph], -np.absolute(L.imag[ind_ph])]])
    
    wn_ph   = (lambda_ph[0,0] ** 2 + lambda_ph[0,1] ** 2) ** 0.5
    zeta_ph = -lambda_ph[0,0] /((lambda_ph[0,0] ** 2 + lambda_ph[0,1] ** 2) ** 0.5)
    t_12_ph = -0.693/lambda_ph[0,0]
    T_ph = 2.0*math.pi / lambda_ph[0,1]
    
    #SHORT PERIOD
    lambda_sp = np.array([[L.real[ind_sp],  np.absolute(L.imag[ind_sp])],
                          [L.real[ind_sp], -np.absolute(L.imag[ind_sp])]])
    
    wn_sp   = (lambda_sp[0,0] ** 2 + lambda_sp[0,1] ** 2) ** 0.5
    zeta_sp = -lambda_sp[0,0] /((lambda_sp[0,0] ** 2 + lambda_sp[0,1] ** 2) ** 0.5)
    t_12_sp = -0.693/lambda_sp[0,0]
    T_sp = 2.0*math.pi / lambda_sp[0,1]
    
    """ 
    CONTROL
    
    """ 
    if input_data['plot'] == 1: 
        # Num
        Adj = []
        ind = np.arange(3)
        for i in range(3):
            for j in range(3):
                if (i+1 + j+1) % 2 == 0:
                    sign = 1
                else:
                    sign = -1
                    
                aux = As[ind != i]
                aux = aux[:, ind != j]
                pol = Poly(sign*(aux[0,0] * aux[1,1] - aux[0,1] * aux[1,0]), s)
                
                Adj.append(pol)
                
        Adj = np.reshape(Adj,[3,3])
        
        N_u_delta_e     = adf_data['Cx_delta_e'] * Adj[0,0] + adf_data['Cz_delta_e'] * Adj[1,0] + (adf_data['Cm_delta_e_dot']*s + adf_data['Cm_delta_e']) * Adj[2,0]
        N_alpha_delta_e = adf_data['Cx_delta_e'] * Adj[0,1] + adf_data['Cz_delta_e'] * Adj[1,1] + (adf_data['Cm_delta_e_dot']*s + adf_data['Cm_delta_e']) * Adj[2,1]
        N_theta_delta_e = adf_data['Cx_delta_e'] * Adj[0,2] + adf_data['Cz_delta_e'] * Adj[1,2] + (adf_data['Cm_delta_e_dot']*s + adf_data['Cm_delta_e']) * Adj[2,2]
        
        N_u_delta_e     = [float(N_u_delta_e.all_coeffs()[k]) for k in range(len(N_u_delta_e.all_coeffs()))]
        N_alpha_delta_e = [float(N_alpha_delta_e.all_coeffs()[k]) for k in range(len(N_alpha_delta_e.all_coeffs()))]
        N_theta_delta_e = [float(N_theta_delta_e.all_coeffs()[k]) for k in range(len(N_theta_delta_e.all_coeffs()))]
        
        # Den
        Ds = [A, B, C, D, E]
        
        # Transfer function
        G_u_delta_e     = (N_u_delta_e, Ds)
        G_alpha_delta_e = (N_alpha_delta_e, Ds)
        G_theta_delta_e = (N_theta_delta_e, Ds)
        
        #CHARTS
        n_points = 500
        height_in = 11.69
        width_in  = 8.27
        font_size = 15
            
        tf = [ G_u_delta_e,    G_alpha_delta_e,    G_theta_delta_e]
        yl = ['$\\Delta$û', '$\\Delta\\alpha$', '$\\Delta\\theta$']
        xl = '$\\^t$'
                        
        # Response to impulse
        tit = 'Response to impulse\n $\\delta$e'
        fig_name = [out_data['file'][:-4] + '_dynamic_response_delta_e_impulse_u' + '.png',
                    out_data['file'][:-4] + '_dynamic_response_delta_e_impulse_alpha' + '.png',
                    out_data['file'][:-4] + '_dynamic_response_delta_e_impulse_theta' + '.png']
        for i in range(len(tf)):
            t, y = signal.impulse(tf[i], N=n_points)
            fig, ax1 = plt.subplots(figsize=(height_in, width_in))  
            ax1.plot(t, y) 
            ax1.set_xlabel(xl, fontsize=font_size)
            ax1.set_ylabel(yl[i], fontsize=font_size)
            ax1.set_title(tit, fontsize=font_size)
            ax1.grid(b=True, which='major', color='#666666', linestyle='-')
            ax1.minorticks_on()
            ax1.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)
            ax1.tick_params(labelsize=font_size)
            fig.tight_layout()
            fig.savefig(fig_name[i])
            #plt.show()
            
        # Response to step
        tit = 'Response to step\n $\\delta$e'
        fig_name = [out_data['file'][:-4] + '_dynamic_response_delta_e_step_u' + '.png',
                    out_data['file'][:-4] + '_dynamic_response_delta_e_step_alpha' + '.png',
                    out_data['file'][:-4] + '_dynamic_response_delta_e_step_theta' + '.png']
        for i in range(len(tf)):
            t, y = signal.step(tf[i], N=n_points)
            fig, ax1 = plt.subplots(figsize=(height_in, width_in))  
            ax1.plot(t, y) 
            ax1.set_xlabel(xl, fontsize=font_size)
            ax1.set_ylabel(yl[i], fontsize=font_size)
            ax1.set_title(tit, fontsize=font_size)
            ax1.grid(b=True, which='major', color='#666666', linestyle='-')
            ax1.minorticks_on()
            ax1.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)
            ax1.tick_params(labelsize=font_size)
            fig.tight_layout()
            fig.savefig(fig_name[i])
            #plt.show()
    
    return lambda_ph, wn_ph, zeta_ph, t_12_ph, T_ph, lambda_sp, wn_sp, zeta_sp, t_12_sp, T_sp

def calculate_dynamic_latdir(adf_data, input_data):
    """Function to calculate the dynamic lateral-directional stability

    Args:
        adf_data (dict): Dictionary with the aircraft data for the selected condition.
        input_data (dict): Dictionary with the input file variables.

    Returns:
        lambda_rs (float): Eigenvalues in roll subsidence mode (-).
        t_12_rs (float): Time constant in roll subsidence mode (s).
        lambda_spi (float): Eigenvalues in spiral mode (-).
        t_12_spi (float): Time constant in spiral mode (s).
        lambda_dr (float): Eigenvalues in dutch roll mode (-).
        wn_dr (float): Natural frequency in dutch roll mode (1/s).
        zeta_dr (float): Dumping coefficient in dutch roll mode (-).
        t_12_dr (float): Time constant in dutch roll mode (s).
        T_dr (float): Period of the dutch roll mode (s).
        Control response charts
    """
    
    rho = atmos.density
    
    mu = adf_data['m'] / (0.5 * rho * adf_data['S'] *  adf_data['b'])
    
    Ix_ad  = adf_data['Ix'] / (rho * adf_data['S'] * (adf_data['b'] / 2) ** 3)
    Iz_ad  = adf_data['Iz'] / (rho * adf_data['S'] * (adf_data['b'] / 2) ** 3)
    Jxz_ad = adf_data['Jxz'] / (rho * adf_data['S'] * (adf_data['b'] / 2) ** 3)
    
    s= symbols('s')

    # STABILITY QUARTIC
    As = np.array([[ 2*mu*s - adf_data['Cy_beta'], -(adf_data['Cy_p']*s - adf_data['Cz_s']),        2*mu - adf_data['Cy_r']],
                   [         -adf_data['Cl_beta'],          Ix_ad*s**2 - adf_data['Cl_p']*s, -(Jxz_ad*s + adf_data['Cl_r'])],
                   [         -adf_data['Cn_beta'],      -(Jxz_ad*s**2 + adf_data['Cn_p']*s),     Iz_ad*s - adf_data['Cn_r']]])
    
    stability_quartic = Poly(As[0,0]*As[1,1]*As[2,2] + As[0,1]*As[1,2]*As[2,0] + As[0,2]*As[1,0]*As[2,1] - (As[0,2]*As[1,1]*As[2,0] + As[0,1]*As[1,0]*As[2,2] + As[0,0]*As[1,2]*As[2,1]),s)
    
    [A,B,C,D,E] = [float(stability_quartic.all_coeffs()[k]) for k in range(len(stability_quartic.all_coeffs()))]
    
    L_ad = np.roots([A, B, C, D, E])
    L = L_ad / (adf_data['b'] / (2 * kt_to_ms * adf_data['KTAS']))
      
    ind_dr  = np.argmax(np.absolute(L.imag))
    
    L_aux = L.real[L.imag== 0.0]
    ind_rs  = np.argmax(np.absolute(L_aux))
    ind_spi = np.argmin(np.absolute(L_aux))

    # ROLL SUBSIDENCE
    lambda_rs = np.array([L_aux[ind_rs],  0.0])

    t_12_rs = -0.693/lambda_rs[0]
    
    # SPIRAL
    lambda_spi = np.array([L_aux[ind_spi],  0.0])

    t_12_spi = -0.693/lambda_spi[0]
    
    # DUTCH ROLL
    lambda_dr = np.array([[L.real[ind_dr],  np.absolute(L.imag[ind_dr])],
                          [L.real[ind_dr], -np.absolute(L.imag[ind_dr])]])
    
    wn_dr   = (lambda_dr[0,0] ** 2 + lambda_dr[0,1] ** 2) ** 0.5
    zeta_dr = -lambda_dr[0,0] /((lambda_dr[0,0] ** 2 + lambda_dr[0,1] ** 2) ** 0.5)
    t_12_dr = -0.693/lambda_dr[0,0]
    T_dr = 2.0*math.pi / lambda_dr[0,1]
    
    """ 
    CONTROL
    
    """
    if input_data['plot'] == 1:
        # Num
        Adj = []
        ind = np.arange(3)
        for i in range(3):
            for j in range(3):
                if (i+1 + j+1) % 2 == 0:
                    sign = 1
                else:
                    sign = -1
                    
                aux = As[ind != i]
                aux = aux[:, ind != j]
                pol = Poly(sign*(aux[0,0] * aux[1,1] - aux[0,1] * aux[1,0]), s)
                
                Adj.append(pol)
                
        Adj = np.reshape(Adj,[3,3])
        
        N_beta_delta_a = (adf_data['Cl_delta_a_dot']*s + adf_data['Cl_delta_a']) * Adj[1,0] + adf_data['Cn_delta_a'] * Adj[2,0]
        N_phi_delta_a  = (adf_data['Cl_delta_a_dot']*s + adf_data['Cl_delta_a']) * Adj[1,1] + adf_data['Cn_delta_a'] * Adj[2,1]
        N_r_delta_a    = (adf_data['Cl_delta_a_dot']*s + adf_data['Cl_delta_a']) * Adj[1,2] + adf_data['Cn_delta_a'] * Adj[2,2]
        
        N_beta_delta_a = [float(N_beta_delta_a.all_coeffs()[k]) for k in range(len(N_beta_delta_a.all_coeffs()))]
        N_phi_delta_a  = [float(N_phi_delta_a.all_coeffs()[k]) for k in range(len(N_phi_delta_a.all_coeffs()))]
        N_r_delta_a    = [float(N_r_delta_a.all_coeffs()[k]) for k in range(len(N_r_delta_a.all_coeffs()))]
        
        
        N_beta_delta_r = adf_data['Cy_delta_r'] * Adj[0,0] + adf_data['Cl_delta_r'] * Adj[1,0] + (adf_data['Cn_delta_r_dot']*s + adf_data['Cn_delta_r']) * Adj[2,0]
        N_phi_delta_r  = adf_data['Cy_delta_r'] * Adj[0,1] + adf_data['Cl_delta_r'] * Adj[1,1] + (adf_data['Cn_delta_r_dot']*s + adf_data['Cn_delta_r']) * Adj[2,1]
        N_r_delta_r    = adf_data['Cy_delta_r'] * Adj[0,2] + adf_data['Cl_delta_r'] * Adj[1,2] + (adf_data['Cn_delta_r_dot']*s + adf_data['Cn_delta_r']) * Adj[2,2]
        
        N_beta_delta_r = [float(N_beta_delta_r.all_coeffs()[k]) for k in range(len(N_beta_delta_r.all_coeffs()))]
        N_phi_delta_r  = [float(N_phi_delta_r.all_coeffs()[k]) for k in range(len(N_phi_delta_r.all_coeffs()))]
        N_r_delta_r    = [float(N_r_delta_r.all_coeffs()[k]) for k in range(len(N_r_delta_r.all_coeffs()))]
        
        # Den
        Ds = [A, B, C, D, E]
        
        # Transfer function
        G_beta_delta_a = (N_beta_delta_a, Ds)
        G_phi_delta_a  = (N_phi_delta_a, Ds)
        G_r_delta_a    = (N_r_delta_a, Ds)
        
        G_beta_delta_r = (N_beta_delta_r, Ds)
        G_phi_delta_r  = (N_phi_delta_r, Ds)
        G_r_delta_r    = (N_r_delta_r, Ds)
        
        # CHARTS
        n_points = 500
        height_in = 11.69
        width_in  = 8.27
        font_size = 15
            
        tf_a = [   G_beta_delta_a,    G_phi_delta_a,    G_r_delta_a]
        tf_r = [   G_beta_delta_r,    G_phi_delta_r,    G_r_delta_r]
        yl   = ['$\\Delta\\beta$', '$\\Delta\\Phi$', '$\\Delta\\^r$']
        xl   = '$\\^t$'
                    
        # Response to impulse
        tit = 'Response to impulse\n $\\delta$a'
        fig_name = [out_data['file'][:-4] + '_dynamic_response_delta_a_impulse_beta' + '.png',
                    out_data['file'][:-4] + '_dynamic_response_delta_a_impulse_phi' + '.png',
                    out_data['file'][:-4] + '_dynamic_response_delta_a_impulse_r' + '.png']
        for i in range(len(tf_a)):
            t, y = signal.impulse(tf_a[i], N=n_points)
            fig, ax1 = plt.subplots(figsize=(height_in, width_in))  
            ax1.plot(t, y) 
            ax1.set_xlabel(xl, fontsize=font_size)
            ax1.set_ylabel(yl[i], fontsize=font_size)
            ax1.set_title(tit, fontsize=font_size)
            ax1.grid(b=True, which='major', color='#666666', linestyle='-')
            ax1.minorticks_on()
            ax1.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)
            ax1.tick_params(labelsize=font_size)
            fig.tight_layout()
            fig.savefig(fig_name[i])
            #plt.show()
            
        tit = 'Response to impulse\n $\\delta$r'
        fig_name = [out_data['file'][:-4] + '_dynamic_response_delta_r_impulse_beta' + '.png',
                    out_data['file'][:-4] + '_dynamic_response_delta_r_impulse_phi' + '.png',
                    out_data['file'][:-4] + '_dynamic_response_delta_r_impulse_r' + '.png']    
        for i in range(len(tf_r)):
            t, y = signal.impulse(tf_r[i], N=n_points)
            fig, ax1 = plt.subplots(figsize=(height_in, width_in))  
            ax1.plot(t, y) 
            ax1.set_xlabel(xl, fontsize=font_size)
            ax1.set_ylabel(yl[i], fontsize=font_size)
            ax1.set_title(tit, fontsize=font_size)
            ax1.grid(b=True, which='major', color='#666666', linestyle='-')
            ax1.minorticks_on()
            ax1.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)
            ax1.tick_params(labelsize=font_size)
            fig.tight_layout()
            fig.savefig(fig_name[i])
            #plt.show()
            
        # Response to step
        tit = 'Response to step\n $\\delta$a'
        fig_name = [out_data['file'][:-4] + '_dynamic_response_delta_a_step_beta' + '.png',
                    out_data['file'][:-4] + '_dynamic_response_delta_a_step_phi' + '.png',
                    out_data['file'][:-4] + '_dynamic_response_delta_a_step_r' + '.png']
        for i in range(len(tf_a)):
            t, y = signal.step(tf_a[i], N=n_points)
            fig, ax1 = plt.subplots(figsize=(height_in, width_in))  
            ax1.plot(t, y) 
            ax1.set_xlabel(xl, fontsize=font_size)
            ax1.set_ylabel(yl[i], fontsize=font_size)
            ax1.set_title(tit, fontsize=font_size)
            ax1.grid(b=True, which='major', color='#666666', linestyle='-')
            ax1.minorticks_on()
            ax1.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)
            ax1.tick_params(labelsize=font_size)
            fig.tight_layout()
            fig.savefig(fig_name[i])
            #plt.show()
            
        tit = 'Response to step\n $\\delta$r'
        fig_name = [out_data['file'][:-4] + '_dynamic_response_delta_r_step_beta' + '.png',
                    out_data['file'][:-4] + '_dynamic_response_delta_r_step_phi' + '.png',
                    out_data['file'][:-4] + '_dynamic_response_delta_r_step_r' + '.png']    
        for i in range(len(tf_r)):
            t, y = signal.step(tf_r[i], N=n_points)
            fig, ax1 = plt.subplots(figsize=(height_in, width_in))  
            ax1.plot(t, y) 
            ax1.set_xlabel(xl, fontsize=font_size)
            ax1.set_ylabel(yl[i], fontsize=font_size)
            ax1.set_title(tit, fontsize=font_size)
            ax1.grid(b=True, which='major', color='#666666', linestyle='-')
            ax1.minorticks_on()
            ax1.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)
            ax1.tick_params(labelsize=font_size)
            fig.tight_layout()
            fig.savefig(fig_name[i])
            #plt.show()

    return lambda_rs, t_12_rs, lambda_spi, t_12_spi, lambda_dr, wn_dr, zeta_dr, t_12_dr, T_dr

def calculate_vmcg(adf_data):
    """Function to calculate the minimum control speed on the ground (VMCG)
    For the Take-off thrust at SL and ISA conditions

    Args:
        adf_data (dict): Dictionary with the aircraft data for the selected condition.

    Returns:
        VMCG (float): Minimum control speed on the groung (kt).
    """
    # Check that take-off condition exists
    if adf_data['TO_condition'] == 0:
        return

    # Variable definition
    altitude = 0.0
    ISA = 0.0
    rho = Atmos.get_density(altitude, ISA)
    n_lg = len(adf_data['CG_lg_x'])
    Vsr = get_stall_speed(adf_data['TOW_min'], rho, adf_data['S'], adf_data['CL_max_TO'])

    # Control limit
    delta_r_max = -np.radians(adf_data['delta_r_max']) # Negative due to right engine failure

    # Propulsion
    critical_eng = get_critical_engine(adf_data)
    operative_eng = [eng for eng in range(adf_data['n_eng']) if eng != critical_eng]
    
    Ft_y = sum(adf_data['TOGA'][i] * math.cos(np.radians(adf_data['epsilon'][i])) * math.sin(np.radians(adf_data['ni'][i])) for i in operative_eng)
    Ft_z = sum(-adf_data['TOGA'][i] * math.sin(np.radians(adf_data['epsilon'][i]))  for i in operative_eng)

    Lt = sum(adf_data['TOGA'][i] * (math.cos(np.radians(adf_data['epsilon'][i])) * math.sin(np.radians(adf_data['ni'][i])) * adf_data['CG_eng_z'][i] +
             math.sin(np.radians(adf_data['epsilon'][i])) * adf_data['CG_eng_y'][i]) for i in operative_eng)
    Mt = sum(adf_data['TOGA'][i] * (math.cos(np.radians(adf_data['epsilon'][i])) * math.cos(np.radians(adf_data['ni'][i])) * adf_data['CG_eng_z'][i] +
             math.sin(np.radians(adf_data['epsilon'][i])) * adf_data['CG_eng_x'][i]) for i in operative_eng)
    Nt = sum(adf_data['TOGA'][i] * (- math.cos(np.radians(adf_data['epsilon'][i])) * math.cos(np.radians(adf_data['ni'][i])) * adf_data['CG_eng_y'][i] +
             math.cos(np.radians(adf_data['epsilon'][i])) * math.sin(np.radians(adf_data['ni'][i])) * adf_data['CG_eng_x'][i]) for i in operative_eng)
    
    # First aproximation (without LG effect)
    def equations_aprox(var):
        V, beta = var

        eq2 = Ft_y + 0.5*rho*V**2*adf_data['S']*(adf_data['Cy_beta_TO']*beta + adf_data['Cy_delta_r_TO']*delta_r_max)
        eq6 = Nt + 0.5*rho*V**2*adf_data['S']*adf_data['b']*(adf_data['Cn_beta_TO']*beta + adf_data['Cn_delta_r_TO']*delta_r_max)

        return [eq2, eq6]

    ini_aprox = [Vsr, 0.]
    VMCG_ini, beta_ini = fsolve(equations_aprox, ini_aprox)

    # Final system (with LG effect)
    # x = [VMCG, beta, K_lg, W]
    def equations(var):
        x = var

        eq2 = Ft_y + 0.5*rho*x[0]**2*adf_data['S']*(adf_data['Cy_beta_TO']*x[1] + adf_data['Cy_delta_r_TO']*delta_r_max) + x[2]*x[1]*sum(x[3+i] for i in range(n_lg))
        eq3 = Ft_z + adf_data['TOW_min']*g + 0.5*rho*x[0]**2*adf_data['S']*adf_data['Cz_0_TO'] + sum(x[3+i] for i in range(n_lg))
        eq4 = Lt + 0.5*rho*x[0]**2*adf_data['S']*adf_data['b']*(adf_data['Cl_beta_TO']*x[1] + adf_data['Cl_delta_r_TO']*delta_r_max) + sum(x[3+i]*adf_data['CG_lg_y'][i] for i in range(n_lg)) - x[2]*x[1]*sum(x[3+i]*adf_data['CG_lg_z'][i] for i in range(n_lg))
        eq5 = Mt + 0.5*rho*x[0]**2*adf_data['S']*adf_data['c']*adf_data['Cm_0_TO'] - sum(x[3+i]*adf_data['CG_lg_x'][i] for i in range(n_lg)) + adf_data['mu_gr']*sum(x[3+i]*adf_data['CG_lg_z'][i] for i in range(n_lg))
        eq6 = Nt + 0.5*rho*x[0]**2*adf_data['S']*adf_data['b']*(adf_data['Cn_beta_TO']*x[1] + adf_data['Cn_delta_r_TO']*delta_r_max) - adf_data['mu_gr']*sum(x[3+i]*adf_data['CG_lg_y'][i] for i in range(n_lg)) + x[2]*x[1]*sum(x[3+i]*adf_data['CG_lg_x'][i] for i in range(n_lg))
        eq7 = -x[2] + 2369.074348*abs(x[1])**5 - 1268.0219021*abs(x[1])**4 + 24.6344543*abs(x[1])**3 + 3.1599044*abs(x[1])**2 + 5.2489678*abs(x[1]) + 0.0086458

        return [eq2, eq3, eq4, eq5, eq6, eq7]

    ini    = np.empty(3 + n_lg)
    ini[:] = np.nan

    ini[0]  = VMCG_ini
    ini[1]  = beta_ini
    ini[2]  = 2369.074348*abs(beta_ini)**5 - 1268.0219021*abs(beta_ini)**4 + 24.6344543*abs(beta_ini)**3 + 3.1599044*abs(beta_ini)**2 + 5.2489678*abs(beta_ini) + 0.0086458
    ini[3:] = [-(adf_data['TOW_min']*g+0.5*rho*VMCG_ini**2*adf_data['S']*adf_data['Cz_0_TO'])/n_lg for i in range(n_lg)]
    
    x = fsolve(equations, ini)
    
    VMCG = Speed.tas_to_cas(x[0], altitude, ISA)/kt_to_ms
    beta = np.degrees(x[1])
    K_lg = x[2]
    W = x[3:]

    return VMCG

def calculate_vmca(adf_data, input_data):
    """Function to calculate the minimum control speed in the air (VMCA)
    For the Take-off thrust at SL and ISA conditions

    Args:
        adf_data (dict): Dictionary with the aircraft data for the selected condition.
        input_data (dict): Dictionary with the input file variables.

    Returns:
        VMCA (float): Minimum control speed in the air for each weight (kg / kt).
        VMCA_m (float): Minimum control speed in the air for the current weight (kt).
    """

    # Check that take-off condition exists
    if adf_data['TO_condition'] == 0:
        return ['None']*2

    # Regulation conditions
    fixed_phi = -5.0 # Negative due to right engine failure

    # Variable definition
    altitude = 0.0
    ISA = 0.0
    rho = Atmos.get_density(altitude, ISA)
    di  = 2.0
    Cd_eng  = (0.1934*(di**2.0))/adf_data['S']
    
    beta          = []
    delta_a       = []
    delta_r       = []

    weight_vector = []
    VMCA          = []
    VMCA_m        = []
    
    # Control limit
    delta_a_max = -np.radians(adf_data['delta_a_max']) # Negative due to right engine failure
    delta_r_max = -np.radians(adf_data['delta_r_max']) # Negative due to right engine failure
    
    # Propulsion
    critical_eng = get_critical_engine(adf_data)
    operative_eng = [eng for eng in range(adf_data['n_eng']) if eng != critical_eng]
    
    Ft_y = sum(adf_data['TOGA'][i] * math.cos(np.radians(adf_data['epsilon'][i])) * math.sin(np.radians(adf_data['ni'][i])) for i in operative_eng)

    Lt   = sum(adf_data['TOGA'][i] * (math.cos(np.radians(adf_data['epsilon'][i])) * math.sin(np.radians(adf_data['ni'][i])) * adf_data['CG_eng_z'][i] +
             math.sin(np.radians(adf_data['epsilon'][i])) * adf_data['CG_eng_y'][i]) for i in operative_eng)

    Nt   = sum(adf_data['TOGA'][i] * (- math.cos(np.radians(adf_data['epsilon'][i])) * math.cos(np.radians(adf_data['ni'][i])) * adf_data['CG_eng_y'][i] +
             math.cos(np.radians(adf_data['epsilon'][i])) * math.sin(np.radians(adf_data['ni'][i])) * adf_data['CG_eng_x'][i]) for i in operative_eng)

    # Beta, VMCA_cross and W_cross for delta_a_max and delta_r_max   (1)
    A = np.array([[adf_data['Cl_beta_TO'], Lt / (0.5 * rho * adf_data['S'] * adf_data['b'])],
                  [adf_data['Cn_beta_TO'], Nt / (0.5 * rho * adf_data['S'] * adf_data['b'])]])

    b = np.array([[-adf_data['Cl_delta_r_TO'] * delta_r_max - adf_data['Cl_delta_a_TO'] * delta_a_max],
                  [-adf_data['Cn_delta_r_TO'] * delta_r_max - adf_data['Cn_delta_a_TO'] * delta_a_max - (Cd_eng * adf_data['CG_eng_y'][critical_eng] / adf_data['b'])]])

    X = np.linalg.solve(A,b)

    betacross     = np.degrees(X[0])
    VMCAcross     = ((X[1])**(-0.5))
    Cy            = (adf_data['Cy_beta_TO'] * X[0] + adf_data['Cy_delta_a_TO'] * delta_a_max + adf_data['Cy_delta_r_TO'] * delta_r_max)
    Wcross        = -(Ft_y + 0.5 * rho * X[1]**(-1) * adf_data['S'] * Cy) / (g * math.sin(np.radians(fixed_phi)))

    if adf_data['TOW_min'] <= Wcross <= adf_data['TOW_max']:
        beta          = np.append(beta, betacross)
        delta_a       = np.append(delta_a, np.degrees(delta_a_max))
        delta_r       = np.append(delta_r, np.degrees(delta_r_max))
        VMCA          = np.append(VMCA, VMCAcross)
        weight_vector = np.append(weight_vector, Wcross)    

    # VMCA for a given weight delta_r_max or delta_a_max
    for W in np.arange(adf_data['TOW_min'], adf_data['TOW_max'] + 1, adf_data['dW']):
        if W < Wcross: # delta_r_max
            A = np.array([[adf_data['Cy_beta_TO'], adf_data['Cy_delta_a_TO'], (W * g * math.sin(np.radians(fixed_phi)) + Ft_y) / (0.5 * rho * adf_data['S'])],
                          [adf_data['Cl_beta_TO'], adf_data['Cl_delta_a_TO'], Lt / (0.5 * rho * adf_data['S'] * adf_data['b'])],
                          [adf_data['Cn_beta_TO'], adf_data['Cn_delta_a_TO'], Nt / (0.5 * rho * adf_data['S'] * adf_data['b'])]])

            b = np.array([[-adf_data['Cy_delta_r_TO'] * delta_r_max],
                          [-adf_data['Cl_delta_r_TO'] * delta_r_max],
                          [-adf_data['Cn_delta_r_TO'] * delta_r_max - (Cd_eng * adf_data['CG_eng_y'][critical_eng] / adf_data['b'])]])

            X = np.linalg.solve(A,b)
            
            weight_vector = np.append(weight_vector, W)
            beta          = np.append(beta, np.degrees(X[0]))
            delta_a       = np.append(delta_a, np.degrees(X[1]))
            delta_r       = np.append(delta_r, np.degrees(delta_r_max))
            VMCA          = np.append(VMCA, ((X[2])**(-0.5)))

        if W > Wcross: # delta_a_max
            A = np.array([[adf_data['Cy_beta_TO'], adf_data['Cy_delta_r_TO'], ((W * g * math.sin(np.radians(fixed_phi)) + Ft_y) / (0.5 * rho * adf_data['S']))],
                          [adf_data['Cl_beta_TO'], adf_data['Cl_delta_r_TO'], Lt / (0.5 * rho * adf_data['S'] * adf_data['b'])],
                          [adf_data['Cn_beta_TO'], adf_data['Cn_delta_r_TO'], Nt / (0.5 * rho * adf_data['S'] * adf_data['b'])]])

            b = np.array([[-adf_data['Cy_delta_a_TO'] * delta_a_max],
                          [-adf_data['Cl_delta_a_TO'] * delta_a_max],
                          [-adf_data['Cn_delta_a_TO'] * delta_a_max - (Cd_eng * adf_data['CG_eng_y'][critical_eng] / adf_data['b'])]])

            X = np.linalg.solve(A,b)

            weight_vector = np.append(weight_vector, W)
            beta          = np.append(beta, np.degrees(X[0]))
            delta_a       = np.append(delta_a, np.degrees(delta_a_max))
            delta_r       = np.append(delta_r, np.degrees(X[1]))
            VMCA          = np.append(VMCA, ((X[2])**(-0.5)))
    
    # Sort vectors by weight
    ind           = np.argsort(weight_vector)
    weight_vector = weight_vector[ind]
    VMCA          = VMCA[ind]
    beta          = beta[ind]
    delta_a       = delta_a[ind]
    delta_r       = delta_r[ind]

    # Convert VMCA from TAS to KCAS
    for i in range(len(VMCA)):
        VMCA[i] = Speed.tas_to_cas(VMCA[i], altitude, ISA)/kt_to_ms 

    # VMCA for given initial weight (m) (4)
    f       = interp1d(weight_vector, VMCA)
    VMCA_m  = f(adf_data['m'])

    # Calculate VSR
    Vsr  = [Speed.tas_to_cas(get_stall_speed(W, rho, adf_data['S'], adf_data['CL_max_TO']), altitude, ISA)/kt_to_ms for W in weight_vector]
    KVsr = np.multiply(Vsr, 1.13)

    # Plot VMCA and 1.13Vsr for each weight
    if input_data['plot'] == 1: 
        height_in = 11.69
        width_in  = 8.27
        font_size = 15  

        fig_name = out_data['file'][:-4] + '_VMCA.png'
        xl  = 'Weight (kg)'
        yl  = 'Speed (kt)'
        tit = '$V_{MCA}$'

        fig, ax1 = plt.subplots(figsize=(height_in, width_in))
        ax1.plot(weight_vector, VMCA, color = 'tab:blue', label = '$V_{MCA}$')
        ax1.plot(weight_vector, KVsr, color = 'tab:red', label = '$1.13V_{SR}$')

        ax1.set_xlabel(xl, fontsize=font_size)
        ax1.set_ylabel(yl, fontsize=font_size)
        ax1.set_title(tit, fontsize=font_size)
        
        ax1.grid(b=True, which='major', color='#666666', linestyle='-')
        ax1.minorticks_on()
        ax1.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)
        ax1.tick_params(labelsize=font_size)
        ax1.legend(loc='best')

        fig.tight_layout()
        fig.savefig(fig_name)
    
    # Merge weight, VMCA and 1.13VSR in one variable
    VMCA = [weight_vector, VMCA, KVsr]
   
    return VMCA, VMCA_m

def calculate_vmcl(adf_data, input_data):
    """Function to calculate the minimum control speed during aproach and langing (VMCL)
    For the Go-around thrust at SL and ISA conditions

    Args:
        adf_data (dict): Dictionary with the aircraft data for the selected condition.
        input_data (dict): Dictionary with the input file variables.

    Returns:
        VMCL (float): Minimum control speed during aproach and langing for each weight (kg / kt).
        VMCL_m (float): Minimum control speed during aproach and langing for the current weight (kt).
    """

    # Check that take-off condition exists
    if adf_data['LD_condition'] == 0:
        return ['None']*2

    # Regulation conditions
    fixed_phi = -5.0 # Negative due to right engine failure

    # Variable definition
    altitude = 0.0
    ISA = 0.0
    rho = Atmos.get_density(altitude, ISA)
    di  = 2.0
    Cd_eng  = (0.1934*(di**2.0))/adf_data['S']

    beta          = []
    delta_a       = []
    delta_r       = []

    weight_vector = []
    VMCL          = []
    VMCL_m        = []

    # Control limit
    delta_a_max = -np.radians(adf_data['delta_a_max']) # Negative due to right engine failure
    delta_r_max = -np.radians(adf_data['delta_r_max']) # Negative due to right engine failure

    # Propulsion
    critical_eng = get_critical_engine(adf_data)
    operative_eng = [eng for eng in range(adf_data['n_eng']) if eng != critical_eng]
    
    Ft_y = sum(adf_data['TOGA'][i] * math.cos(np.radians(adf_data['epsilon'][i])) * math.sin(np.radians(adf_data['ni'][i])) for i in operative_eng)

    Lt   = sum(adf_data['TOGA'][i] * (math.cos(np.radians(adf_data['epsilon'][i])) * math.sin(np.radians(adf_data['ni'][i])) * adf_data['CG_eng_z'][i] +
             math.sin(np.radians(adf_data['epsilon'][i])) * adf_data['CG_eng_y'][i]) for i in operative_eng)

    Nt   = sum(adf_data['TOGA'][i] * (- math.cos(np.radians(adf_data['epsilon'][i])) * math.cos(np.radians(adf_data['ni'][i])) * adf_data['CG_eng_y'][i] +
             math.cos(np.radians(adf_data['epsilon'][i])) * math.sin(np.radians(adf_data['ni'][i])) * adf_data['CG_eng_x'][i]) for i in operative_eng)

    # Beta, VMCL_cross and W_cross for delta_a_max and delta_r_max   (1)
    A = np.array([[adf_data['Cl_beta_LD'], Lt / (0.5 * rho * adf_data['S'] * adf_data['b'])],
                  [adf_data['Cn_beta_LD'], Nt / (0.5 * rho * adf_data['S'] * adf_data['b'])]])

    b = np.array([[-adf_data['Cl_delta_r_LD'] * delta_r_max - adf_data['Cl_delta_a_LD'] * delta_a_max],
                  [-adf_data['Cn_delta_r_LD'] * delta_r_max - adf_data['Cn_delta_a_LD'] * delta_a_max - (Cd_eng * adf_data['CG_eng_y'][critical_eng] / adf_data['b'])]])

    X = np.linalg.solve(A,b)

    betacross     = np.degrees(X[0])
    VMCLcross     = ((X[1])**(-0.5))
    Cy            = (adf_data['Cy_beta_LD'] * X[0] + adf_data['Cy_delta_a_LD'] * delta_a_max + adf_data['Cy_delta_r_LD'] * delta_r_max)
    Wcross        = -(Ft_y + 0.5 * rho * X[1]**(-1) * adf_data['S'] * Cy) / (g * math.sin(np.radians(fixed_phi)))

    if adf_data['LW_min'] <= Wcross <= adf_data['LW_max']:
        weight_vector = np.append(weight_vector, Wcross)
        VMCL          = np.append(VMCL, VMCLcross)
        beta          = np.append(beta, betacross)
        delta_a       = np.append(delta_a, np.degrees(delta_a_max))
        delta_r       = np.append(delta_r, np.degrees(delta_r_max))

    # VMCL for a given weight delta_r_max or delta_a_max
    for W in np.arange(adf_data['LW_min'], adf_data['LW_max'] + 1, adf_data['dW']):
        if W < Wcross: # delta_r_max
            A = np.array([[adf_data['Cy_beta_LD'], adf_data['Cy_delta_a_LD'], (W * g * math.sin(np.radians(fixed_phi)) + Ft_y) / (0.5 * rho * adf_data['S'])],
                          [adf_data['Cl_beta_LD'], adf_data['Cl_delta_a_LD'], Lt / (0.5 * rho * adf_data['S'] * adf_data['b'])],
                          [adf_data['Cn_beta_LD'], adf_data['Cn_delta_a_LD'], Nt / (0.5 * rho * adf_data['S'] * adf_data['b'])]])

            b = np.array([[-adf_data['Cy_delta_r_LD'] * delta_r_max],
                          [-adf_data['Cl_delta_r_LD'] * delta_r_max],
                          [-adf_data['Cn_delta_r_LD'] * delta_r_max - (Cd_eng * adf_data['CG_eng_y'][critical_eng] / adf_data['b'])]])

            X = np.linalg.solve(A,b)
            
            weight_vector = np.append(weight_vector, W)
            VMCL          = np.append(VMCL, ((X[2])**(-0.5)))
            beta          = np.append(beta, np.degrees(X[0]))
            delta_a       = np.append(delta_a, np.degrees(X[1]))
            delta_r       = np.append(delta_r, np.degrees(delta_r_max))

        if W > Wcross: # delta_a_max
            A = np.array([[adf_data['Cy_beta_LD'], adf_data['Cy_delta_r_LD'], ((W * g * math.sin(np.radians(fixed_phi)) + Ft_y) / (0.5 * rho * adf_data['S']))],
                          [adf_data['Cl_beta_LD'], adf_data['Cl_delta_r_LD'], Lt / (0.5 * rho * adf_data['S'] * adf_data['b'])],
                          [adf_data['Cn_beta_LD'], adf_data['Cn_delta_r_LD'], Nt / (0.5 * rho * adf_data['S'] * adf_data['b'])]])

            b = np.array([[-adf_data['Cy_delta_a_LD'] * delta_a_max],
                          [-adf_data['Cl_delta_a_LD'] * delta_a_max],
                          [-adf_data['Cn_delta_a_LD'] * delta_a_max - (Cd_eng * adf_data['CG_eng_y'][critical_eng] / adf_data['b'])]])

            X = np.linalg.solve(A,b)

            weight_vector = np.append(weight_vector, W)
            VMCL          = np.append(VMCL, ((X[2])**(-0.5)))
            beta          = np.append(beta, np.degrees(X[0]))
            delta_a       = np.append(delta_a, np.degrees(delta_a_max))
            delta_r       = np.append(delta_r, np.degrees(X[1]))

    # Sort vectors by weight
    ind           = np.argsort(weight_vector)
    weight_vector = weight_vector[ind]
    VMCL          = VMCL[ind]
    beta          = beta[ind]
    delta_a       = delta_a[ind]
    delta_r       = delta_r[ind]

    # Convert VMCL from TAS to KCAS
    for i in range(len(VMCL)):
        VMCL[i] = Speed.tas_to_cas(VMCL[i], altitude, ISA)/kt_to_ms 

    # VMCL for given initial weight (m) (4)
    f       = interp1d(weight_vector, VMCL)
    VMCL_m  = f(adf_data['m'])

    # Plot VMCL for each weight
    if input_data['plot'] == 1: 
        height_in = 11.69
        width_in  = 8.27
        font_size = 15  

        fig_name = out_data['file'][:-4] + '_VMCL.png'
        xl  = 'Weight (kg)'
        yl  = 'Speed (kt)'
        tit = '$V_{MCL}$'

        fig, ax1 = plt.subplots(figsize=(height_in, width_in))
        ax1.plot(weight_vector, VMCL, color = 'tab:blue', label = '$V_{MCL}$')

        ax1.set_xlabel(xl, fontsize=font_size)
        ax1.set_ylabel(yl, fontsize=font_size)
        ax1.set_title(tit, fontsize=font_size)
        
        ax1.grid(b=True, which='major', color='#666666', linestyle='-')
        ax1.minorticks_on()
        ax1.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.5)
        ax1.tick_params(labelsize=font_size)
        ax1.legend(loc='best')

        fig.tight_layout()
        fig.savefig(fig_name)
   
    
    # Merge weight and VMCL in one variable
    VMCL = [weight_vector, VMCL]

    return VMCL, VMCL_m

def write_output_file(out_data, input_data, adf_data, aircraft_config):
    """Function to write the calculated results in an output file

    Args:
        out_data (dict): Dictionary with the output data.
        input_data (dict): Dictionary with the input file variables.

    Returns:
        Output file with the calculated results
    """
    
    lines = '*' * 80 + '\n'
    
    with open(out_data['file'], 'w') as f:

        f.write('\n')
        
        # AIRCRAFT AND CONDITIONS
        f.write(' Aircraft     : ' + input_data['aircraft'] + '\n')
        f.write(' Configuration: ' + '\n')
        
        variable = ['Fuselage', 'Wing', 'Tail', 'Prop_dist', 'Prop_tip', 'Prop_tail','Air_inlet', 'Canard']
        for var in variable:
            f.write(' '*5 + '{:<9s}'.format(var) + ': ' + aircraft_config[var] + '\n')
        if 'Version' in aircraft_config:
            f.write(' '*5 + '{:<9s}'.format('Version') + ': ' + aircraft_config['Version'] + '\n')
        
        f.write(' Condition    : ' + '{:.0f}'.format(input_data['condition']) + '\n')
        f.write('     ISA      : ' + '{:.1f}'.format(adf_data['ISA']) + '\n')
        f.write('     Altitude : ' + '{:.1f}'.format(adf_data['altitude']) + ' ft\n')
        if adf_data['speed'] < 5.0:
            f.write('     Mach     : ' + '{:.2f}'.format(adf_data['Mach']) + '\n')
        else:
            f.write('     KCAS     : ' + '{:.1f}'.format(adf_data['KCAS']) + ' kt\n')
        f.write(' CG position  : ' + '{:.2f}'.format(adf_data['CG']) + '\n')

        f.write('\n')
        
        # STATIC - LONGITUDINAL
        if input_data['type'] == 1 or input_data['type'] == 0:
            f.write(lines)
            f.write('{:^80s}'.format('STATIC - LONGITUDINAL') + '\n')
            f.write(lines)
            
            # Static margin and neutral point
            f.write('{:>15s}'.format('Neutral point :') + '{:8.3f}'.format(out_data['neutral_point']) + '\n')
            f.write('{:>15s}'.format('Static margin :') + '{:8.3f}'.format(out_data['static_margin']) + '\n')
            f.write('\n')

            # Trim
            f.write('Trim:\n')        
            f.write('{:>11s}'.format('Alpha :') + '{:8.3f}'.format(out_data['alpha_trim']) + ' deg \n')
            f.write('{:>11s}'.format('Delta_e :') + '{:8.3f}'.format(out_data['delta_e_trim']) + ' deg \n')
            f.write('\n')
            
            # Control sensitivity
            if 'alpha' in input_data:
                variable = ['alpha', 'delta_e']
                header   = ['Alpha', 'Delta_e']
                units    = [  'deg',     'deg']
                
                f.write('Control sensitivity:\n') 
                
                # Header
                string = ['{:>11s}'.format(s) for s in header]
                string = ''.join(string) + '\n'
                f.write(string)
                
                # Units
                string = ['{:>11s}'.format(s) for s in units]
                string = ''.join(string) + '\n'
                f.write(string)
                
                # Line
                string = '    ' + '-' * (11 * len(header) - 4) + '\n'
                f.write(string)
        
                # Variables
                for i in range(len(input_data['alpha'])):
                    string = input_data[variable[0]][i]
                    string = np.append(string, [out_data[var][i] for var in variable[1:]])
                    string = ['{:>11.3f}'.format(s) for s in string]
                    string = ''.join(string) + '\n'
                    f.write(string) 
                f.write('\n')
                
            f.write('\n')

        # STATIC - LATERAL/DIRECTIONAL     
        if input_data['type'] == 2 or input_data['type'] == 0:
            f.write(lines)
            f.write('{:^80s}'.format('STATIC - LATERAL/DIRECTIONAL') + '\n')
            f.write(lines)
            if 'beta' in input_data:
                variable = ['beta', 'phi_b', 'delta_a_b', 'delta_r_b']
                header   = ['Beta',   'Phi',   'Delta_a',   'Delta_r']
                units    = [ 'deg',   'deg',       'deg',       'deg']
                
                f.write('Trim - sideslip angle:\n') 
                
                # Header
                string = ['{:>11s}'.format(s) for s in header]
                string = ''.join(string) + '\n'
                f.write(string)
                
                # Units
                string = ['{:>11s}'.format(s) for s in units]
                string = ''.join(string) + '\n'
                f.write(string)
                
                # Line
                string = '    ' + '-' * (11 * len(header) - 4) + '\n'
                f.write(string)
        
                # Variables
                for i in range(len(input_data['beta'])):
                    string = input_data[variable[0]][i]
                    string = np.append(string, [out_data[var][i] for var in variable[1:]])
                    string = ['{:>11.3f}'.format(s) for s in string]
                    string = ''.join(string) + '\n'
                    f.write(string) 
                
                f.write('\n')                
                
            if 'phi' in input_data:
                variable = [ 'phi', 'n_p', 'beta_p', 'delta_a_p', 'delta_r_p']
                header   = [ 'Phi',   'n',   'Beta',   'Delta_a',   'Delta_r']
                units    = [ 'deg',   '-',    'deg',       'deg',       'deg']
                
                f.write('Trim - stationary turn:\n') 
                
                # Header
                string = ['{:>11s}'.format(s) for s in header]
                string = ''.join(string) + '\n'
                f.write(string)
                
                # Units
                string = ['{:>11s}'.format(s) for s in units]
                string = ''.join(string) + '\n'
                f.write(string)
                
                # Line
                string = '    ' + '-' * (11 * len(header) - 4) + '\n'
                f.write(string)
        
                # Variables
                for i in range(len(input_data['phi'])):
                    string = input_data[variable[0]][i]
                    string = np.append(string, [out_data[var][i] for var in variable[1:]])
                    string = ['{:>11.3f}'.format(s) for s in string]
                    string = ''.join(string) + '\n'
                    f.write(string) 
                
                f.write('\n')
                
            f.write('\n')
            
        # DYNAMIC - LONGITUDINAL 
        if input_data['type'] == 3 or input_data['type'] == 0:
            f.write(lines)
            f.write('{:^80s}'.format('DYNAMIC - LONGITUDINAL') + '\n')
            f.write(lines)
            
            # Modes
            f.write('Modes:\n')  
            f.write('{:10.5f}'.format(out_data['lambda_ph'][0,0]) + '{:+8.5f}'.format(out_data['lambda_ph'][0,1]) + 'j\n')
            f.write('{:10.5f}'.format(out_data['lambda_ph'][1,0]) + '{:+8.5f}'.format(out_data['lambda_ph'][1,1]) + 'j\n')
            f.write('{:10.5f}'.format(out_data['lambda_sp'][0,0]) + '{:+8.5f}'.format(out_data['lambda_sp'][0,1]) + 'j\n')
            f.write('{:10.5f}'.format(out_data['lambda_sp'][1,0]) + '{:+8.5f}'.format(out_data['lambda_sp'][1,1]) + 'j\n')
            f.write('\n')
            
            # Phugoid
            f.write('Phugoid:\n')
            f.write('{:>11s}'.format('wn :') + '{:8.5f}'.format(out_data['wn_ph']) + ' s\n')
            f.write('{:>11s}'.format('zeta :') + '{:8.5f}'.format(out_data['zeta_ph']) + '\n')
            f.write('{:>11s}'.format('t_1/2 :') + '{:8.3f}'.format(out_data['t_12_ph']) + ' s\n')
            f.write('{:>11s}'.format('T :') + '{:8.3f}'.format(out_data['T_ph']) + ' s\n')
            f.write('\n')
            
            # Short period
            f.write('Short period:\n')
            f.write('{:>11s}'.format('wn :') + '{:8.5f}'.format(out_data['wn_sp']) + ' s\n')
            f.write('{:>11s}'.format('zeta :') + '{:8.5f}'.format(out_data['zeta_sp']) + '\n')
            f.write('{:>11s}'.format('t_1/2 :') + '{:8.3f}'.format(out_data['t_12_sp']) + ' s\n')
            f.write('{:>11s}'.format('T :') + '{:8.3f}'.format(out_data['T_sp']) + ' s\n')
            f.write('\n')
            
            f.write('\n')
            
        # DYNAMIC - LATERAL/DIRECTIONAL 
        if input_data['type'] == 4 or input_data['type'] == 0:
            f.write(lines)
            f.write('{:^80s}'.format('DYNAMIC - LATERAL/DIRECTIONAL') + '\n')
            f.write(lines)
            
            # Modes
            f.write('Modes:\n')  
            f.write('{:10.5f}'.format(out_data['lambda_rs'][0]) + '{:+8.5f}'.format(out_data['lambda_rs'][1]) + 'j\n')
            f.write('{:10.5f}'.format(out_data['lambda_spi'][0]) + '{:+8.5f}'.format(out_data['lambda_spi'][1]) + 'j\n')
            f.write('{:10.5f}'.format(out_data['lambda_dr'][0,0]) + '{:+8.5f}'.format(out_data['lambda_dr'][0,1]) + 'j\n')
            f.write('{:10.5f}'.format(out_data['lambda_dr'][1,0]) + '{:+8.5f}'.format(out_data['lambda_dr'][1,1]) + 'j\n')
            f.write('\n')
            
            # Roll subsidence
            f.write('Roll subsidence:\n')
            f.write('{:>11s}'.format('lambda :') + '{:8.5f}'.format(out_data['lambda_rs'][0]) + ' s\n')
            f.write('{:>11s}'.format('t_1/2 :') + '{:8.3f}'.format(out_data['t_12_rs']) + ' s\n')
            f.write('\n')
            
            # Spiral
            f.write('Spiral:\n')
            f.write('{:>11s}'.format('lambda :') + '{:8.5f}'.format(out_data['lambda_spi'][0]) + ' s\n')
            f.write('{:>11s}'.format('t_1/2 :') + '{:8.3f}'.format(out_data['t_12_spi']) + ' s\n')
            f.write('\n')
            
            # Dutch roll
            f.write('Dutch roll:\n')
            f.write('{:>11s}'.format('wn :') + '{:8.5f}'.format(out_data['wn_dr']) + ' s\n')
            f.write('{:>11s}'.format('zeta :') + '{:8.5f}'.format(out_data['zeta_dr']) + ' s\n')
            f.write('{:>11s}'.format('t_1/2 :') + '{:8.3f}'.format(out_data['t_12_dr']) + ' s\n')
            f.write('{:>11s}'.format('T :') + '{:8.3f}'.format(out_data['T_dr']) + ' s\n')
            f.write('\n')
            
            f.write('\n')

        # MINIMUM SPEEDS
        if input_data['type'] == 5 or input_data['type'] == 0:
            f.write(lines)
            f.write('{:^80s}'.format('MINIMUM SPEEDS') + '\n')
            f.write(lines)

            # Values for current weight
            f.write('Values for current weight\n')

            # VSR
            f.write('{:>11s}'.format('VSR :') + '{:8.3f}'.format(out_data['Vsr']/kt_to_ms) + ' kt\n')

            # VMCG
            if adf_data['TO_condition'] == 0:
                f.write('      VMCG not caclculated: There is no data for take-off condition\n')
            else:
                f.write('{:>11s}'.format('VMCG :') + '{:8.3f}'.format(out_data['VMCG']) + ' kt\n')
            
            # VMCA
            if adf_data['TO_condition'] == 0:
                f.write('      VMCA not caclculated: There is no data for take-off condition\n')
            else:
                f.write('{:>11s}'.format('VMCA :') + '{:8.3f}'.format(out_data['VMCA_m']) + ' kt\n')
            
            # VMCL
            if adf_data['LD_condition'] == 0:
                f.write('      VMCL not caclculated: There is no data for landing condition\n')
            else:
                f.write('{:>11s}'.format('VMCL :') + '{:8.3f}'.format(out_data['VMCL_m']) + ' kt\n')

            f.write('\n')

            # VMCA table
            if adf_data['TO_condition'] != 0:
                header   = ['Weight', 'VMCA', '1.13VSR']
                units    = [    'kg',   'kt',       '-']

                f.write('VMCA vs weight:\n')

                # Header
                string = ['{:>11s}'.format(s) for s in header]
                string = ''.join(string) + '\n'
                f.write(string)

                # Units
                string = ['{:>11s}'.format(s) for s in units]
                string = ''.join(string) + '\n'
                f.write(string)

                # Line
                string = ' ' + '-' * (11 * len(header) - 1) + '\n'
                f.write(string)

                # Variables
                n_row = len(out_data['VMCA'][0])
                n_col = len(out_data['VMCA'])
                for row in range(n_row):
                    string = []
                    string = np.append(string, [out_data['VMCA'][col][row] for col in range(n_col)])
                    string = ['{:>11.3f}'.format(s) for s in string]
                    string = ''.join(string) + '\n'
                    f.write(string) 

                f.write('\n')

            # VMCL table
            if adf_data['LD_condition'] != 0:
                header   = [        'Weight', 'VMCL']
                units    = [            'kg',   'kt']

                f.write('VMCL vs weight:\n')

                # Header
                string = ['{:>11s}'.format(s) for s in header]
                string = ''.join(string) + '\n'
                f.write(string)

                # Units
                string = ['{:>11s}'.format(s) for s in units]
                string = ''.join(string) + '\n'
                f.write(string)

                # Line
                string = ' ' + '-' * (11 * len(header) - 1) + '\n'
                f.write(string)

                # Variables
                n_row = len(out_data['VMCL'][0])
                n_col = len(out_data['VMCL'])
                for row in range(n_row):
                    string = []
                    string = np.append(string, [out_data['VMCL'][col][row] for col in range(n_col)])
                    string = ['{:>11.3f}'.format(s) for s in string]
                    string = ''.join(string) + '\n'
                    f.write(string) 

                f.write('\n')
                                  
"""
MAIN PROGRAM

"""

"""
 Global variables:
     feet_to_meters: Feet to meters conversion factor
     kt_to_ms: Knots to meters per second conversion factor
     g: Acceleration of gravity

"""
feet_to_meters = 0.3048
kt_to_ms       = 0.514444
g              = 9.80665

################################################################################
#                                    INPUT                                     #
################################################################################
# Define input file path
file = './input_file_template.dat'

################################################################################

# Read input and aircraft data
input_data, adf_data, aircraft_config = read_input_file(file)

# Calculate atmosphere for the study conditions
atmos = Atmos(feet_to_meters*adf_data['altitude'], adf_data['ISA'], 1)

# Initialize output dictionary
out_data = {}
out_data['file']     = file[:-4] + '.out'

# CALCULATE
# Static - longitudinal
if input_data['type'] == 1 or input_data['type'] == 0:
    out_data['static_margin'], out_data['neutral_point'], out_data['alpha_trim'], out_data['delta_e_trim'], out_data['delta_e'] = calculate_static_long(adf_data, input_data)
            
# Static - lateral-directional   
if input_data['type'] == 2 or input_data['type'] == 0:
    out_data['phi_b'], out_data['delta_a_b'], out_data['delta_r_b'], \
    out_data['beta_p'], out_data['delta_a_p'], out_data['delta_r_p'], out_data['n_p'] = calculate_static_latdir(adf_data, input_data)

# Dynamic - longitudinal 
if input_data['type'] == 3 or input_data['type'] == 0:
    out_data['lambda_ph'], out_data['wn_ph'], out_data['zeta_ph'], out_data['t_12_ph'], out_data['T_ph'], \
    out_data['lambda_sp'], out_data['wn_sp'], out_data['zeta_sp'], out_data['t_12_sp'], out_data['T_sp']  = calculate_dynamic_long(adf_data, input_data)    

# Dynamic - longitudinal 
if input_data['type'] == 4 or input_data['type'] == 0:
    out_data['lambda_rs'], out_data['t_12_rs'], \
    out_data['lambda_spi'], out_data['t_12_spi'],\
    out_data['lambda_dr'], out_data['wn_dr'], out_data['zeta_dr'], out_data['t_12_dr'], out_data['T_dr'] = calculate_dynamic_latdir(adf_data, input_data)    

# Minimum  speeds
if input_data['type'] == 5 or input_data['type'] == 0: 
    out_data['Vsr']  = get_stall_speed(adf_data['m'], atmos.density, adf_data['S'], adf_data['CL_max'])
    out_data['VMCG'] = calculate_vmcg(adf_data)
    out_data['VMCA'], out_data['VMCA_m'] = calculate_vmca(adf_data, input_data)
    out_data['VMCL'], out_data['VMCL_m'] = calculate_vmcl(adf_data, input_data)
      
# PRINT RESULTS   
write_output_file(out_data, input_data, adf_data, aircraft_config)

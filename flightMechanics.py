# -*- coding: utf-8 -*-
"""
   FLIGHT MECHANICS
"""
import math
import numpy as np
import matplotlib.pyplot as plt


def read_input_file(file):
    """ 
    This function reads the input file into a dict
    
    Inputs:
        file: Path to the input file
        
    Outputs:
        input_data: Dictionary with the input file variables
        adf_data: Dictionary with the aircraft data
    
    """
    with open(file,'r') as f:
        lines = [line.strip().split(':') for line in f.readlines() if (line.strip() and not line.startswith("#"))]
        
    variable_name = [line[0].strip() for line in lines]
    variable_data = [line[1].replace(" ","").split(',') for line in lines]
    
    input_data = {}   
    for index, title in enumerate(variable_name):    
        if title == 'aircraft':
            input_data[title] = variable_data[index][0]
        elif title == 'type':
            input_data[title] = np.array([int(var) for var in variable_data[index]])[0]
        else:
            input_data[title] = np.array([float(var) for var in variable_data[index]])
                
    adf_data = read_adf_file(input_data['aircraft'])
    
    return input_data, adf_data

def read_adf_file(file):
    """ 
    This function reads the adf file into a dict
    
    Inputs:
        file: Path to the adf file
        
    Outputs:
        adf_data: Dictionary with the aircraft data
    
    """          
    with open(file,'r') as f:
        lines = [line.strip().split(':') for line in f.readlines() if (line.strip() and not line.startswith("#"))]
        
    variable_name = [line[0].strip() for line in lines]
    variable_data = [[float(var) for var in line[1].replace(" ","").split(',')] for line in lines]

    adf_data = {}      
    for index, title in enumerate(variable_name):
        if (title == 'T' or title == 'epsilon' or title == 'ni' or 
            title == 'd_CG_x' or title == 'd_CG_y' or title == 'd_CG_z'):
            adf_data[title] = np.array(variable_data[index])
        else:
            adf_data[title] = np.array(variable_data[index])[0]
    
    # Get speed in Mach and KTAS
    a = get_speed_of_sound(feet_to_meters * adf_data['altitude'])
    if adf_data['speed'] < 5.0 :
        adf_data['Mach'] = adf_data['speed']
        adf_data['KTAS'] = (adf_data['Mach'] * a) / kt_to_ms
    else:
        adf_data['KTAS'] = adf_data['speed']
        adf_data['Mach'] = (kt_to_ms * adf_data['KTAS']) / a 
        
    # Get Cz_s for the stability condition
    rho = get_density(feet_to_meters * adf_data['altitude'])
    p_d = get_dynamic_pressure(rho, kt_to_ms * adf_data['KTAS'])
    adf_data['Cz_s'] = -(adf_data['m'] * g) / (p_d * adf_data['S'])
        
    
    return adf_data

def get_temperature (altitude):
    """
    This function calculates the temperature for a given altitude
    
    Inputs:
        altitude: Flight altitude (m)
        
    Outputs:
        T: Temperature (K)
    
    """
    alpha_T  = 6.5e-3
    T_0     = 288.15

    if altitude <= 11000.0:    
        T = T_0 - alpha_T * altitude
    else:
        T = T_0 - alpha_T * 11000.0
    
    return T

def get_density (altitude):
    """
    This function calculates the air density for a given altitude
    
    Inputs:
        altitude: Flight altitude (m)
        
    Outputs:
        rho: Air density (kg/m3)
    
    """
    rho_0   = 1.225
    alpha_T  = 6.5e-3
    T_0     = 288.15
    Ra      = 287.05
    
    if altitude <= 11000.0:    
        rho = rho_0 * ((1 - (alpha_T * altitude)/T_0) ** ((g/(Ra * alpha_T)) - 1))
    else:
        rho_11 = rho_0 * ((1 - (alpha_T * 11000.0)/T_0) ** ((g/(Ra * alpha_T)) - 1))
        T_11   = 216.65
        rho    = rho_11 * math.e ** (-g * (altitude - 11000.0) / (Ra * T_11))
    
    return rho

def get_speed_of_sound (altitude):
    """
    This function calculates the speed of sound for a given altitude
    
    Inputs:
        altitude: Flight altitude (m)
        
    Outputs:
        a: speed of sound (m/s)
    
    """
    gamma = 1.4
    R     = 286.0
    
    T = get_temperature(altitude)
    
    a = (gamma * R * T) ** 0.5
    
    return a

def get_dynamic_pressure(rho,speed):
    """
    This function calculates the dynamic pressure for a given density and speed
    
    Inputs:
        rho: Air density (kg/m3)
        speed : Speed (m/s)
        
    Outputs:
        p_d: Dynamic pressure (N/m2)
    
    """
    
    p_d = 0.5 * rho * speed ** 2
    
    return p_d

def calculate_static_long(adf_data, input_data):
    """
    This function calculates the static longitudinal stability and control
    
    Inputs:
        adf_data: dictionary with the aircraft data
        input_data: dictionary with the input data
        
    Outputs:
        alpha_trim, delta_e_trim : Trim values (deg)
        alpha, delta_e : Elevator deflection vs. alpha (deg)   
        
    """  

    # Intermediate calculations   
    rho = get_density(feet_to_meters * adf_data['altitude'])
    p_d = get_dynamic_pressure(rho, kt_to_ms * adf_data['KTAS'])
    
    """
    # Propulsion
    
    n_engines = len(adf_data['T'])
    
    Ft_z = sum(-adf_data['T'][i] * math.sin(np.radians(adf_data['epsilon'][i]))  for i in range(n_engines))

    
    Mt_y = sum(adf_data['T'][i] * (math.cos(np.radians(adf_data['epsilon'][i])) * math.cos(np.radians(adf_data['ni'][i])) * adf_data['d_CG_z'][i] +
               math.sin(np.radians(adf_data['epsilon'][i])) * adf_data['d_CG_x'][i]) for i in range(n_engines))
    """
    
    # TRIM
    CL_alpha   = - adf_data['Cz_alpha']
    CL_delta_e = - adf_data['Cz_delta_e']
    CL_0       = - adf_data['Cz_0']
    
    A = np.array([[            CL_alpha,             CL_delta_e],
                  [adf_data['Cm_alpha'], adf_data['Cm_delta_e']]])
    
    """
    B = np.array([[(adf_data['m'] * g + Ft_z) / (p_d * adf_data['S']) - adf_data['Cz_0']],
                  [-Mt_y / (p_d * adf_data['S'] * adf_data['c']) -adf_data['Cm_0']]])
    """
    B = np.array([[(adf_data['m'] * g) / (p_d * adf_data['S']) - CL_0],
                  [-adf_data['Cm_0']]])

    X = np.linalg.solve(A,B)
    
    alpha_trim   = np.degrees(X[0])[0]
    delta_e_trim = np.degrees(X[1])[0]
    
    # CONTROL
    # Cinematic momentum in Y axis
    #M_y = 0.0

    # Calculate elevator deflection for balance
    if 'alpha' in input_data:
        delta_e = np.array([])
        for alpha in input_data['alpha']:
            #d_e = ((M_y - Mt_y) / (p_d * adf_data['S'] * adf_data['c']) - adf_data['Cm_0'] - adf_data['Cm_alpha'] * (np.radians(alpha))) / adf_data['Cm_delta_e']
            d_e = (- adf_data['Cm_0'] - adf_data['Cm_alpha'] * (np.radians(alpha))) / adf_data['Cm_delta_e']
    
            delta_e = np.append(delta_e, np.degrees(d_e))
    else:
        delta_e = np.nan
        
    return alpha_trim, delta_e_trim, delta_e

def calculate_static_latdir(adf_data, input_data):
    """
    This function calculates the static lateral-directionsl stability 
    and control
    
    Inputs:
        adf_data: dictionary with the aircraft data
        input_data: dictionary with the input data
        
    Outputs:
        phi_b, delta_a_b, delta_r_b: Trim values for a given beta
        beta_p, delta_a_p, delta_r_p, n_p: Trim values for a stationary turn (fixed phi)
    
    """  
    
    # Intermediate calculations   
    rho = get_density(feet_to_meters * adf_data['altitude'])
    p_d = get_dynamic_pressure(rho, kt_to_ms * adf_data['KTAS'])
    
    # Propulsion
    n_engines = len(adf_data['T'])
    
    Ft_y = sum(input_data['T_rate'][i] * adf_data['T'][i] * math.cos(np.radians(adf_data['epsilon'][i])) * math.sin(np.radians(adf_data['ni'][i])) for i in range(n_engines))
    
    Lt   = sum(-input_data['T_rate'][i] * adf_data['T'][i] * (math.cos(np.radians(adf_data['epsilon'][i])) * math.sin(np.radians(adf_data['ni'][i])) * adf_data['d_CG_z'][i] +
             math.sin(np.radians(adf_data['epsilon'][i])) * adf_data['d_CG_y'][i]) for i in range(n_engines))
    
    Nt   = sum(input_data['T_rate'][i] * adf_data['T'][i] * (- math.cos(np.radians(adf_data['epsilon'][i])) * math.cos(np.radians(adf_data['ni'][i])) * adf_data['d_CG_y'][i] +
             math.cos(np.radians(adf_data['epsilon'][i])) * math.sin(np.radians(adf_data['ni'][i])) * adf_data['d_CG_x'][i]) for i in range(n_engines))
    
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

def calculate_dynamic_long(adf_data):
    """
    This function calculates the dynamic longitudinal stability
    
    Inputs:
        adf_data: dictionary with the aircraft data

    Outputs:
        wn_ph, zeta_ph, lambda_ph: results for phugoide mode
        wn_sp, zeta_sp, lambda_sp: results for short period mode
    
    """  
    
    rho = get_density(feet_to_meters * adf_data['altitude'])
    
    mu = adf_data['m'] / (0.5 * rho * adf_data['S'] *  adf_data['c'])
    
    Iy_ad = adf_data['Iy'] / (rho * adf_data['S'] * (adf_data['c'] / 2) ** 3)
    
    # PHUGOID
    wn_ph_ad = ((adf_data['Cz_s'] * (2 * adf_data['Cz_s'] + adf_data['Cz_u'])) / (2 * mu * (2 * mu + adf_data['Cz_q']))) ** 0.5
    wn_ph    = wn_ph_ad / (adf_data['c'] / (2 * kt_to_ms * adf_data['KTAS']))
    
    zeta_ph  = -(adf_data['Cx_u'] * (2 * mu + adf_data['Cz_q'])) / (2 * (adf_data['Cz_s'] * (2 * adf_data['Cz_s'] + adf_data['Cz_u']) * 2 * mu * (2 * mu + adf_data['Cz_q'])) ** 0.5)
    
    lambda_r  = -zeta_ph * wn_ph
    lambda_i  = wn_ph * (1 - zeta_ph ** 2) ** 0.5
    lambda_ph = np.array([[lambda_r, lambda_i], [lambda_r, -lambda_i]])
    
    #wn2 = 2 ** 0.5 * g / ( kt_to_ms * adf_data['KTAS']) 
    #zeta2 = adf_data['Cx_u']/(2 * 2 ** 0.5 * adf_data['Cz_s'])
    
    # SHORT PERIOD
    wn_sp_ad = ((adf_data['Cz_alpha'] * adf_data['Cm_q'] - (2 * mu + adf_data['Cz_q']) * adf_data['Cm_alpha']) / ((2 * mu - adf_data['Cz_alpha_dot']) * Iy_ad)) ** 0.5
    wn_sp    = wn_sp_ad / (adf_data['c'] / (2 * kt_to_ms * adf_data['KTAS']))
    
    zeta_sp = -((2 * mu - adf_data['Cz_alpha_dot']) * adf_data['Cm_q'] + adf_data['Cz_alpha'] * Iy_ad + (2 * mu + adf_data['Cz_q']) * adf_data['Cm_alpha_dot']) / (2 * ((2 * mu - adf_data['Cz_alpha_dot']) * Iy_ad * (adf_data['Cz_alpha'] * adf_data['Cm_q'] - (2 * mu + adf_data['Cz_q']) * adf_data['Cm_alpha'])) ** 0.5)
    
    lambda_r  = -zeta_sp * wn_sp
    lambda_i  = wn_sp * (1 - zeta_sp ** 2) ** 0.5
    lambda_sp = np.array([[lambda_r, lambda_i], [lambda_r, -lambda_i]])

    return wn_ph, zeta_ph, lambda_ph, wn_sp, zeta_sp, lambda_sp 

def calculate_dynamic_latdir(adf_data):
    """
    This function calculates the dynamic lateral-directional stability
    
    Inputs:
        adf_data: dictionary with the aircraft data

    Outputs:
        wn_rs, lambda_rs: results for roll subsidence mode
        wn_spi, lambda_spi: results for spiral mode
        wn_dr, zeta_dr, lambda_dr: results for dutch roll mode
    
    """
    
    rho = get_density(feet_to_meters * adf_data['altitude'])
    
    mu = adf_data['m'] / (0.5 * rho * adf_data['S'] *  adf_data['b'])
    
    Ix_ad = adf_data['Ix'] / (rho * adf_data['S'] * (adf_data['b'] / 2) ** 3)
    Iz_ad = adf_data['Iz'] / (rho * adf_data['S'] * (adf_data['b'] / 2) ** 3)
    
    # ROLL SUBSIDENCE
    lambda_rs_ad = adf_data['Cl_p'] / Ix_ad 
    lambda_rs    = lambda_rs_ad/ (adf_data['b'] / (2 * kt_to_ms * adf_data['KTAS']))
    
    wn_rs = lambda_rs
    
    # SPIRAL
    lambda_spi_ad = -(adf_data['Cz_s'] * (adf_data['Cl_beta'] * adf_data['Cn_r'] - adf_data['Cn_beta'] * adf_data['Cl_r'])) / ((2 * mu) * (adf_data['Cn_beta'] * adf_data['Cl_p'] - adf_data['Cl_beta'] * adf_data['Cn_p']) + adf_data['Cy_beta'] * (adf_data['Cl_p'] * adf_data['Cn_r'] - adf_data['Cn_p'] * adf_data['Cl_r']))
    lambda_spi    = lambda_spi_ad/ (adf_data['b'] / (2 * kt_to_ms * adf_data['KTAS']))
    
    wn_spi = lambda_spi
    
    # DUTCH ROLL
    wn_dr_ad = (adf_data['Cn_beta'] / Iz_ad) ** 0.5
    wn_dr    = wn_dr_ad / (adf_data['b'] / (2 * kt_to_ms * adf_data['KTAS']))
    
    zeta_dr = -adf_data['Cn_r'] / (2 * (Iz_ad * adf_data['Cn_beta']) ** 0.5)
    
    lambda_r  = -zeta_dr * wn_dr
    lambda_i  = wn_dr * (1 - zeta_dr ** 2) ** 0.5
    lambda_dr = np.array([[lambda_r, lambda_i], [lambda_r, -lambda_i]])

    return wn_rs, lambda_rs, wn_spi, lambda_spi, wn_dr, zeta_dr, lambda_dr
    
def plot_out_static(out_data):
    """ 
    This function plot the calculated outputs
    
    Inputs:
        out_data: Dictionary containing the results
        
    Outputs:
        plots
    
    """  

    """
    TO BE UPDATED
    """
    
    for alt in np.unique(out_data['altitude']): # Plot for each altitude
        condition_alt = out_data['altitude'] == alt
        for spd in np.unique(out_data['V']):    # Plot for each speed
            condition_spd = out_data['V'] == spd
            
            condition = np.logical_and(condition_alt, condition_spd)
            
            fig, ax1 = plt.subplots()
            
            # Longitudinal
            if out_data['type'] == 1: 
                alpha   = np.extract(condition, out_data['alpha'])
                delta_e = np.extract(condition, out_data['delta_e'])
                                             
                fig_name = out_data['file'][:-4] + '_alt' + str(int(alt)) + '_spd' + str(int(spd)) + '.png'
                xl  = 'Alpha (deg)'
                yl  = 'Elevator deflection (deg)'
                tit = 'Altitude: ' + str(int(alt)) + ' ft\n' + 'Speed: ' + str(int(spd)) + ' kt'
                
                ax1.plot(alpha, delta_e)
                
                ax1.set(xlabel = xl, ylabel = yl, title = tit)
                ax1.grid()
                
                fig.tight_layout()
                fig.savefig(fig_name)
                plt.show()
                
            # Lateral-directional
            if out_data['type'] == 2: 
                beta    = np.extract(condition, out_data['beta'])
                delta_a = np.extract(condition, out_data['delta_a'])
                delta_r = np.extract(condition, out_data['delta_r'])
                
                fig_name = out_data['file'][:-4] + '_alt' + str(int(alt)) + '_spd' + str(int(spd)) + '.png'
                xl  = 'Beta (deg)'
                yl  = 'Control deflection (deg)'
                tit = 'Altitude: ' + str(int(alt)) + ' ft\n' + 'Speed: ' + str(int(spd)) + ' kt'
                
                ax1.plot(beta, delta_a, color = 'tab:blue', label = 'Aileron')
                ax1.plot(beta, delta_r, color = 'tab:red', label = 'Rudder')
                
                ax1.set(xlabel = xl, ylabel = yl, title = tit)
                ax1.grid()
                ax1.legend(loc='best')
                
                fig.tight_layout()
                fig.savefig(fig_name)
                plt.show()

def write_output_file(out_data, input_data):
    """ 
    This function writes an output file
    ​
    Inputs: 
    out_data: dictionary with the variables already calculated
    ​
    Outputs:
    output_file: file with the calculated values
    ​
    """
    
    lines = '*' * 80 + '\n'
    
    with open(out_data['file'], 'w') as f:

        f.write('\n')
        
        # AIRCRAFT AND CONDITIONS
        f.write(' Aircraft : ' + input_data['aircraft'] + '\n')
        f.write(' Altitude : ' + '{:.1f}'.format(out_data['altitude']) + ' ft\n')
        if out_data['speed'] < 5.0:
            f.write('    Speed : M' + '{:.2f}'.format(out_data['speed']) + '\n')
        else:
            f.write('    Speed : ' + '{:.1f}'.format(out_data['speed']) + ' kt\n')
        
        f.write('\n')
        
        # STATIC - LONGITUDINAL
        if input_data['type'] == 1 or input_data['type'] == 0:
            f.write(lines)
            f.write('{:^80s}'.format('STATIC - LONGITUDINAL') + '\n')
            f.write(lines)
            
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
            
            # Phugoid
            f.write('Modes:\n')  
            f.write('{:10.4f}'.format(out_data['lambda_ph'][0,0]) + '{:+8.4f}'.format(out_data['lambda_ph'][0,1]) + 'j\n')
            f.write('{:10.4f}'.format(out_data['lambda_ph'][1,0]) + '{:+8.4f}'.format(out_data['lambda_ph'][1,1]) + 'j\n')
            f.write('{:10.4f}'.format(out_data['lambda_sp'][0,0]) + '{:+8.4f}'.format(out_data['lambda_sp'][0,1]) + 'j\n')
            f.write('{:10.4f}'.format(out_data['lambda_sp'][1,0]) + '{:+8.4f}'.format(out_data['lambda_sp'][1,1]) + 'j\n')
            f.write('\n')
            
            f.write('Phugoid:\n')
            f.write('{:>11s}'.format('wn :') + '{:8.5f}'.format(out_data['wn_ph']) + '\n')
            f.write('{:>11s}'.format('zeta :') + '{:8.5f}'.format(out_data['zeta_ph']) + '\n')
            f.write('\n')
            
            f.write('Short period:\n')
            f.write('{:>11s}'.format('wn :') + '{:8.5f}'.format(out_data['wn_sp']) + '\n')
            f.write('{:>11s}'.format('zeta :') + '{:8.5f}'.format(out_data['zeta_sp']) + '\n')
            f.write('\n')
            
            f.write('\n')
            
        # DYNAMIC - LONGITUDINAL 
        if input_data['type'] == 4 or input_data['type'] == 0:
            f.write(lines)
            f.write('{:^80s}'.format('DYNAMIC - LATERAL/DIRECTIONAL') + '\n')
            f.write(lines)
            
            # Phugoid
            f.write('Modes:\n')  
            f.write('{:10.4f}'.format(out_data['lambda_rs']) + '{:+8.4f}'.format(0.0) + 'j\n')
            f.write('{:10.4f}'.format(out_data['lambda_spi']) + '{:+8.4f}'.format(0.0) + 'j\n')
            f.write('{:10.4f}'.format(out_data['lambda_dr'][0,0]) + '{:+8.4f}'.format(out_data['lambda_dr'][0,1]) + 'j\n')
            f.write('{:10.4f}'.format(out_data['lambda_dr'][1,0]) + '{:+8.4f}'.format(out_data['lambda_dr'][1,1]) + 'j\n')
            f.write('\n')
            
            f.write('Roll subsidence:\n')
            f.write('{:>11s}'.format('wn :') + '{:8.5f}'.format(out_data['wn_rs']) + '\n')
            f.write('\n')
            
            f.write('Spiral:\n')
            f.write('{:>11s}'.format('wn :') + '{:8.5f}'.format(out_data['wn_spi']) + '\n')
            f.write('\n')
            
            f.write('Dutch roll:\n')
            f.write('{:>11s}'.format('wn :') + '{:8.5f}'.format(out_data['wn_dr']) + '\n')
            f.write('{:>11s}'.format('zeta :') + '{:8.5f}'.format(out_data['zeta_dr']) + '\n')
            f.write('\n')
            
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

# Define input file path
file = './study/HS.dat'

# Read input and aircraft data
input_data, adf_data = read_input_file(file)

# Initialize output dictionary
out_data = {}
out_data['file']     = file[:-4] + '.out'
out_data['altitude'] = adf_data['altitude']
out_data['speed']    = adf_data['speed']

# CALCULATE
# Static - longitudinal
if input_data['type'] == 1 or input_data['type'] == 0:
    out_data['alpha_trim'], out_data['delta_e_trim'], out_data['delta_e'] = calculate_static_long(adf_data, input_data)
            
# Static - lateral-directional   
if input_data['type'] == 2 or input_data['type'] == 0:
    out_data['phi_b'], out_data['delta_a_b'], out_data['delta_r_b'], out_data['beta_p'], out_data['delta_a_p'], out_data['delta_r_p'], out_data['n_p'] = calculate_static_latdir(adf_data, input_data)

# Dynamic - longitudinal 
if input_data['type'] == 3 or input_data['type'] == 0:
    out_data['wn_ph'], out_data['zeta_ph'], out_data['lambda_ph'], out_data['wn_sp'], out_data['zeta_sp'], out_data['lambda_sp'] = calculate_dynamic_long(adf_data)    

# Dynamic - longitudinal 
if input_data['type'] == 4 or input_data['type'] == 0:
    out_data['wn_rs'], out_data['lambda_rs'], out_data['wn_spi'], out_data['lambda_spi'], out_data['wn_dr'], out_data['zeta_dr'], out_data['lambda_dr'] = calculate_dynamic_latdir(adf_data)    
    
      
# Print results   
write_output_file(out_data, input_data)

"""
# Plot results
if input_data['type'] == 1 or input_data['type'] == 2: # Static    
    plot_out_static(out_data)
else: # Dynamic
    print('Dynamic stability and control not available for plotting')
"""

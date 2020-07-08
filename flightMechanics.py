# -*- coding: utf-8 -*-
"""
   FLIGHT MECHANICS
"""
import math
import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, Poly
from scipy import signal
from scipy.optimize import fsolve


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
            title == 'eng_x' or title == 'eng_y' or title == 'eng_z' or
            title == 'lg_x' or title == 'lg_y' or title == 'lg_z'):
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
    
    # Temporary 
    adf_data['Cm_delta_e_dot'] = 0.0 
    adf_data['Cl_delta_a_dot'] = 0.0 
    adf_data['Cn_delta_r_dot'] = 0.0 
    
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

def get_stall_speed(adf_data):
    """
    This function calculates the reference stall speed (Vsr)
    
    Inputs:
        adf_data: dictionary with the aircraft data
        
    Outputs:
        Vsr: reference stall speed (m/s)

    """
    rho = get_density(feet_to_meters*adf_data['altitude'])
    Vsr = ((adf_data['m']*g)/(0.5*rho*adf_data['S']*adf_data['CL_max']))**0.5

    return Vsr

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

    
    Mt_y = sum(adf_data['T'][i] * (math.cos(np.radians(adf_data['epsilon'][i])) * math.cos(np.radians(adf_data['ni'][i])) * adf_data['eng_z'][i] +
               math.sin(np.radians(adf_data['epsilon'][i])) * adf_data['eng_x'][i]) for i in range(n_engines))
    """
    
    # TRIM
    CL_alpha   = - adf_data['Cz_alpha']
    CL_delta_e = - adf_data['Cz_delta_e']
    CL_0       = - adf_data['Cz_0']
    
    A = np.array([[            CL_alpha,             CL_delta_e],
                  [adf_data['Cm_alpha'], adf_data['Cm_delta_e']]])
    
    """
    B = np.array([[(adf_data['m'] * g + Ft_z) / (p_d * adf_data['S']) - adf_data['CL_0']],
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
    
    Lt   = sum(-input_data['T_rate'][i] * adf_data['T'][i] * (math.cos(np.radians(adf_data['epsilon'][i])) * math.sin(np.radians(adf_data['ni'][i])) * adf_data['eng_z'][i] +
             math.sin(np.radians(adf_data['epsilon'][i])) * adf_data['eng_y'][i]) for i in range(n_engines))
    
    Nt   = sum(input_data['T_rate'][i] * adf_data['T'][i] * (- math.cos(np.radians(adf_data['epsilon'][i])) * math.cos(np.radians(adf_data['ni'][i])) * adf_data['eng_y'][i] +
             math.cos(np.radians(adf_data['epsilon'][i])) * math.sin(np.radians(adf_data['ni'][i])) * adf_data['eng_x'][i]) for i in range(n_engines))
    
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
        
        control response charts
    
    """  
    
    rho = get_density(feet_to_meters * adf_data['altitude'])
    
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
    
    #SHORT PERIOD
    lambda_sp = np.array([[L.real[ind_sp],  np.absolute(L.imag[ind_sp])],
                          [L.real[ind_sp], -np.absolute(L.imag[ind_sp])]])
    
    wn_sp   = (lambda_sp[0,0] ** 2 + lambda_sp[0,1] ** 2) ** 0.5
    zeta_sp = -lambda_sp[0,0] /((lambda_sp[0,0] ** 2 + lambda_sp[0,1] ** 2) ** 0.5)
    
    """ 
    CONTROL
    
    """  
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
        
        control response charts
    
    """
    
    rho = get_density(feet_to_meters * adf_data['altitude'])
    
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
    
    wn_rs   = lambda_rs[0]
    
    # SPIRAL
    lambda_spi = np.array([L_aux[ind_spi],  0.0])
    
    wn_spi   = lambda_spi[0]
    
    # DUTCH ROLL
    lambda_dr = np.array([[L.real[ind_dr],  np.absolute(L.imag[ind_dr])],
                          [L.real[ind_dr], -np.absolute(L.imag[ind_dr])]])
    
    wn_dr   = (lambda_dr[0,0] ** 2 + lambda_dr[0,1] ** 2) ** 0.5
    zeta_dr = -lambda_dr[0,0] /((lambda_dr[0,0] ** 2 + lambda_dr[0,1] ** 2) ** 0.5)
    
    """ 
    CONTROL
    
    """
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


    return wn_rs, lambda_rs, wn_spi, lambda_spi, wn_dr, zeta_dr, lambda_dr

def calculate_vmcg(adf_data):
    """
    This function calculates the minimum control speed on ground (VMCG)
    For the Take-off thrust at SL and ISA conditions
    
    Inputs:
        adf_data: dictionary with the aircraft data
        
    Outputs:
        V: VMC (kt)

    """
    delta_r = -np.radians(adf_data['delta_r_max']) # Right engine out
    rho     = get_density(feet_to_meters*adf_data['altitude'])
    n_lg    = len(adf_data['lg_x'])
    Vsr     = get_stall_speed(adf_data)

    # Propulsion
    n_engines = len(adf_data['T']) - 1 # Right engine out
    
    Ft_y = sum(adf_data['T'][i] * math.cos(np.radians(adf_data['epsilon'][i])) * math.sin(np.radians(adf_data['ni'][i])) for i in range(n_engines))
    Ft_z = sum(-adf_data['T'][i] * math.sin(np.radians(adf_data['epsilon'][i]))  for i in range(n_engines))

    Lt = sum(adf_data['T'][i] * (math.cos(np.radians(adf_data['epsilon'][i])) * math.sin(np.radians(adf_data['ni'][i])) * adf_data['eng_z'][i] +
             math.sin(np.radians(adf_data['epsilon'][i])) * adf_data['eng_y'][i]) for i in range(n_engines))
    Mt = sum(adf_data['T'][i] * (math.cos(np.radians(adf_data['epsilon'][i])) * math.cos(np.radians(adf_data['ni'][i])) * adf_data['eng_z'][i] +
             math.sin(np.radians(adf_data['epsilon'][i])) * adf_data['eng_x'][i]) for i in range(n_engines))
    Nt = sum(adf_data['T'][i] * (- math.cos(np.radians(adf_data['epsilon'][i])) * math.cos(np.radians(adf_data['ni'][i])) * adf_data['eng_y'][i] +
             math.cos(np.radians(adf_data['epsilon'][i])) * math.sin(np.radians(adf_data['ni'][i])) * adf_data['eng_x'][i]) for i in range(n_engines))
    
    # First aproximation
    def equations_aprox(var):
        V, beta = var

        eq2 = Ft_y + 0.5*rho*V**2*adf_data['S']*(adf_data['Cy_beta']*beta + adf_data['Cy_delta_r']*delta_r)
        eq6 = Nt + 0.5*rho*V**2*adf_data['S']*adf_data['b']*(adf_data['Cn_beta']*beta + adf_data['Cn_delta_r']*delta_r)

        return [eq2, eq6]

    ini_aprox = [0.8*Vsr, 0.]
    V_ini, beta_ini = fsolve(equations_aprox, ini_aprox)

    # Final system (with LG effect)
    # x = [V, beta, K_lg, W]
    def equations(var):
        x = var

        eq2 = Ft_y + 0.5*rho*x[0]**2*adf_data['S']*(adf_data['Cy_beta']*x[1] + adf_data['Cy_delta_r']*delta_r) + x[2]*x[1]*sum(x[3+i] for i in range(n_lg))
        eq3 = Ft_z + adf_data['m']*g + 0.5*rho*x[0]**2*adf_data['S']*adf_data['Cz_0'] + sum(x[3+i] for i in range(n_lg))
        eq4 = Lt + 0.5*rho*x[0]**2*adf_data['S']*adf_data['b']*(adf_data['Cl_beta']*x[1] + adf_data['Cl_delta_r']*delta_r) + sum(x[3+i]*adf_data['lg_y'][i] for i in range(n_lg)) - x[2]*x[1]*sum(x[3+i]*adf_data['lg_z'][i] for i in range(n_lg))
        eq5 = Mt + 0.5*rho*x[0]**2*adf_data['S']*adf_data['c']*adf_data['Cm_0'] - sum(x[3+i]*adf_data['lg_x'][i] for i in range(n_lg)) + adf_data['mu_gr']*sum(x[3+i]*adf_data['lg_z'][i] for i in range(n_lg))
        eq6 = Nt + 0.5*rho*x[0]**2*adf_data['S']*adf_data['b']*(adf_data['Cn_beta']*x[1] + adf_data['Cn_delta_r']*delta_r) - adf_data['mu_gr']*sum(x[3+i]*adf_data['lg_y'][i] for i in range(n_lg)) + x[2]*x[1]*sum(x[3+i]*adf_data['lg_x'][i] for i in range(n_lg))
        eq7 = -x[2] + 2369.074348*abs(x[1])**5 - 1268.0219021*abs(x[1])**4 + 24.6344543*abs(x[1])**3 + 3.1599044*abs(x[1])**2 + 5.2489678*abs(x[1]) + 0.0086458

        return [eq2, eq3, eq4, eq5, eq6, eq7]

    ini    = np.empty(3 + n_lg)
    ini[:] = np.nan

    ini[0]  = V_ini
    ini[1]  = beta_ini
    ini[2]  = 2369.074348*abs(beta_ini)**5 - 1268.0219021*abs(beta_ini)**4 + 24.6344543*abs(beta_ini)**3 + 3.1599044*abs(beta_ini)**2 + 5.2489678*abs(beta_ini) + 0.0086458
    ini[3:] = [adf_data['m']*g/n_lg for i in range(n_lg)]
    
    x = fsolve(equations, ini)
    
    V = x[0]/kt_to_ms
    beta = np.degrees(x[1])
    K_lg = x[2]
    W = x[3:]

    return V

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
            
            # Modes
            f.write('Modes:\n')  
            f.write('{:10.4f}'.format(out_data['lambda_ph'][0,0]) + '{:+8.4f}'.format(out_data['lambda_ph'][0,1]) + 'j\n')
            f.write('{:10.4f}'.format(out_data['lambda_ph'][1,0]) + '{:+8.4f}'.format(out_data['lambda_ph'][1,1]) + 'j\n')
            f.write('{:10.4f}'.format(out_data['lambda_sp'][0,0]) + '{:+8.4f}'.format(out_data['lambda_sp'][0,1]) + 'j\n')
            f.write('{:10.4f}'.format(out_data['lambda_sp'][1,0]) + '{:+8.4f}'.format(out_data['lambda_sp'][1,1]) + 'j\n')
            f.write('\n')
            
            # Phugoid
            f.write('Phugoid:\n')
            f.write('{:>11s}'.format('wn :') + '{:8.5f}'.format(out_data['wn_ph']) + '\n')
            f.write('{:>11s}'.format('zeta :') + '{:8.5f}'.format(out_data['zeta_ph']) + '\n')
            f.write('\n')
            
            # Short period
            f.write('Short period:\n')
            f.write('{:>11s}'.format('wn :') + '{:8.5f}'.format(out_data['wn_sp']) + '\n')
            f.write('{:>11s}'.format('zeta :') + '{:8.5f}'.format(out_data['zeta_sp']) + '\n')
            f.write('\n')
            
            f.write('\n')
            
        # DYNAMIC - LATERAL/DIRECTIONAL 
        if input_data['type'] == 4 or input_data['type'] == 0:
            f.write(lines)
            f.write('{:^80s}'.format('DYNAMIC - LATERAL/DIRECTIONAL') + '\n')
            f.write(lines)
            
            # Modes
            f.write('Modes:\n')  
            f.write('{:10.4f}'.format(out_data['lambda_rs'][0]) + '{:+8.4f}'.format(out_data['lambda_rs'][1]) + 'j\n')
            f.write('{:10.4f}'.format(out_data['lambda_spi'][0]) + '{:+8.4f}'.format(out_data['lambda_spi'][1]) + 'j\n')
            f.write('{:10.4f}'.format(out_data['lambda_dr'][0,0]) + '{:+8.4f}'.format(out_data['lambda_dr'][0,1]) + 'j\n')
            f.write('{:10.4f}'.format(out_data['lambda_dr'][1,0]) + '{:+8.4f}'.format(out_data['lambda_dr'][1,1]) + 'j\n')
            f.write('\n')
            
            # Roll subsidence
            f.write('Roll subsidence:\n')
            f.write('{:>11s}'.format('wn :') + '{:8.5f}'.format(out_data['wn_rs']) + '\n')
            f.write('\n')
            
            # Spiral
            f.write('Spiral:\n')
            f.write('{:>11s}'.format('wn :') + '{:8.5f}'.format(out_data['wn_spi']) + '\n')
            f.write('\n')
            
            # Dutch roll
            f.write('Dutch roll:\n')
            f.write('{:>11s}'.format('wn :') + '{:8.5f}'.format(out_data['wn_dr']) + '\n')
            f.write('{:>11s}'.format('zeta :') + '{:8.5f}'.format(out_data['zeta_dr']) + '\n')
            f.write('\n')
            
            f.write('\n')

        # MINIMUM SPEEDS
        if input_data['type'] == 5 or input_data['type'] == 0:
            f.write(lines)
            f.write('{:^80s}'.format('MINIMUM SPEEDS') + '\n')
            f.write(lines)

            # Vsr
            f.write('{:>11s}'.format('Vsr :') + '{:8.3f}'.format(out_data['Vsr']/kt_to_ms) + ' kt\n')
            f.write('\n')

            # VMCG
            f.write('{:>11s}'.format('VMCG :') + '{:8.3f}'.format(out_data['VMCG']) + ' kt\n')
            f.write('\n')
            
            # VMCA
                                  
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
file = './study/LS.dat'

################################################################################

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

# Minimum control speeds
if input_data['type'] == 5 or input_data['type'] == 0: 
    out_data['Vsr']  = get_stall_speed(adf_data) 
    out_data['VMCG'] = calculate_vmcg(adf_data)
    #out_data['VMCA'] = calculate_vmca(adf_data)
      
# PRINT RESULTS   
write_output_file(out_data, input_data)

"""
# Plot results
if input_data['type'] == 1 or input_data['type'] == 2: # Static    
    plot_out_static(out_data)
else: # Dynamic
    print('Dynamic stability and control not available for plotting')
"""

"""
   FLIGHT MECHANICS
"""

import math
import numpy as np


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
            input_data[title] = variable_data[index]
        elif title == 'type':
            input_data[title] = np.array([int(var) for var in variable_data[index]])
        else:
            input_data[title] = np.array([float(var) for var in variable_data[index]])
                
    adf_data = read_adf_file(input_data['aircraft'][0])
    
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
        adf_data[title] = np.array(variable_data[index])
    
    return adf_data

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
        
    rho = rho_0 * ((1 - (alpha_T * altitude)/T_0) ** ((g/(Ra * alpha_T)) - 1))

    return rho

def get_dynamic_pressure(rho,V):
    """
    This function calculates the dynamic pressure for a given density and speed
    
    Inputs:
        rho: Air density (kg/m3)
        V : Speed (m/s)
        
    Outputs:
        p_d: Dynamic pressure (N/m2)
    
    """
    
    p_d = 0.5 * rho * V ** 2
    
    return p_d

def calculate_static_long(S, c, Ix, Iy, Iz, Jxz, p, q, r, pp, qq, rr,
                          T, epsilon, ni, d_CG_x, d_CG_y, d_CG_z,
                          altitude, V, cm_0, cm_alpha, cm_delta_e, alpha):
    """
    This function calculates the static longitudinal stability and control
    
    Inputs:
        S : Wing surface (m2)
        c : Mean aerodynamic chord (m)
        Ix, Iy, Iz, Jxz : Tensor of inertia (kg*m2)
        p, q, r : Angular velocity (rad/s)
        pp, qq, rr : Angular acceleration (rad/s2)
        T : Engine thrust (N)
        epsilon : Engine angle of atack (deg)
        ni : Engine sideslip angle (deg)
        d_CG_x : X distance from the engine to the CG (m)
        d_CG_y : Y distance from the engine to the CG (m)
        d_CG_z : Z distance from the engine to the CG (m)
        altitude : Flight altitude (ft)
        V : Speed (kt)
        cm_0, cm_alpha, cm_delta_e : longitudinal derivatives (-)
        alpha: Angle of atack (deg)
        
    Output:
        delta_e : Elevator deflection (deg)
    
    """  
    
    # Cinematic momentum in Y axis
    M = Iy * qq - (Iz - Ix) * p * r + Jxz * (p ** 2 - r ** 2)
    
    # Thrust momentum
    n_engines = len(T)
    Mt = 0
    for i in range(n_engines):
        Mt = Mt + T[i] * (math.cos(deg_to_rad * epsilon[i]) * math.cos(deg_to_rad * ni[i]) * d_CG_z[i] +
                          math.sin(deg_to_rad * epsilon[i]) * d_CG_x[i])
    
    # Calculate elevator deflection for balance
    altitude_m = feet_to_meters * altitude
    V_ms = kt_to_ms * V 
    
    rho = get_density(altitude_m)
    p_d = get_dynamic_pressure(rho, V_ms)
    
    delta_e = ((M - Mt) / (p_d * S * c) - cm_0 - cm_alpha * (alpha * deg_to_rad)) / cm_delta_e
    
    delta_e = delta_e / deg_to_rad

    return delta_e

def calculate_static_latdir(S, b, Ix, Iy, Iz, Jxz, p, q, r, pp, qq, rr,
                            T, epsilon, ni, d_CG_x, d_CG_y, d_CG_z, altitude, V,
                            cl_0, cl_beta, cl_delta_a, cl_delta_r,
                            cn_0, cn_beta, cn_delta_a, cn_delta_r,
                            beta):
    """
    This function calculates the static lateral-directionsl stability 
    and control
    
    Inputs:
        S : Wing surface (m2)
        b : Wingspan (m)
        Ix, Iy, Iz, Jxz : Tensor of inertia (kg*m2)
        p, q, r : Angular velocity (rad/s)
        pp, qq, rr : Angular acceleration (rad/s2)
        T : Engine thrust (N)
        epsilon : Engine angle of atack (deg)
        ni : Engine sideslip angle (deg)
        d_CG_x : X distance from the engine to the CG (m)
        d_CG_y : Y distance from the engine to the CG (m)
        d_CG_z : Z distance from the engine to the CG (m)
        altitude : Flight altitude (ft)
        V : Speed (kt)
        cl_0, cl_beta, cl_delta_a, cl_delta_re : lateral derivatives (-)
        cn_0, cn_beta, cn_delta_a, cn_delta_r : Direcional derivatives (-)
        beta: Sideslip angle (deg)
        
    Output:
        delta_a : Aileron deflection (deg)
        delta_r : Rudder deflection (deg)
    
    """  
    
    altitude_m = feet_to_meters * altitude
    V_ms = kt_to_ms * V 
    
    rho = get_density(altitude_m)
    p_d = get_dynamic_pressure(rho, V_ms)
    
    # Cinematic momentum in X axis
    L = Ix * pp - Jxz * rr + (Iz - Iy) * q * r - Jxz * p * q
    
    # Thrust momentum in X asis
    n_engines = len(T)
    Lt = 0
    for i in range(n_engines):
        Lt = Lt + T[i] * (- math.cos(deg_to_rad * epsilon[i]) * math.sin(deg_to_rad * ni[i]) * d_CG_z[i]
                          - math.sin(deg_to_rad * epsilon[i]) * d_CG_y[i])
    
    # Cinematic momentum in Z axis
    N = Iz * rr - Jxz * pp - (Ix - Iy) * p * q + Jxz * q * r
    
    # Thrust momentum in Z asis
    n_engines = len(T)
    Nt = 0
    for i in range(n_engines):
        Nt = Nt + T[i] * (- math.cos(deg_to_rad * epsilon[i]) * math.cos(deg_to_rad * ni[i]) * d_CG_y[i]
                          + math.cos(deg_to_rad * epsilon[i]) * math.sin(deg_to_rad * ni[i]) * d_CG_x[i])
    
    # Calculate ailerond and rudder deflection for balance (Eq. system: A*X=B)
    A = np.array([[cl_delta_a, cl_delta_r], [cn_delta_a, cn_delta_r]])
    B = np.array([[(L - Lt) / (p_d * S * b) - cl_0 - cl_beta * (deg_to_rad * beta)],
                  [(N - Nt) / (p_d * S * b) - cn_0 - cn_beta * (deg_to_rad * beta)]])
    
    X = np.linalg.solve(A[:,:,0],B[:,:,0])
    
    delta_a = X[0] / deg_to_rad
    delta_r = X[1] / deg_to_rad

    return delta_a, delta_r

"""
MAIN PROGRAM

"""

"""
 Global variables:
     feet_to_meters: Feet to meters conversion factor
     kt_to_ms: Knots to meters per second conversion factor
     deg_to_rad: Degrees to radians conversion factor
     g: Acceleration of gravity

"""
feet_to_meters = 0.3048
kt_to_ms = 0.514444
deg_to_rad = math.pi / 180.0
g = 9.80665
    

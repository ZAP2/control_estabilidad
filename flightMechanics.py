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
            input_data[title] = np.array([int(var) for var in variable_data[index]])
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

def plot_out_static(out_data):
    """ 
    This function plot the calculated outputs
    
    Inputs:
        out_data: Dictionary containing the results
        
    Outputs:
        plots
    
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

def output_file(out_data):
    """ 
    this function writes an output file
    ​
    Inputs: 
    out_data: dictionary with the variables already calculated
    ​
    Outputs:
    output_file: file.txt  with the study results
    ​
    """
    variables = ['altitude',     'V',     'p',     'q',     'r',     'pp',     'qq',     'rr', 'alpha', 'beta', 'delta_e', 'delta_a', 'delta_r']
    header    = ['Altitude', 'Speed',     'p',     'q',     'r',     'pp',     'qq',     'rr', 'Alpha', 'Beta', 'Delta_e', 'Delta_a', 'Delta_r']
    units     = [      'ft',    'kt', 'rad/s', 'rad/s', 'rad/s', 'rad/s2', 'rad/s2', 'rad/s2',   'deg',  'deg',     'deg',     'deg',     'deg']

    n_var = len(variables)
    n_row = len(out_data['altitude'])

    with open(out_data['file'], "w") as f:

        f.write("\n")

        f.write("Aircraft: " + out_data["aircraft"]) #escribir 'titulo' indicando el avion que es
        f.write("\n"*3)

        # Headers
        string = ['{:>11s}'.format(s) for s in header]
        string = ''.join(string) + '\n'
        f.write(string) 

        #Units  
        string = ['{:>11s}'.format(s) for s in units]
        string = ''.join(string) + '\n'
        f.write(string)

        #variables
        for i in range(n_var):
            string = [out_data[var][i] for var in variables]
            string = ['{:>11.3f}'.format(s) for s in string]
            string = ''.join(string) + '\n'
            f.write(string)


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

# Define input file path
#file = './longitudinal/static_long.dat'
file = 'input_file_template.dat'

# Read input and aircraft data
input_data, adf_data = read_input_file(file)

# Check alpha and beta input values
if input_data['type'] == 1 and (len(input_data['beta']) > 1 or input_data['beta'][0] != 0):
    print('Longitudinal study: beta should be zero')
    print('Setting beta as 0 deg')
    input_data['beta'] = np.array([0.0])
if input_data['type'] == 2 and len(input_data['alpha']) > 1:
    print('Lateral dirctional study: alpha should be constant, an array is not allowed')
    print('Setting alpha as alpha[0]')
    input_data['alpha'] = np.array([input_data['alpha'][0]])

for key,val in input_data.items():
        exec(key + '=val')
for key,val in adf_data.items():
        exec(key + '=val')

# Calculate stability and control

# Static
if input_data['type'] == 1 or input_data['type'] == 2:
    
    # Initialize output dictionary
    variable_name = ['file']
    for key, val in input_data.items():
        variable_name.append(key)
    variable_name = variable_name + ['delta_e', 'delta_a', 'delta_r']
    out_data = dict([(key, []) for key in variable_name])
    
    out_data['file']     = file[:-4] + '.out'
    out_data['aircraft'] = input_data['aircraft']
    out_data['type']     = input_data['type']
    
    # Calculate    
    for alt in altitude:    # Scan altitude array
        for spd in V:       # Scan speed array
            
            # Static - longitudinal
            if input_data['type'] == 1:
                for alpha_ in alpha: # Scan alpha array
                    delta_e = calculate_static_long(S, c, Ix, Iy, Iz, Jxz, p, q, r, pp, qq, rr,
                                                    T, epsilon, ni, d_CG_x, d_CG_y, d_CG_z,
                                                    alt, spd, cm_0, cm_alpha, cm_delta_e, alpha_)
                    
                    out_data['altitude'] = np.append(out_data['altitude'], alt)
                    out_data['V']        = np.append(out_data['V'], spd)
                    out_data['alpha']    = np.append(out_data['alpha'], alpha_)
                    out_data['beta']     = np.append(out_data['beta'], input_data['beta'][0])
                    out_data['delta_e']  = np.append(out_data['delta_e'], delta_e[0])
                    out_data['delta_a']  = np.append(out_data['delta_a'], 0.0)
                    out_data['delta_r']  = np.append(out_data['delta_r'], 0.0)
                    
            # Static - lateral-directional   
            elif input_data['type'] == 2:
                for beta_ in beta:  # Scan beta array
                    delta_a, delta_r = calculate_static_latdir(S, b, Ix, Iy, Iz, Jxz, p, q, r, pp, qq, rr,
                                                               T, epsilon, ni, d_CG_x, d_CG_y, d_CG_z, alt, spd,
                                                               cl_0, cl_beta, cl_delta_a, cl_delta_r,
                                                               cn_0, cn_beta, cn_delta_a, cn_delta_r,
                                                               beta_)
                    
                    out_data['altitude'] = np.append(out_data['altitude'], alt)
                    out_data['V']        = np.append(out_data['V'], spd)
                    out_data['alpha']    = np.append(out_data['alpha'], input_data['alpha'][0])
                    out_data['beta']     = np.append(out_data['beta'], beta_)
                    out_data['delta_e']  = np.append(out_data['delta_e'], 0.0)
                    out_data['delta_a']  = np.append(out_data['delta_a'], delta_a[0])
                    out_data['delta_r']  = np.append(out_data['delta_r'], delta_r[0])
            
    # Compete dictionary with constant variables            
    n_cases = len(out_data['altitude'])
    out_data['p']  = np.ones(n_cases) * input_data['p']
    out_data['q']  = np.ones(n_cases) * input_data['q']
    out_data['r']  = np.ones(n_cases) * input_data['r']
    out_data['pp'] = np.ones(n_cases) * input_data['pp']
    out_data['qq'] = np.ones(n_cases) * input_data['qq']
    out_data['rr'] = np.ones(n_cases) * input_data['rr']

# Dynamic - TO DO   
else:
    print('Dynamic stability and control not available for calculation')
    
    
# Plot results
if input_data['type'] == 1 or input_data['type'] == 2: # Static    
    plot_out_static(out_data)
else: # Dynamic
    print('Dynamic stability and control not available for plotting')


# Print results (output_file)
if input_data['type'] == 1 or input_data['type'] == 2: # Static    
    output_file(out_data)
else: # Dynamic
    print('Dynamic stability and control not available for printing')

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
from scipy.interpolate import interp1d
from scipy.interpolate import griddata
from scipy.optimize import fsolve
import CoolProp.CoolProp as CP
from scipy.optimize import minimize

# Conversion factor from Barrer to SI units (mol/m^2/s/Pa)
BARRER_TO_SI = 3.3464e-16

# FEED STREAM INPUT DATA FOR THE MEMBRANE in Kmol/hr
# This data is obtained by the Aspen simulations
#The current feed composition
feed_composition = {
    'Hydrogen':137.607,
    'CarbonDioxide': 45.8676,
    'Methanol': 0.747245,
    'Water': 0.147976,
    'CarbonMonoxide': 5.09148,
}
'''
feed_composition = {
    'Hydrogen': 347.877,
    'CarbonDioxide': 10.1493,
    'Methanol': 0.705218,
    'Water': 0.000357253,
    'CarbonMonoxide': 0.00763042,
}
'''
'''
#New feed composition
feed_composition = {
    'Hydrogen': 84.6703,
    'CarbonDioxide': 270.173,
    'Methanol': 1.56457,
    'Water': 0.275612,
    'CarbonMonoxide': 28.5863,
}
'''

# Function to interpolate permeability values based on the data provided.
def interpolate_permeability(x, x_data, y_data):
    return np.interp(x, x_data, y_data)

# Function to calculate the permeablity data for given Temperature, Pressure and mole fraction of the components
def get_permeability_data(temp, pressure, mole_fraction):
    # Data for permeability in Barrer
    # The data is obtained from a research paper and this is for a perticular membrane
    data = {
        100: {
            'Hydrogen': [15, 20, 25],
            'Hydrogen_permeability': [981.23, 981.63, 982.03],
            'CarbonDioxide': [15, 20, 25],
            'CarbonDioxide_permeability': [533.63, 529.55, 527.38],
            'Methanol_activity': [0.0140, 0.0211, 0.0282, 0.0494, 0.0885],
            'Methanol_permeability': [7498.73, 6625.36, 6161.08, 5086.42, 4110.14],
            'Water_activity': [0.0077, 0.0119, 0.0156, 0.0203],
            'Water_permeability': [33049.12, 29433.43, 26703.53, 24944.60],
        },
        150: {
            'Hydrogen': [15, 20, 25],
            'Hydrogen_permeability': [1221.79, 1222.19, 1222.60],
            'CarbonDioxide': [15, 20, 25],
            'CarbonDioxide_permeability': [465.56, 460.54, 459.34],
            'Methanol_activity': [0.0046, 0.0065, 0.0082, 0.0108, 0.0139],
            'Methanol_permeability': [5017.77, 4335.78, 4162.89, 3735.39, 3444.19],
            'Water_activity': [0.0049, 0.0077, 0.0118, 0.0156, 0.0202],
            'Water_permeability': [23967.36, 20791.01, 17417.20, 15574.39, 14541.11],
        },
        200: {
            'Hydrogen': [15, 20, 25],
            'Hydrogen_permeability': [1349.34, 1346.52, 1343.69],
            'CarbonDioxide': [15, 20, 25],
            'CarbonDioxide_permeability': [405.16, 406.85, 405.67],
            'Methanol_activity': [0.0046, 0.0061, 0.0082, 0.0107, 0.0139],
            'Methanol_permeability': [4145.05, 3781.27, 3517.44, 3262.67, 3071.47],
            'Water_activity': [0.0049, 0.0077, 0.0118, 0.0156],
            'Water_permeability': [11104.58, 9339.47, 7941.41, 7227.63],
        },
    }

    if temp not in data:
        raise ValueError("Temperature must be one of [100, 150, 200] °C.")

    temp_data = data[temp]

    # Calculate permeabilities for each component
    hydrogen_perm = interpolate_permeability(pressure, temp_data['Hydrogen'], temp_data['Hydrogen_permeability']) * BARRER_TO_SI
    co2_perm = interpolate_permeability(pressure, temp_data['CarbonDioxide'], temp_data['CarbonDioxide_permeability']) * BARRER_TO_SI
    methanol_perm = interpolate_permeability(
        mole_fraction['Methanol'], temp_data['Methanol_activity'], temp_data['Methanol_permeability']
    ) * BARRER_TO_SI
    water_perm = interpolate_permeability(
        mole_fraction['Water'], temp_data['Water_activity'], temp_data['Water_permeability']
    ) * BARRER_TO_SI
    # Carbon monoxide permeability is assumed to be the same as CO2
    co_perm = co2_perm

    # Construct the output dictionary
    permeability_data = {
        'Hydrogen': hydrogen_perm,
        'CarbonDioxide': co2_perm,
        'Methanol': methanol_perm,
        'Water': water_perm,
        'CarbonMonoxide': co_perm,
    }

    return permeability_data

# This function simulates membrane using Real gas conditions and using fugacity
def simulate_membrane_module_partial_pressure_debug_fugacity(
    feed_composition, feed_pressure, feed_temperature, permeability, membrane_area, stage_cut, membrane_thickness, length, T_D, N_fiber
):
    # Constants
    # Constants and initialization
    R = 8.314  # J/(mol·K)
    T = feed_temperature + 273.15  # Convert Celsius to Kelvin
    n_points = 1000  # Number of discretized steps
    dx = 1 / n_points  # Step size for discretization
    dL = length / n_points # Step size for length
    
    flux_Hydrogen = []
    flux_CO2 = []
    flux_methanol = []
    flux_water = []
    flux_co = []

    # Initialize pressures and flows
    feed_side_pressure = feed_pressure  # bar
    permeate_side_pressure = feed_pressure * (1 - stage_cut)  # bar
    total_feed_flow = sum(feed_composition.values())  # kmol/hr

    # Normalize feed composition to mole fractions
    feed_fraction = {comp: flow / total_feed_flow for comp, flow in feed_composition.items()}

    # Initialize flow rates and compositions
    retentate_flow = total_feed_flow
    permeate_flow = 0
    retentate_composition = feed_fraction.copy()
    permeate_composition = {comp: 0 for comp in feed_composition}
    
     # Data storage for plotting mole fractions
    mole_fraction_retentate = {comp: [] for comp in feed_composition}
    mole_fraction_permeate = {comp: [] for comp in feed_composition}
    
    # Discretized simulation loop
    for i in range(n_points):
        current_flow = retentate_flow
        current_fraction = retentate_composition
        
        permeability = get_permeability_data(feed_temperature, feed_side_pressure, current_fraction)
        feed_side_pressure_drop =  calculate_pressure_drop(current_fraction, feed_side_pressure, T, T_D, length, N_fiber)
        feed_side_pressure = feed_side_pressure-feed_side_pressure_drop
        
        permeation = {}
        for comp, y in current_fraction.items():
            # Partial pressure driving force
            p_feed = feed_side_pressure * y  # Partial pressure on feed side (bar)
            z_feed = PropsSI('Z', 'P', feed_side_pressure * 1e5, 'T', T, comp)
            fugacity_feed = z_feed * p_feed * 1e5  # Pa
            p_perm = permeate_side_pressure * permeate_composition.get(comp, 0)  # Partial pressure on permeate side (bar)
            z_perm = PropsSI('Z', 'P', permeate_side_pressure * 1e5, 'T', T, comp)
            fugacity_perm = z_perm * p_perm * 1e5  # Pa
            driving_force = max(fugacity_feed - fugacity_perm, 0)  # Ensure non-negative driving force

            # Permeation flux calculation
            perm = permeability.get(comp, 0)  # Permeability (m³(STP)/m²/s/Pa)
            flux = (perm / membrane_thickness) * driving_force #df in pa
            permeation[comp] = flux * membrane_area * dx  # Molar permeation rate (kmol/hr)
            

            if comp == 'Hydrogen':
                flux_Hydrogen.append(permeation['Hydrogen'])
            elif comp == 'CarbonDioxide':
                flux_CO2.append(permeation['CarbonDioxide'])
            elif comp == 'CarbonMonoxide':
                flux_co.append(permeation['CarbonMonoxide'])
            elif comp == 'Methanol':
                flux_methanol.append(permeation['Methanol'])
            elif comp == 'Water':
                flux_water.append(permeation['Water'])
            
        # Total permeation flux
        total_permeation = sum(permeation.values())
        if total_permeation > current_flow:
            # Cap total permeation at the current retentate flow
            scale_factor = current_flow / total_permeation
            permeation = {comp: rate * scale_factor for comp, rate in permeation.items()}
            total_permeation = current_flow
        
        
        # Update flows
        permeate_flow += total_permeation
        retentate_flow -= total_permeation

        # Update permeate composition
        for comp, rate in permeation.items():
            permeate_composition[comp] += rate / permeate_flow if permeate_flow > 0 else 0

        # Normalize permeate composition to ensure mole fractions
        if permeate_flow > 0:
            total_permeate_fraction = sum(permeate_composition.values())
            permeate_composition = {comp: x / total_permeate_fraction for comp, x in permeate_composition.items()}

        # Update retentate composition
        for comp in current_fraction:
            retentate_composition[comp] -= permeation[comp] / retentate_flow if retentate_flow > 0 else 0
        # **Normalize retentate composition to ensure mole fractions sum to 1**
        if retentate_flow > 0:
            total_retentate_fraction = sum(retentate_composition.values())
            retentate_composition = {comp: x / total_retentate_fraction for comp, x in retentate_composition.items()}
        # Store mole fractions for plotting
        for comp in feed_composition:
            mole_fraction_retentate[comp].append(retentate_composition[comp])
            mole_fraction_permeate[comp].append(permeate_composition[comp])
    # Final outputs
    return {
        'feed_flow_rate': total_feed_flow,
        'retentate_flow_rate': retentate_flow, #kmol/hr
        'permeate_flow_rate': permeate_flow,
        'retentate_composition': retentate_composition, #mole fraction 
        'permeate_composition': permeate_composition,
        'retentate_pressure': feed_side_pressure, # bar
        'permeate_pressure': permeate_side_pressure,
    }

# It is assumed that the membrane module is a hollow fiber module and the feed is passed through the hollow fiber
# This function is used get the Fiber length and the fiber Diameter 
def calculate_hollow_fiber_dimensions(mem_area, num_fibers, fiber_length=None, fiber_diameter=None):

    if fiber_length and fiber_diameter:
        raise ValueError("Provide either fiber_length or fiber_diameter, not both.")

    if fiber_length:
        # Solve for fiber diameter
        fiber_diameter = mem_area / (num_fibers * math.pi * fiber_length)
    
    elif fiber_diameter:
        # Solve for fiber length
        fiber_length = mem_area / (num_fibers * math.pi * fiber_diameter)

    else:
        raise ValueError("Either fiber_length or fiber_diameter must be provided.")

    return {
        "Fiber Length (m)": fiber_length,
        "Fiber Diameter (mm)": fiber_diameter*1000,
        "Fiber Diameter (m)": fiber_diameter
    }

# This function calculates the pressure drop if the feed/Retentate which through the hollow fiber tube
# It is Assumed and also later observed that the Hagen–Poiseuille equation is used and there is a laminar flow in the tube
# Non Ideal gas conditions are followed to calculate the parameters required in the calculations
def calculate_pressure_drop(feed_composition, pressure, temperature, tube_diameter, tube_length, num_fibers):
    # Universal gas constant in J/(mol*K)
    R = 8.314  

    # Convert pressure from bar to Pa
    pressure_pa = pressure * 1e5  

    # Convert temperature to Kelvin
    temperature_k = temperature + 273.15  

    # Convert feed composition from kmol/hr to mol/s
    total_molar_flow = sum(feed_composition.values()) * 1000 / 3600  # Convert to mol/s

    # Calculate total volumetric flow rate using ideal gas law (Q = nRT/P)
    volumetric_flow_rate = (total_molar_flow * R * temperature_k) / pressure_pa  # m³/s

    # Calculate velocity (v = Q / A)
    tube_cross_section = (math.pi * (tube_diameter**2)) / 4  # m²
    if num_fibers * tube_cross_section == 0:
        raise ValueError("Invalid tube dimensions; denominator is zero.")

    velocity = volumetric_flow_rate / (num_fibers * tube_cross_section)  # m/s

    # Define gas components recognized by CoolProp
    coolprop_components = {
        'Hydrogen': 'H2',
        'CarbonDioxide': 'CO2',
        'Methanol': 'Methanol',
        'Water': 'Water',
        'CarbonMonoxide': 'CO'
    }

    # Normalize feed composition to mole fractions
    total_moles = sum(feed_composition.values())
    mole_fractions = {comp: mol / total_moles for comp, mol in feed_composition.items()}

    # Calculate mixture density using mole-fraction-weighted averages
    density = sum(
        mole_fractions[comp] * CP.PropsSI('D', 'T', temperature_k, 'P', pressure_pa, coolprop_components[comp])
        for comp in feed_composition
    )

    # Viscosity calculation using Wilke's mixing rule
    viscosities = {}
    for comp in feed_composition:
        try:
            viscosities[comp] = CP.PropsSI('V', 'T', temperature_k, 'P', pressure_pa, coolprop_components[comp])
        except ValueError:
            #print(f"Warning: No viscosity data for {comp}. Using estimated values.")
            viscosities[comp] = estimate_viscosity(comp, temperature_k)  # Use empirical estimate

    # Apply Wilke’s mixing rule
    viscosity = wilkes_mixing_rule(viscosities, mole_fractions)

    # Calculate Reynolds number (Re = ρvD / μ)
    reynolds_number = (density * velocity * tube_diameter) / viscosity
    
    # Print intermediate values for debugging

    # Calculate friction factor
    if reynolds_number < 2000:
        # Laminar flow: f = 64/Re
        friction_factor = 64 / max(reynolds_number, 1e-5)  # Avoid division by zero
    elif reynolds_number > 4000:
        # Turbulent flow: Use Haaland approximation for stability
        friction_factor = 1 / (1.8 * math.log10(6.9 / reynolds_number))**2
    else:
        # Transitional flow: Estimate using weighted average
        laminar_f = 64 / reynolds_number
        turbulent_f = 1 / (1.8 * math.log10(6.9 / reynolds_number))**2
        friction_factor = (laminar_f + turbulent_f) / 2  

    # Calculate pressure drop using Darcy-Weisbach equation
    pressure_drop_pa = (friction_factor * tube_length / tube_diameter) * (density * velocity**2 / 2)  # Pa
    pressure_drop_bar = pressure_drop_pa / 1e5  # Convert to bar

    return pressure_drop_bar

def estimate_viscosity(component, temperature_k):
    """Returns estimated viscosity (Pa.s) for missing gases using empirical correlations."""
    viscosity_data = {
        'CarbonMonoxide': 17.2e-6,  # Approximate viscosity (Pa.s) at 300K
        'Methanol': 0.000544        # Approximate viscosity (Pa.s) at 300K
    }
    return viscosity_data.get(component, 1e-5)  # Default fallback

def wilkes_mixing_rule(viscosities, mole_fractions):
    """Applies Wilke's mixing rule for gas viscosity."""
    mu_mix = 0
    for i in mole_fractions:
        sum_x_mu = sum(
            (mole_fractions[j] * (1 + math.sqrt(viscosities[i] / viscosities[j])) ** 2) / 
            (8 * math.sqrt(1 + viscosities[i] / viscosities[j]))
            for j in mole_fractions
        )
        mu_mix += mole_fractions[i] * viscosities[i] / sum_x_mu
    return mu_mix

'''
INPUT SECTION
_______________________________________________________________________________________________________________________
'''
feed_pressure = 65  # bar
feed_temperature = 200  # Celsius
membrane_ar = 103  # m^2
stage_cut = 0.9
num_tubes = 10000   # Number of fibers
fiber_length = 3 # Assume fiber length in meters

'''
CALCULATED INPUT SECTION
________________________________________________________________________________________________________________________
'''
result = calculate_hollow_fiber_dimensions(membrane_ar, num_tubes, fiber_length=fiber_length)
FIBER_LENGTH_M = result['Fiber Length (m)']
FIBER_DIAMETER_M = result['Fiber Diameter (m)'] 
# Calculate the total molar flow rate
total_molar_flow_rate = sum(feed_composition.values())
print(total_molar_flow_rate)
# Calculate the mole fractions
mole_fractions = {component: flow_rate / total_molar_flow_rate for component, flow_rate in feed_composition.items()}

# Get permeability data
permeability = get_permeability_data(feed_temperature, feed_pressure, mole_fractions)

'''
MAIN MEMBRANE SIMULATING SECTION
________________________________________________________________________________________________________________________
'''
# Here in simulating the membrane the thickness of the membrane is assumed to be 1 micron
results_DF_PP_f = simulate_membrane_module_partial_pressure_debug_fugacity(
    feed_composition, feed_pressure, feed_temperature, permeability, membrane_ar, stage_cut, 1e-6, FIBER_LENGTH_M, FIBER_DIAMETER_M, num_tubes
)

print(results_DF_PP_f)


'''
MEMBRANE OPTIMIZATION SECTION USING RESPONSE SURFACE METHEDOLOGY
________________________________________________________________________________________________________________________
'''
def perform_rsm_analysis_2(feed_composition, feed_pressure, feed_temperature, permeability, design_points, function, FIBER_LENGTH_M, FIBER_DIAMETER_M, num_tubes):
    results = []

    for membrane_area, stage_cut in design_points:
        # Run the simulation function
        simulation_output = function(
            feed_composition, feed_pressure, feed_temperature, 
            permeability, membrane_area, stage_cut, 1e-6, FIBER_LENGTH_M, FIBER_DIAMETER_M, num_tubes
        )
        
        # Calculate %recovery
        permeate_flow_rate = simulation_output['permeate_flow_rate']
        perm_compo = simulation_output['permeate_composition']
        retentate_flow_rate = simulation_output['retentate_flow_rate']
        ret_comp = simulation_output['retentate_composition']
        Perm_comp_molar_F = {key: value * permeate_flow_rate for key, value in perm_compo.items()}
        ret_comp_molar_F = {key: value * retentate_flow_rate for key, value in ret_comp.items()}
        feed_flow_rate = simulation_output['feed_flow_rate']
        percent_recovery = (permeate_flow_rate / feed_flow_rate) * 100
        percentage_recovery_Hydrogen =  (Perm_comp_molar_F['Hydrogen']/feed_composition['Hydrogen'])*100
        percentage_loss_Hydrogen =  (ret_comp_molar_F['Hydrogen']/feed_composition['Hydrogen'])*100

        # Append results
        results.append({
            'Membrane Area': membrane_area,
            'Stage Cut': stage_cut,
            '%Recovery_H2': percentage_recovery_Hydrogen,
            '%Loss_H2':percentage_loss_Hydrogen,
            'ret_flowrate':retentate_flow_rate
        })
    
    # Convert results to a DataFrame
    results_df = pd.DataFrame(results)
    return results_df

# Define the objective function to minimize hydrogen loss
def objective_function(params, x_col, y_col, z_col, data):
    membrane_area = params[0]
    stage_cut = params[1]
    
    # Extract data columns
    x = data[x_col].values
    y = data[y_col].values
    z = data[z_col].values  # Hydrogen loss column

    # Interpolate the loss at (membrane_area, stage_cut)
    loss = griddata((x, y), z, (membrane_area, stage_cut), method='cubic')
    
    # If interpolation is NaN (out of bounds), return a large number
    if np.isnan(loss):
        return np.inf

    return loss  # We minimize this function

# Function to plot contour and find minimum loss of hydrogen
def plot_contour_and_find_minimum_loss(data, x_col, y_col, z_col):
    # Initial guess (middle of the design space)
    initial_guess = [np.mean(data[x_col].values), np.mean(data[y_col].values)]

    # Bounds for optimization
    bounds = [(min(data[x_col].values), max(data[x_col].values)), 
              (min(data[y_col].values), max(data[y_col].values))]

    # Perform minimization
    result = minimize(objective_function, initial_guess, args=(x_col, y_col, z_col, data), bounds=bounds, method='L-BFGS-B')

    # Extract the optimal values
    optimal_membrane_area, optimal_stage_cut = result.x
    optimal_loss = result.fun
    
    print(f"Optimal Membrane Area: {optimal_membrane_area}")
    print(f"Optimal Stage Cut: {optimal_stage_cut}")
    print(f"Minimum Hydrogen Loss: {optimal_loss}%")

    # Create grid for contour plot
    x = data[x_col].values
    y = data[y_col].values
    z = data[z_col].values
    xi = np.linspace(min(x), max(x), 100)
    yi = np.linspace(min(y), max(y), 100)
    xi, yi = np.meshgrid(xi, yi)
    zi = griddata((x, y), z, (xi, yi), method='cubic')

    # Plot contour
    plt.figure(figsize=(8, 6))
    contour = plt.contourf(xi, yi, zi, levels=20, cmap='coolwarm')
    plt.colorbar(contour, label=z_col)
    plt.scatter(x, y, c='black', label='Design Points')
    #plt.scatter(optimal_membrane_area, optimal_stage_cut, c='blue', marker='x', label='Optimal Point', s=100)
    plt.title(f'Contour Plot of {z_col}')
    plt.xlabel(x_col)
    plt.ylabel(y_col)
    plt.legend()
    plt.show()

    return {
        'Membrane Area': optimal_membrane_area,
        'Stage Cut': optimal_stage_cut,
        'Min %Loss_H2': optimal_loss
    }

# NEW FUNCTION FOR OPTIMIZATION
def find_optimal_point_and_plot(data, x_col, y_col, z_col, ret_flowrate_col, area_weight=0, stage_cut_weight=0.1):
    # Filter out rows where retentate flow rate is zero
    valid_data = data[data[ret_flowrate_col] > 0].copy()

    # Define the objective function
    def objective_function(params):
        membrane_area, stage_cut = params
        # Interpolate the hydrogen loss for the given membrane area and stage cut
        points = valid_data[[x_col, y_col]].values
        values = valid_data[z_col].values
        loss = griddata(points, values, (membrane_area, stage_cut), method='linear')
        # Add penalties for larger membrane area and higher stage cut
        penalty = area_weight * membrane_area + stage_cut_weight * stage_cut
        return loss + penalty

    # Initial guess (middle of the design space)
    initial_guess = [np.mean(valid_data[x_col].values), np.mean(valid_data[y_col].values)]

    # Bounds for optimization
    bounds = [(min(valid_data[x_col].values), max(valid_data[x_col].values)), 
              (min(valid_data[y_col].values), max(valid_data[y_col].values))]

    # Perform minimization
    result = minimize(objective_function, initial_guess, bounds=bounds, method='L-BFGS-B')

    # Extract the optimal values
    optimal_membrane_area, optimal_stage_cut = result.x
    optimal_loss = result.fun

    # Clamp the optimal loss to the minimum and maximum values in the valid data
    optimal_loss = np.clip(optimal_loss, np.min(valid_data[z_col].values), np.max(valid_data[z_col].values))

    # Print the optimal values
    print(f"Optimal Membrane Area: {optimal_membrane_area}")
    print(f"Optimal Stage Cut: {optimal_stage_cut}")
    print(f"Minimum Hydrogen Loss: {optimal_loss}%")

    # Create grid for contour plot
    x = valid_data[x_col].values
    y = valid_data[y_col].values
    z = valid_data[z_col].values
    xi = np.linspace(min(x), max(x), 100)
    yi = np.linspace(min(y), max(y), 100)
    xi, yi = np.meshgrid(xi, yi)
    zi = griddata((x, y), z, (xi, yi), method='linear')

    # Clamp the interpolated values to the minimum and maximum values in the valid data
    zi = np.clip(zi, np.min(z), np.max(z))

    # Plot contour
    plt.figure(figsize=(10, 8))
    contour = plt.contourf(xi, yi, zi, levels=20, cmap='coolwarm')
    
    # Increase font size for colorbar label
    cbar = plt.colorbar(contour)
    cbar.set_label('Hydrogen Loss (%)', fontsize=18)
    
    plt.scatter(x, y, c='black', label='Design Points')
    # plt.scatter(optimal_membrane_area, optimal_stage_cut, c='blue', marker='x', label='Optimal Point', s=100)
    
    # Increase font sizes
    plt.title('Contour Plot of Hydrogen Loss', fontsize=18)
    plt.xlabel('Membrane Area (m²)', fontsize=18)
    plt.ylabel('Stage Cut', fontsize=18)
    
    plt.legend(fontsize=18)  # Increase legend font size
    plt.show()

    return {
        'Membrane Area': optimal_membrane_area,
        'Stage Cut': optimal_stage_cut,
        'Min %Loss_H2': optimal_loss
    }

'''
RSM INPUT SECTION
_______________________________________________________________________________________________________________________
'''
# Define ranges and levels
membrane_area_range = (5, 120) #Currently it says the range is in between 5m2 to 2000m2
stage_cut_range = (0.1, 0.9)

num_levels_membrane_area = 8  # Number of points for membrane area
num_levels_stage_cut = 5      # Number of points for stage cut

# Generate levels
membrane_area_points = np.linspace(*membrane_area_range, num_levels_membrane_area)
stage_cut_points = np.linspace(*stage_cut_range, num_levels_stage_cut)

# Create a full factorial design
design_points_700 = [(area, cut) for area in membrane_area_points for cut in stage_cut_points]

'''
RSM OUTPUT SECTION
________________________________________________________________________________________________________________________
'''

results_df_2 = perform_rsm_analysis_2(
    feed_composition, feed_pressure, feed_temperature, 
    permeability, design_points_700, simulate_membrane_module_partial_pressure_debug_fugacity, FIBER_LENGTH_M, FIBER_DIAMETER_M, num_tubes
)

'''
# Run the function to find the minimum hydrogen loss
'''


# NEW Optimization
# Find the optimal point
result_optmum = find_optimal_point_and_plot(results_df_2, 'Membrane Area', 'Stage Cut', '%Loss_H2', 'ret_flowrate', area_weight=0.1, stage_cut_weight=0.1)
print(result_optmum)

'''
OUTPUT SECTION
________________________________________________________________________________________________________________________
'''

# Extracting the Relevant outputs
Feed_comp_Real = results_DF_PP_f['feed_flow_rate']
retentate_flow_rate_real = results_DF_PP_f['retentate_flow_rate']
permeate_flow_rate_real = results_DF_PP_f['permeate_flow_rate']
retentate_composition_real = results_DF_PP_f['retentate_composition']
permeate_composition_real = results_DF_PP_f['permeate_composition']
Retentate_pressure = results_DF_PP_f['retentate_pressure']
Permeate_pressure = results_DF_PP_f['permeate_pressure']

# Printing the outputs
print('Membrane specs:')
print('Fiber Length in meter')
print(FIBER_LENGTH_M)
print('Fiber Inside Diameter in mm:')
print(FIBER_DIAMETER_M*1000)
print('Membrane Outlet:')
print('Feed Flow Rate in Kmol/hr:')
print(Feed_comp_Real)
print('Retentate Flow Rate in Kmol/hr:')
print(retentate_flow_rate_real)
print('Permeate Flow Rate in Kmol/hr:')
print(permeate_flow_rate_real)
print('Retentate Composition in mole fraction')
print(retentate_composition_real)    
print('Permeate Composition in mole fraction')
print(permeate_composition_real)
print('Retentate Pressure in Bar')
print(Retentate_pressure)
print('Permeate Pressure in Bar:')
print(Permeate_pressure)


#Printing the data frame from the RSM analysis
print('RSM ANALYSIS DATA FRAME:')
print(results_df_2)

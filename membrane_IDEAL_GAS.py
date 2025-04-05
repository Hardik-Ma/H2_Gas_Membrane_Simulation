import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
from scipy.interpolate import interp1d

# Conversion factor from Barrer to SI units (mol/m^2/s/Pa)
BARRER_TO_SI = 3.3464e-16

# FEED STREAM INPUT DATA FOR THE MEMBRANE
feed_composition = {
    'Hydrogen': 111.971,
    'CarbonDioxide': 181.941,
    'Methanol': 1.04958,
    'Water': 0.201152,
    'CarbonMonoxide': 15.8056,
}

def interpolate_permeability(x, x_data, y_data):
    """
    Function to interpolate permeability values based on the data provided.
    """
    return np.interp(x, x_data, y_data)

def get_permeability_data(temp, pressure, mole_fraction):
    """
    Function to calculate the permeability data for given feed temperature,
    feed pressure, and mole fractions of components.
    
    Parameters:
    temp (int): Temperature in Celsius (100, 150, or 200).
    pressure (float): Pressure in bars.
    mole_fraction (dict): Mole fraction of components (keys: 'Hydrogen', 'CarbonDioxide', 'Methanol', 'Water', 'CarbonMonoxide').

    Returns:
    dict: Permeability data in SI units (mol/m^2/s/Pa).
    """
    # Data for permeability in Barrer
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

def simulate_membrane_module_partial_pressure_debug(
    feed_composition, feed_pressure, feed_temperature, permeability, membrane_area, stage_cut, membrane_thickness
):
    """
    Simulates a membrane module with simplified driving force calculations using partial pressure.

    Parameters:
    - feed_composition: dict of component molar flow rates in kmol/hr.
    - feed_pressure: feed-side pressure in bar.
    - feed_temperature: feed temperature in Celsius.
    - permeability: dict of permeability values in SI units (m³(STP)/m²/s/Pa).
    - membrane_area: membrane area in m².
    - stage_cut: ratio of permeate flow to feed flow.
    - membrane_thickness: thickness of the membrane in meters.

    Returns:
    - dict containing retentate and permeate flows, compositions, and pressures.
    """
    # Constants
    n_points = 1000  # Number of discretized steps
    dx = 1 / n_points  # Step size for discretization

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

        permeation = {}
        for comp, y in current_fraction.items():
            # Partial pressure driving force
            p_feed = feed_side_pressure * y  # Partial pressure on feed side (bar)
            p_perm = permeate_side_pressure * permeate_composition.get(comp, 0)  # Partial pressure on permeate side (bar)
            driving_force = max(p_feed - p_perm, 0)  # Ensure non-negative driving force

            # Permeation flux calculation
            perm = permeability.get(comp, 0)  # Permeability (m³(STP)/m²/s/Pa)
            flux = (perm / membrane_thickness) * driving_force * 1e5  # Convert bar to Pa
            permeation[comp] = flux * membrane_area * dx  # Molar permeation rate (kmol/hr)

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

    # Plotting mole fractions
    x = [i * dx for i in range(n_points)]  # Membrane position from 0 to 1
    plt.figure(figsize=(12, 8))
    for comp in feed_composition:
        plt.plot(x, mole_fraction_retentate[comp], label=f'{comp} (Retentate)', linestyle='--')
        plt.plot(x, mole_fraction_permeate[comp], label=f'{comp} (Permeate)', linestyle='-')
    plt.xlabel('Membrane Position')
    plt.ylabel('Mole Fraction')
    plt.title('Mole Fractions Across the Membrane')
    plt.legend()
    plt.grid(True)
    plt.show()

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

feed_pressure = 65  # bar
feed_temperature = 200  # Celsius
membrane_area = 600  # m^2
stage_cut = 0.47

# Calculate the total molar flow rate
total_molar_flow_rate = sum(feed_composition.values())
print(total_molar_flow_rate)
# Calculate the mole fractions
mole_fractions = {component: flow_rate / total_molar_flow_rate for component, flow_rate in feed_composition.items()}

# Get permeability data
permeability = get_permeability_data(feed_temperature, feed_pressure, mole_fractions)

results_DF_PP = simulate_membrane_module_partial_pressure_debug(
    feed_composition, feed_pressure, feed_temperature, permeability, membrane_area, stage_cut, 1e-6
)

# Extract relevant data
feed_flow_rate = results_DF_PP['feed_flow_rate']
retentate_flow_rate = results_DF_PP['retentate_flow_rate']
permeate_flow_rate = results_DF_PP['permeate_flow_rate']
retentate_composition = results_DF_PP['retentate_composition']
permeate_composition = results_DF_PP['permeate_composition']

print('feed_flow_rate')
print(feed_flow_rate)
print('ret_flow rate')
print(retentate_flow_rate)
print('perm_flow rate')
print(permeate_flow_rate)
print('ret composition')
print(retentate_composition)
print('perm comp')
print(permeate_composition)

def perform_rsm_analysis(feed_composition, feed_pressure, feed_temperature, permeability, design_points, function):
    """
    Perform RSM analysis using the provided function and design points.

    Parameters:
        feed_composition (dict): Feed composition (e.g., {'O2': 0.21, 'N2': 0.79}).
        feed_pressure (float): Feed pressure in bar.
        feed_temperature (float): Feed temperature in Kelvin.
        permeability (dict): Permeabilities of gases (e.g., {'O2': 1.0, 'N2': 0.1}).
        design_points (list of tuples): Design points [(membrane_area, stage_cut), ...].
        function (callable): The simulation function.

    Returns:
        pd.DataFrame: DataFrame with design points and %recovery as the response.
    """
    results = []

    for membrane_area, stage_cut in design_points:
        # Run the simulation function
        simulation_output = function(
            feed_composition, feed_pressure, feed_temperature, 
            permeability, membrane_area, stage_cut, 1e-6
        )
        
        # Calculate %recovery
        permeate_flow_rate = simulation_output['permeate_flow_rate']
        feed_flow_rate = simulation_output['feed_flow_rate']
        percent_recovery = (permeate_flow_rate / feed_flow_rate) * 100

        # Append results
        results.append({
            'Membrane Area': membrane_area,
            'Stage Cut': stage_cut,
            '%Recovery': percent_recovery
        })
    
    # Convert results to a DataFrame
    results_df = pd.DataFrame(results)
    return results_df

design_points = [
    (5, 0.1), (50, 0.1), (5, 0.9), (50, 0.9), (27.5, 0.1), 
    (27.5, 0.9), (5, 0.5), (50, 0.5), (27.5, 0.5), (15, 0.5), 
    (40, 0.5), (27.5, 0.3), (27.5, 0.7)
]
'''
results_df = perform_rsm_analysis(
    feed_composition, feed_pressure, feed_temperature, 
    permeability, design_points, simulate_membrane_module_partial_pressure_debug
)

# Display the results
#print(results_df)
'''
# Membrane Gas Separation Simulation

A Python-based simulation of a hollow fiber membrane module for gas separation, with optimization using Response Surface Methodology (RSM).

## Overview

This project models a membrane separation process for recovering hydrogen from a multi-component gas mixture (H2, CO2, CH3OH, H2O, CO). The simulation incorporates:
- Real gas behavior using fugacity calculations
- Hollow fiber membrane geometry
- Pressure drop calculations
- Permeability interpolation
- RSM-based optimization

## Installation

### Prerequisites
- Python 3.8+
- CoolProp (for thermodynamic properties)

### Dependencies
- Python 3.x
- Required libraries:
  - NumPy
  - Pandas
  - Matplotlib
  - SciPy
  - CoolProp

  ## Key Features

| Feature | Description |
|---------|-------------|
| Real Gas Modeling | Uses fugacity instead of partial pressures |
| Hollow Fiber Geometry | Calculates fiber dimensions from membrane area |
| Pressure Drop | Accounts for laminar/turbulent flow in fibers |
| Permeability Interpolation | Uses experimental data at multiple temperatures |
| RSM Optimization | Finds optimal membrane area and stage cut |

## Usage

1. Set the feed composition in the `feed_composition` dictionary
2. Adjust operating conditions in the INPUT SECTION
3. Run the script to:
   - Calculate fiber dimensions
   - Simulate the membrane separation
   - Perform RSM optimization
   - Generate output reports and plots

**Important Note:**
Ensure the retentate molar flow rate (`retentate_flow_rate`) never reaches zero during optimization. A zero value indicates unrealistic physical conditions (complete permeation) and will lead to false optimal values for stage cut and membrane area. Check the `results_df_2` DataFrame for valid retentate flows (> 0) before interpreting RSM results.

## Input Parameters

### Feed Composition (kmol/hr)
feed_composition = {
    'Hydrogen': 137.607,
    'CarbonDioxide': 45.8676,
    'Methanol': 0.747245,
    'Water': 0.147976,
    'CarbonMonoxide': 5.09148,
}

### Operating Conditions
- Feed pressure: 65 bar
- Feed temperature: 200°C
- Membrane area: 103 m²
- Stage cut: 0.9
- Number of fibers: 10,000
- Fiber length: 3 m

## Outputs

The simulation provides:
- Flow rates (feed, retentate, permeate)
- Compositions (mole fractions)
- Pressures
- Optimal membrane area and stage cut
- Hydrogen recovery and loss percentages

## Optimization Approach

The RSM analysis varies:
- Membrane area (5-120 m²)
- Stage cut (0.1-0.9)

The optimization minimizes hydrogen loss while considering:
- Membrane area (capital cost)
- Stage cut (energy requirement)

## Caveats

- **Zero Retentate Flow**  
The optimization assumes non-zero retentate flow. If `retentate_flow_rate = 0`, the model violates mass balance, leading to:
  - Physically meaningless stage cuts (e.g., > 100%)
  - Overestimated membrane areas
  - Numerical instabilities in RSM interpolation

**Solution:**
- Filter RSM results (`results_df_2`) to exclude points where `ret_flowrate ≤ 0`
- Use the provided `find_optimal_point_and_plot()` function, which automatically handles this check
# Ideal Gas Membrane Separation Simulation

Python implementation for membrane gas separation using partial pressure driving force (ideal gas assumption).

## Key Features

- **Ideal Gas Modeling**: Uses partial pressure instead of fugacity
- **Component Permeability**: Interpolates experimental data for H₂, CO₂, CH₃OH, H₂O, CO
- **Visualization**: Plots composition profiles across membrane
- **Response Surface Methodology**: Analyzes performance across parameter space

## Theory

### Driving Force Calculation
ΔPᵢ = P_feed × yᵢ - P_perm × xᵢ [bar]

Where:
- ΔPᵢ = Partial pressure difference for component i
- yᵢ = Mole fraction in feed
- xᵢ = Mole fraction in permeate

### Flux Equation
Jᵢ = (Pᵢ / δ) × ΔPᵢ × 10⁵ [kmol/m²/hr]

Where:
- Pᵢ = Permeability [mol/m²/s/Pa]
- δ = Membrane thickness [m]

## Usage

### Input Parameters
```python
feed_composition = {
    'Hydrogen': 111.971,      # kmol/hr
    'CarbonDioxide': 181.941, # kmol/hr
    'Methanol': 1.04958,      # kmol/hr
    'Water': 0.201152,        # kmol/hr
    'CarbonMonoxide': 15.8056 # kmol/hr
}

# Operating Conditions
feed_pressure = 65      # bar
feed_temperature = 200  # °C
membrane_area = 600     # m²
stage_cut = 0.47        # -
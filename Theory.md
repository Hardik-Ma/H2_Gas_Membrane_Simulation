# Theoretical Foundations

## 1. Permeability Calculations

### Barrer to SI Unit Conversion
$$ P_{SI} \text{ [mol/m²/s/Pa]} = P_{Barrer} \times 3.3464 \times 10^{-16} $$

### Permeability Interpolation
For component $i$ at temperature $T$:
$$ P_i(x) = y_0 + \frac{(y_1 - y_0)}{(x_1 - x_0)} \times (x - x_0) $$

## 2. Fugacity Driving Force

For each component $i$:
$$ f_{feed} = z_{feed} \times p_{feed} \times 10^5 \text{ [Pa]} $$
$$ f_{perm} = z_{perm} \times p_{perm} \times 10^5 \text{ [Pa]} $$
$$ \Delta f_i = \max(f_{feed} - f_{perm}, 0) $$

## 3. Permeation Flux

$$ J_i = \frac{P_i}{\delta} \times \Delta f_i \text{ [mol/m²/s]} $$

## 4. Hollow Fiber Geometry

### Total Membrane Area
$$ A_{total} = N \times \pi \times d \times L $$

### Fiber Diameter/Length
$$ d = \frac{A_{total}}{N \times \pi \times L} \quad \text{[m]} $$
$$ L = \frac{A_{total}}{N \times \pi \times d} \quad \text{[m]} $$

## 5. Pressure Drop Calculations

### Reynolds Number
$$ Re = \frac{\rho \times v \times d}{\mu} $$

### Friction Factor
**Laminar flow** ($Re < 2000$):
$$ f = \frac{64}{Re} $$

**Turbulent flow** ($Re > 4000$):
$$ \frac{1}{\sqrt{f}} = -1.8 \times \log_{10}\left[\left(\frac{6.9}{Re}\right)^{1.11} + \frac{\epsilon/d}{3.7}\right] $$

### Pressure Drop
$$ \Delta P = f \times \frac{L}{d} \times \frac{\rho v^2}{2} \text{ [Pa]} $$

## 6. Wilke's Mixing Rule

$$ \mu_{mix} = \sum_{i=1}^n \frac{x_i \mu_i}{\sum_{j=1}^n x_j \Phi_{ij}} $$

Where:
$$ \Phi_{ij} = \frac{\left[1 + \left(\frac{\mu_i}{\mu_j}\right)^{1/2} \left(\frac{M_j}{M_i}\right)^{1/4}\right]^2}{\sqrt{8}\left(1 + \frac{M_i}{M_j}\right)^{1/2}} $$

## 7. RSM Optimization

Objective function:
$$ \text{minimize: } \%H_2 loss + w_1 A + w_2 \theta $$

Constraints:
$$ Q_{retentate} > 0 $$
$$ 0.1 \leq \theta \leq 0.9 $$

## 8. Stage Cut
$$ \theta = \frac{Q_{perm}}{Q_{feed}} $$

## 9. Component Recovery
$$ \%Recovery_i = \frac{Q_{perm} \times y_i}{Q_{feed} \times x_i} \times 100 $$

## Nomenclature

| Symbol | Description | Units |
|--------|-------------|-------|
| $P$ | Permeability | mol/m²/s/Pa |
| $f$ | Fugacity | Pa |
| $J$ | Flux | mol/m²/s |
| $d$ | Fiber diameter | m |
| $Re$ | Reynolds number | - |
| $\mu$ | Viscosity | Pa·s |
| $\Delta P$ | Pressure drop | Pa |
| $\theta$ | Stage cut | - |
# Monte Carlo and Sensitivity Analysis for Geothermal Power Potential

This repository provides MATLAB scripts for rigorous estimation and analysis of geothermal power plant capacity using Monte Carlo simulations and global sensitivity analysis. The models are adaptable for various geothermal reservoirs and have been validated against the La Palma field and applied to the Al-Lisi reservoir in Yemen.

## Contents

- **geothermal_montecarlo_alisi_yemen.m**:  
  Main script for Monte Carlo simulation and convergence testing for the Al-Lisi geothermal reservoir (Yemen).  
  - Calculates distributions for power capacity, thermal energy, and energy breakdown.
  - Includes convergence diagnostics, statistics, and publication-ready visualizations.
  - All results in MWe.

- **geothermal_montecarlo_plot_pub.m**:  
  Publication-ready plotting utility for Monte Carlo simulation results (MWe and energy).
  - Generates PDF/CDF/thermal plots with labeled confidence intervals, modes, and statistics.
  - All annotations are formatted for figure inclusion in scientific publications.

- **geothermal_montecarlo_la_palma_validated.m**:  
  Monte Carlo simulation script validated against the La Palma geothermal field.  
  - Uses parameter distributions from La Palma studies.
  - Reports statistics and visualizes harvested power and energy distributions.

- **global_sensitivity_morris_alisi_yemen.m**:  
  Morris (elementary effects) global sensitivity analysis for the Al-Lisi model.
  - Focuses on the top 5 most influential parameters.
  - Includes ranking, μ* vs σ plot, and publication-quality tables/figures.

## Usage

1. **Monte Carlo Simulation**
   - Run `geothermal_montecarlo_main_ALisi_YEMEN.m` for a full simulation with convergence testing.
   - The script prints all key results and saves arrays for further plotting or analysis.

2. **Plotting**
   - Use `geothermal_montecarlo_allblueplots_mwe.m` to create high-quality plots from simulation outputs.

3. **Sensitivity Analysis**
   - Run `Global_sensitivity_Model_Alisi_Yemen.m` to assess parameter importance and produce sensitivity plots.

4. **Validation**
   - Use `MonteCarlo_La_Palma_VALIDATED.m` to reproduce results for the La Palma field or validate the model against other published studies.

## Author

All scripts and models were developed by **Mugahed Kolaibi** (moKolaibi).  
Please cite the author if using these codes in scientific publications.

## File Naming Convention

- Scripts follow the format:  
  - `geothermal_montecarlo_<case>.m` for Monte Carlo simulation.
  - `geothermal_montecarlo_plot_pub.m` for plotting utilities.
  - `global_sensitivity_morris_<case>.m` for sensitivity analysis.
  - The La Palma validated model retains the `geothermal_montecarlo_la_palma_validated.m` name for clarity.

## License

[MIT License](LICENSE)

---

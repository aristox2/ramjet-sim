# Ramjet Engine Simulation 🚀

A **comprehensive Python module** for simulating and analyzing the performance of a **ramjet engine**. This module models the full thermodynamic cycle—including inlet, combustion chamber, and nozzle—allowing for thrust, efficiency, and performance calculations.

---

![Ramjet Diagram](https://upload.wikimedia.org/wikipedia/commons/thumb/6/6c/Ramjet_operation.svg/1920px-Ramjet_operation.svg.png)
*Simplified Ramjet Operation Diagram (Wikipedia)*

---

## Features

- **Inlet Analysis**
  - Handles subsonic and supersonic flow
  - Normal shock relations for supersonic Mach numbers
  - Calculates pressure and temperature ratios

- **Combustion Chamber**
  - Fuel-air mixing and stoichiometry
  - Heat addition and combustion efficiency
  - Temperature rise and pressure losses

- **Nozzle Expansion**
  - Isentropic expansion
  - Choked flow detection
  - Exit velocity and temperature calculations

- **Performance Metrics**
  - Thrust and specific impulse
  - Thrust-specific fuel consumption (TSFC)
  - Thermal, propulsive, and overall efficiencies

- **User-Friendly Interface**
  - Input parameters via Python dataclass
  - Descriptive output results
  - Sample analysis included

---

## Installation

Requires Python 3.8+ and the following packages:

```bash
pip install numpy matplotlib
git clone (https://github.com/aristox2/ramjet-sim/)
cd ramjet-engine-sim

```

# Run a sample ramjet analysis

from ramjet_engine import create_sample_analysis
engine, results = create_sample_analysis()



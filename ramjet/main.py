import numpy as np
from typing import Dict, Tuple, Optional
from dataclasses import dataclass
import matplotlib.pyplot as plt

# Physical constants
GAMMA_AIR = 1.4  # Specific heat ratio for air
GAMMA_COMBUSTION = 1.3  # Specific heat ratio for combustion products
R_AIR = 287.0  # Gas constant for air (J/kg·K)
CP_AIR = 1005.0  # Specific heat at constant pressure for air (J/kg·K)
CP_COMBUSTION = 1156.0  # Specific heat at constant pressure for combustion products (J/kg·K)

@dataclass
class RamjetParameters:
    """Input parameters for ramjet engine analysis."""
    # Flight conditions
    altitude: float  # meters
    mach_number: float  # flight Mach number
    ambient_temperature: float  # Kelvin
    ambient_pressure: float  # Pa
    
    # Engine geometry
    inlet_area_ratio: float  # A1/A0 (inlet to free stream area ratio)
    combustion_chamber_area_ratio: float  # A3/A1
    nozzle_exit_area_ratio: float  # A4/A3
    
    # Fuel properties
    fuel_heating_value: float  # J/kg
    fuel_air_ratio: float  # stoichiometric fuel-air ratio
    actual_fuel_air_ratio: float  # actual fuel-air ratio
    
    # Component efficiencies
    inlet_efficiency: float = 0.95
    combustion_efficiency: float = 0.98
    nozzle_efficiency: float = 0.98

class RamjetEngine:
    """
    Ramjet engine performance analysis.
    
    This class implements the complete thermodynamic cycle analysis
    for a ramjet engine including inlet compression, combustion,
    and nozzle expansion.
    """
    
    def __init__(self, params: RamjetParameters):
        """
        Initialize the ramjet engine model.
        
        Args:
            params: RamjetParameters object containing all input parameters
        """
        self.params = params
        self.results = {}
        
    def calculate_ambient_conditions(self) -> Dict[str, float]:
        """
        Calculate ambient atmospheric conditions based on altitude.
        
        Returns:
            Dictionary containing ambient temperature and pressure
        """
        # Standard atmosphere model (simplified)
        if self.params.altitude <= 11000:  # Troposphere
            ambient_temperature = 288.15 - 0.0065 * self.params.altitude
            ambient_pressure = 101325 * (ambient_temperature / 288.15) ** 5.256
        else:  # Stratosphere
            ambient_temperature = 216.65
            ambient_pressure = 22632 * np.exp(-0.0001577 * (self.params.altitude - 11000))
            
        return {
            'temperature': ambient_temperature,
            'pressure': ambient_pressure,
            'density': ambient_pressure / (R_AIR * ambient_temperature)
        }
    
    def inlet_analysis(self) -> Dict[str, float]:
        """
        Analyze inlet compression and shock wave formation.
        
        Returns:
            Dictionary containing inlet performance parameters
        """
        flight_mach = self.params.mach_number
        ambient_conditions = self.calculate_ambient_conditions()
        
        # Normal shock relations
        if flight_mach > 1:
            # Normal shock relations
            inlet_exit_mach = np.sqrt((flight_mach**2 + 2/(GAMMA_AIR-1)) / (2*GAMMA_AIR/(GAMMA_AIR-1) * flight_mach**2 - 1))
            inlet_pressure_ratio = 1 + 2*GAMMA_AIR/(GAMMA_AIR+1) * (flight_mach**2 - 1)
            inlet_temperature_ratio = (1 + 2*GAMMA_AIR/(GAMMA_AIR+1) * (flight_mach**2 - 1)) * (2 + (GAMMA_AIR-1)*flight_mach**2) / ((GAMMA_AIR+1)*flight_mach**2)
        else:
            # Isentropic compression for subsonic flow
            inlet_exit_mach = flight_mach
            inlet_pressure_ratio = (1 + (GAMMA_AIR-1)/2 * flight_mach**2)**(GAMMA_AIR/(GAMMA_AIR-1))
            inlet_temperature_ratio = (1 + (GAMMA_AIR-1)/2 * flight_mach**2)
        
        # Apply inlet efficiency
        inlet_pressure_ratio *= self.params.inlet_efficiency
        inlet_temperature_ratio *= self.params.inlet_efficiency
        
        return {
            'mach_1': inlet_exit_mach,
            'pressure_ratio': inlet_pressure_ratio,
            'temperature_ratio': inlet_temperature_ratio,
            'pressure_1': ambient_conditions['pressure'] * inlet_pressure_ratio,
            'temperature_1': ambient_conditions['temperature'] * inlet_temperature_ratio
        }
    
    def combustion_analysis(self, inlet_results: Dict[str, float]) -> Dict[str, float]:
        """
        Analyze combustion chamber thermodynamics.
        
        Args:
            inlet_results: Results from inlet analysis
            
        Returns:
            Dictionary containing combustion performance parameters
        """
        # Combustion chamber conditions
        combustor_inlet_temperature = inlet_results['temperature_1']
        combustor_inlet_pressure = inlet_results['pressure_1']
        
        # Fuel-air mixture properties
        actual_fuel_air_ratio = self.params.actual_fuel_air_ratio
        stoichiometric_fuel_air_ratio = self.params.fuel_air_ratio
        
        # Heat addition in combustion chamber
        heat_addition_rate = actual_fuel_air_ratio * self.params.fuel_heating_value * self.params.combustion_efficiency
        
        # Temperature rise in combustion chamber
        # Using energy balance: m_dot * cp * (T3 - T2) = m_dot_fuel * HV * eta_comb
        combustor_exit_temperature = combustor_inlet_temperature + heat_addition_rate / CP_COMBUSTION
        
        # Pressure loss in combustion chamber (simplified)
        combustor_exit_pressure = combustor_inlet_pressure * 0.95  # 5% pressure loss
        
        return {
            'temperature_3': combustor_exit_temperature,
            'pressure_3': combustor_exit_pressure,
            'heat_addition': heat_addition_rate,
            'combustion_temperature_rise': combustor_exit_temperature - combustor_inlet_temperature
        }
    
    def nozzle_analysis(self, combustion_results: Dict[str, float]) -> Dict[str, float]:
        """
        Analyze nozzle expansion and thrust calculation.
        
        Args:
            combustion_results: Results from combustion analysis
            
        Returns:
            Dictionary containing nozzle performance parameters
        """
        combustor_exit_temperature = combustion_results['temperature_3']
        combustor_exit_pressure = combustion_results['pressure_3']
        ambient_conditions = self.calculate_ambient_conditions()
        
        # Isentropic expansion in nozzle
        nozzle_pressure_ratio = ambient_conditions['pressure'] / combustor_exit_pressure
        
        if nozzle_pressure_ratio < (2/(GAMMA_COMBUSTION+1))**(GAMMA_COMBUSTION/(GAMMA_COMBUSTION-1)):
            # Choked flow
            nozzle_exit_mach = 1.0
            nozzle_exit_pressure = combustor_exit_pressure * (2/(GAMMA_COMBUSTION+1))**(GAMMA_COMBUSTION/(GAMMA_COMBUSTION-1))
        else:
            # Subsonic expansion
            nozzle_exit_mach = np.sqrt(2/(GAMMA_COMBUSTION-1) * ((combustor_exit_pressure/ambient_conditions['pressure'])**((GAMMA_COMBUSTION-1)/GAMMA_COMBUSTION) - 1))
            nozzle_exit_pressure = ambient_conditions['pressure']
        
        # Temperature at nozzle exit
        nozzle_exit_temperature = combustor_exit_temperature * (nozzle_exit_pressure/combustor_exit_pressure)**((GAMMA_COMBUSTION-1)/GAMMA_COMBUSTION)
        
        # Exit velocity
        nozzle_exit_velocity = nozzle_exit_mach * np.sqrt(GAMMA_COMBUSTION * R_AIR * nozzle_exit_temperature)
        
        return {
            'mach_4': nozzle_exit_mach,
            'temperature_4': nozzle_exit_temperature,
            'pressure_4': nozzle_exit_pressure,
            'exit_velocity': nozzle_exit_velocity,
            'nozzle_efficiency': self.params.nozzle_efficiency
        }
    
    def calculate_thrust(self, inlet_results: Dict[str, float], 
                        nozzle_results: Dict[str, float]) -> Dict[str, float]:
        """
        Calculate thrust and specific impulse.
        
        Args:
            inlet_results: Results from inlet analysis
            nozzle_results: Results from nozzle analysis
            
        Returns:
            Dictionary containing thrust performance parameters
        """
        ambient_conditions = self.calculate_ambient_conditions()
        
        # Mass flow rate (simplified - assuming constant area)
        reference_area = 1.0  # Reference area
        inlet_density = inlet_results['pressure_1'] / (R_AIR * inlet_results['temperature_1'])
        inlet_velocity = inlet_results['mach_1'] * np.sqrt(GAMMA_AIR * R_AIR * inlet_results['temperature_1'])
        air_mass_flow_rate = inlet_density * inlet_velocity * reference_area
        
        # Thrust calculation
        # F = m_dot * (V4 - V0) + (P4 - P0) * A4
        flight_velocity = self.params.mach_number * np.sqrt(GAMMA_AIR * R_AIR * ambient_conditions['temperature'])
        nozzle_exit_area = reference_area * self.params.inlet_area_ratio * self.params.combustion_chamber_area_ratio * self.params.nozzle_exit_area_ratio
        
        momentum_thrust_component = air_mass_flow_rate * (nozzle_results['exit_velocity'] - flight_velocity)
        pressure_thrust_component = (nozzle_results['pressure_4'] - ambient_conditions['pressure']) * nozzle_exit_area
        total_thrust = momentum_thrust_component + pressure_thrust_component
        
        # Specific impulse
        fuel_mass_flow_rate = air_mass_flow_rate * self.params.actual_fuel_air_ratio
        specific_impulse = total_thrust / (fuel_mass_flow_rate * 9.81)  # seconds
        
        # Thrust specific fuel consumption
        thrust_specific_fuel_consumption = fuel_mass_flow_rate / total_thrust  # kg/N·s
        
        return {
            'thrust': total_thrust,
            'specific_impulse': specific_impulse,
            'tsfc': thrust_specific_fuel_consumption,
            'mass_flow_rate': air_mass_flow_rate,
            'fuel_flow_rate': fuel_mass_flow_rate,
            'exit_velocity': nozzle_results['exit_velocity']
        }
    
    def calculate_efficiency(self, inlet_results: Dict[str, float],
                           combustion_results: Dict[str, float],
                           thrust_results: Dict[str, float]) -> Dict[str, float]:
        """
        Calculate various efficiency parameters.
        
        Args:
            inlet_results: Results from inlet analysis
            combustion_results: Results from combustion analysis
            thrust_results: Results from thrust calculation
            
        Returns:
            Dictionary containing efficiency parameters
        """
        ambient_conditions = self.calculate_ambient_conditions()
        
        # Thermal efficiency
        combustor_exit_temperature = combustion_results['temperature_3']
        combustor_inlet_temperature = inlet_results['temperature_1']
        ambient_temperature = ambient_conditions['temperature']
        
        thermal_efficiency = 1 - ambient_temperature/combustor_exit_temperature
        
        # Propulsive efficiency
        flight_velocity = self.params.mach_number * np.sqrt(GAMMA_AIR * R_AIR * ambient_conditions['temperature'])
        exhaust_velocity = thrust_results['exit_velocity']
        
        propulsive_efficiency = 2 * flight_velocity / (flight_velocity + exhaust_velocity)
        
        # Overall efficiency
        overall_efficiency = thermal_efficiency * propulsive_efficiency
        
        return {
            'thermal_efficiency': thermal_efficiency,
            'propulsive_efficiency': propulsive_efficiency,
            'overall_efficiency': overall_efficiency
            
        }
    
    def run_analysis(self) -> Dict[str, Dict[str, float]]:
        """
        Run complete ramjet engine analysis.
        
        Returns:
            Dictionary containing all analysis results
        """
        # Step 1: Inlet analysis
        inlet_results = self.inlet_analysis()
        
        # Step 2: Combustion analysis
        combustion_results = self.combustion_analysis(inlet_results)
        
        # Step 3: Nozzle analysis
        nozzle_results = self.nozzle_analysis(combustion_results)
        
        # Step 4: Thrust calculation
        thrust_results = self.calculate_thrust(inlet_results, nozzle_results)
        
        # Step 5: Efficiency calculation
        efficiency_results = self.calculate_efficiency(inlet_results, combustion_results, thrust_results)
        
        # Compile all results
        self.results = {
            'inlet': inlet_results,
            'combustion': combustion_results,
            'nozzle': nozzle_results,
            'thrust': thrust_results,
            'efficiency': efficiency_results,
            'ambient': self.calculate_ambient_conditions()
        }
        
        return self.results
    
    def print_results(self):
        """Print formatted analysis results."""
        if not self.results:
            print("No analysis results available. Run run_analysis() first.")
            return
        
        print("=" * 60)
        print("RAMJET ENGINE PERFORMANCE ANALYSIS")
        print("=" * 60)
        
        print(f"\nFlight Conditions:")
        print(f"  Altitude: {self.params.altitude:.0f} m")
        print(f"  Mach Number: {self.params.mach_number:.2f}")
        print(f"  Ambient Temperature: {self.results['ambient']['temperature']:.1f} K")
        print(f"  Ambient Pressure: {self.results['ambient']['pressure']:.0f} Pa")
        
        print(f"\nInlet Performance:")
        print(f"  Inlet Mach Number: {self.results['inlet']['mach_1']:.3f}")
        print(f"  Pressure Ratio: {self.results['inlet']['pressure_ratio']:.3f}")
        print(f"  Temperature Ratio: {self.results['inlet']['temperature_ratio']:.3f}")
        
        print(f"\nCombustion Performance:")
        print(f"  Combustion Temperature: {self.results['combustion']['temperature_3']:.1f} K")
        print(f"  Combustion Pressure: {self.results['combustion']['pressure_3']:.0f} Pa")
        print(f"  Heat Addition: {self.results['combustion']['heat_addition']:.0f} J/kg")
        
        print(f"\nNozzle Performance:")
        print(f"  Exit Mach Number: {self.results['nozzle']['mach_4']:.3f}")
        print(f"  Exit Velocity: {self.results['nozzle']['exit_velocity']:.0f} m/s")
        print(f"  Exit Temperature: {self.results['nozzle']['temperature_4']:.1f} K")
        
        print(f"\nThrust Performance:")
        print(f"  Thrust: {self.results['thrust']['thrust']:.0f} N")
        print(f"  Specific Impulse: {self.results['thrust']['specific_impulse']:.0f} s")
        print(f"  TSFC: {self.results['thrust']['tsfc']*1000:.2f} mg/N·s")
        print(f"  Mass Flow Rate: {self.results['thrust']['mass_flow_rate']:.2f} kg/s")
        
        print(f"\nEfficiency:")
        print(f"  Thermal Efficiency: {self.results['efficiency']['thermal_efficiency']:.3f}")
        print(f"  Propulsive Efficiency: {self.results['efficiency']['propulsive_efficiency']:.3f}")
        print(f"  Overall Efficiency: {self.results['efficiency']['overall_efficiency']:.3f}")
        print("=" * 60)

def create_sample_analysis():
    """Create and run a sample ramjet analysis."""
    
    # Define sample parameters
    params = RamjetParameters(
        altitude=10000,  # 10 km altitude
        mach_number=2.5,  # Mach 2.5 flight
        ambient_temperature=223.15,  # K
        ambient_pressure=26436,  # Pa
        
        inlet_area_ratio=1.2,
        combustion_chamber_area_ratio=1.5,
        nozzle_exit_area_ratio=2.0,
        
        fuel_heating_value=43e6,  # J/kg (typical for hydrocarbon fuel)
        fuel_air_ratio=0.067,  # stoichiometric F/A ratio
        actual_fuel_air_ratio=0.04,  # actual F/A ratio
        
        inlet_efficiency=0.95,
        combustion_efficiency=0.98,
        nozzle_efficiency=0.98
    )
    
    # Create engine model and run analysis
    engine = RamjetEngine(params)
    results = engine.run_analysis()
    
    # Print results
    engine.print_results()
    
    return engine, results

if __name__ == "__main__":
    # Run sample analysis

    engine, results = create_sample_analysis()

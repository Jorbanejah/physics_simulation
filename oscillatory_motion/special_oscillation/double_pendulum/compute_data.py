import numpy as np
import os
from typing import Sequence, Dict, Any
from double_pendulum import DoublePendulumSimulator

def compute_initial_angles(params, theta1: Sequence[float], theta2: Sequence[float], linearized: bool = False, directory: str = os.getcwd(), filename: str = "compute_initial_angles.npz") -> Dict[str, Any]:
    """
    Compute energy drift for varying initial conditions.

    Parameters
    ----------
    params : Params
        Physical parameters of the pendulum.
    theta1 : Sequence[float]
        List of initial angles for rod 1 (radians).
    theta2 : Sequence[float]
        List of initial angles for rod 2 (radians).
    linearized : bool
        If True, use the small-angle approximation equations.
    directory : str
        Path to save the output file.
    filename : str
        Name of the output file.

    Returns
    -------
    Dict
        Nested dictionary containing scan results and grid results.
    """
    if len(theta1) == 0 or len(theta2) == 0:
        raise ValueError("theta1 and theta2 must each contain at least one angle")

    # Initialize storage structures
    theta1_scan = {th: {} for th in theta1}
    theta2_scan = {th: {} for th in theta2}
    grid = {th1: {th2: {} for th2 in theta2} for th1 in theta1}

    def run_simulation(q0: tuple) -> dict:
        """Helper to run simulation and extract data."""
        params.theta1_0, params.theta2_0 = q0
        sim = DoublePendulumSimulator(params=params)
        # Pass linearized flag if supported by your DoublePendulumSimulator
        result = sim.run(linearized=linearized) 
        
        try:
            theta1_sol = result.y[0] if hasattr(result, 'y') else result.theta1
            theta2_sol = result.y[1] if hasattr(result, 'y') else result.theta2
            omega1_sol = result.y[2] if hasattr(result, 'y') else result.omega1
            omega2_sol = result.y[3] if hasattr(result, 'y') else result.omega2
        except:
            theta1_sol = result.theta1
            theta2_sol = result.theta2
            omega1_sol = result.omega1
            omega2_sol = result.omega2

        # Energy calculation
        total_energy = result.total # Access total energy array
        energy_drift = np.max(total_energy) - np.min(total_energy)
        
        return {
            "theta1": theta1_sol, "theta2": theta2_sol, 
            "omega1": omega1_sol, "omega2": omega2_sol, 
            "drift": energy_drift
        }

    # --- Process Scans ---
    
    # 1. Scan Theta1 (Theta2 = 0)
    print("=" * 60)
    mode = "Linearized" if linearized else "Non-linear"
    print(f"Scanning Theta 1 ({mode} model)")
    print("=" * 60)
    
    for i, th in enumerate(theta1):
        print(f"Progress: {100*(i+1)/len(theta1):.1f}%", end="\r")
        theta1_scan[th] = run_simulation((th, 0.0))

    # 2. Scan Theta2 (Theta1 = 0)
    print("\n" + "=" * 60)
    print(f"Scanning Theta 2 ({mode} model)")
    print("=" * 60)
    
    for i, th in enumerate(theta2):
        # Fixed bug: was using len(theta1) instead of len(theta2)
        print(f"Progress: {100*(i+1)/len(theta2):.1f}%", end="\r") 
        theta2_scan[th] = run_simulation((0.0, th))

    # 3. Full Grid (Theta1 x Theta2)
    print("\n" + "=" * 60)
    print(f"Building Full Grid ({mode} model)")
    print("=" * 60)
    
    total_steps = len(theta1) * len(theta2)
    count = 0
    
    for th1 in theta1:
        for th2 in theta2:
            count += 1
            print(f"Grid Progress: {100*count/total_steps:.1f}%", end="\r")
            grid[th1][th2] = run_simulation((th1, th2))

    print("\nSaving data...")
    route = os.path.join(directory, filename)
    # Save as numpy compressed archive
    np.savez(route, theta1_scan=theta1_scan, theta2_scan=theta2_scan, grid=grid)
    
    return {"theta1_scan": theta1_scan, "theta2_scan": theta2_scan, "grid": grid}
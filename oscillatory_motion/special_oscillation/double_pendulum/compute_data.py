import numpy as np
import os
from typing import Sequence, Dict, Any
from double_pendulum import DoublePendulumSimulator

"""
This code is only use in main_error.py file
"""
def compute_initial_angles(params, theta1: Sequence[float], theta2: Sequence[float], directory: str = os.getcwd(), filename: str = "compute_intial_angles.npz") -> Dict[str, Any]:
    """
    Compute energy drift for:
    - varying theta1 with theta2 = 0
    - varying theta2 with theta1 = 0
    - full grid (theta1, theta2)

    Returns a dictionary with three entries:
        "theta1_scan" - "theta1": [], "theta2": [], "omega1": [], "omega_2": [], "kinetic": [], "potencial": [], "Drift": []
        "theta2_scan" - "theta1": [], "theta2": [], "omega1": [], "omega_2": [], "kinetic": [], "potencial": [], "Drift": []
        "grid" - "theta1": [], "theta2": [], "omega1": [], "omega_2": [], "kinetic": [], "potencial": [], "Drift": []
    Uses:
    - Fractal
    - Errors calculus (main_error.py)

    """
    #Validation
    if len(theta1) == 0 or len(theta2) == 0:
        raise ValueError("theta1 and theta2 must each contain at least one angle")

    #Define dictionaries
    theta1_scan = {th: {"theta1": [], "theta2": [], "omega1": [], "omega_2": [], "kinetic": [], "potencial": [], "Drift": []} for th in theta1}
    theta2_scan = {th: {"theta1": [], "theta2": [], "omega1": [], "omega_2": [], "kinetic": [], "potencial": [], "Drift": []} for th in theta2}
    grid = {th1: {th2: {"theta1": [], "theta2": [], "omega1": [], "omega_2": [], "kinetic": [], "potencial": [], "Drift": []} for th2 in theta2} for th1 in theta1}

    def _compute(q0: Sequence[float])->dict:
        params.theta1_0, params.theta2_0 = q0
        sim = DoublePendulumSimulator(params=params)
        results = sim.run()
        T, U, Et = results.kinetic, results.potential, results.total
        theta1, theta2, omega1, omega2 = results.y

        drift = _compute_energy_drift(Et)
        return {"theta1":theta1, "theta2": theta2, "omega1": omega1, "omega2": omega2, "kinetic": T, "potencial":U, "drift": drift}

    def _compute_energy_drift(Et: Sequence[float])->float:

        energy_drift = max(Et) - min(Et)

        return energy_drift

    total = len(theta1)
    bar_len = 20

    # Scan theta1 (theta2 = 0)
    print("=" * 60)
    print("Starting with theta1_0")
    print("=" * 60)
    for i, th in enumerate(theta1):
        progress = (i+1)/len(theta1)
        filled = int(progress * 20)
        bar = "█" * filled + "-" * (bar_len - filled)
        print(f"[{bar}]  {progress*100:5.1f}%   θ₁ = {th:.4f}", end="\r", flush=True)
        theta1_scan[th] = _compute((th, 0.0))

    #Scan theta2 (theta1 = 0)
    print("=" * 60)
    print("Starting with theta2_0")
    print("=" * 60)
    for i, th in enumerate(theta2):
        progress = (i+1)/len(theta1)
        filled = int(progress * 20)
        bar = "█" * filled + "-" * (bar_len - filled)
        print(f"[{bar}]  {progress*100:5.1f}%   θ₁ = {th:.4f}", end="\r", flush=True)
        theta2_scan[th] = _compute((0.0, th))

    # Full grid (theta1, theta2)
    print("=" * 60)
    print("Starting full grid")
    print("=" * 60)
    for i, th1 in enumerate(theta1):
        progress = (i + 1) / total
        filled = int(progress * bar_len)
        bar = "█" * filled + "-" * (bar_len - filled)
        print(f"[{bar}]  {progress*100:5.1f}%   θ₁ = {th1:.4f}", end="\r", flush=True)

        for th2 in theta2:
            grid[th1][th2] = _compute((th1, th2))

    print()  # newline after progress bar

    # -----------------------------
    # 6. Save results
    # -----------------------------
    route = os.path.join(directory, filename)
    np.savez(route, theta1_scan=theta1_scan, theta2_scan=theta2_scan, grid=grid)

    return {"theta1_scan": theta1_scan, "theta2_scan": theta2_scan, "grid": grid}
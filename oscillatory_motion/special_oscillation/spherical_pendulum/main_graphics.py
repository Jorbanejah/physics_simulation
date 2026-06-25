"""
Spherical pendulum - main graphics:

Main graphics will be run by both: linearized and normal simulation.

- Time vs omega, Time vs phi.
- Phase space (4x4) mix (theta, phi, dtheta, dphi)
- Varying initial condition:
    - Phase portrait: nutation, precession
    - Resonance ratio vs dtheta/dphi

Numerical graphics:
- Comparison method.
- Energy drift colormap.
- Runtime.
- Convergence/stability.

Special grahics:
- Invariant-torus reconstruction. Use delay embedding or Fourier decomposition to visualize the torus in 3D.
- Frequency-map analysis (Laskar). Plot frequency drift to detect weak chaos. (MIRA EN FAVORITOS DE GOOGLE)
- Action-angle coordinate plots. If you compute approximate actions, J_theta,J_phi, plot trajectories in action space.
"""

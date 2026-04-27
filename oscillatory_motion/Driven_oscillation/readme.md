Hello, and welcome to this new class of oscillation. Today: DRIVEN OSCILLATION.

This project explores driven oscillations, the final topic in the family of single-degree-of-freedom oscillators that I studied during my physics degree.

The repository is organized into the following sections:

###  Linear
This folder covers the linear driven oscillator, such as the driven spring–mass system or the small-angle approximation of the pendulum.  
Here we discuss the well-known phenomenon of RESONANCE and how it emerges in linear systems. Moreover, we are going to deep into new numerical method call: Verlet integration.

###  Nonlinear
This folder focuses on the nonlinear driven pendulum, where the restoring force is no longer approximated - That is scary, isn't it?  
We explore how the system behaves when the driving force interacts with the full nonlinear dynamics. Although “nonlinear resonance” is not a standard term, we examine what happens when we attempt to reproduce the same analysis used in the linear case. Also, we see how well the numerical methods works through energies.
However, the driven nonlinear pendulum is also the first system in this series that exhibits chaotic motion, so we study:

- Bifurcation diagrams  
- Poincaré sections  
- The Lyapunov exponent  

(These graphics are extremely gorgeous, so I have considerd make a proper animation)
### Figures
Contains the figures generated throughout the project.

### Media + Animation
Stores all video animations produced during the simulations and the code which produce it.

Additionally, the repository includes the data files (.gitignore) used for the numerical computations.

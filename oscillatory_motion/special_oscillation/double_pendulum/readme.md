# Double pendulum: A fateway to chaos

The double pendulum is a classic example of nonlinear chaotic system. Despite its simple structure - two masses connected by rigid rods - it exhibbits highly sensitive depndence on initial conditions, making long-term prediction nearly impossible without precise computation/

- Lagrangian formulation:

Our friend Lagrange make the way a little bit easy. Identifying the system's degrees of freedom (our two horinzontal-angles $\theta_1$ and $\theta_2$), we derive the equations of motion using the Euler-Lagrange equations.

The results - they could be a little odd or horrous -is a set of coupled nonlinar second-order differential equations that cannot be solved analytically in general.

For the first pendulum:

For the second pendulum:

These marvellous equation describe the full nonlinear dynamics of the system, and are typically solved using numerical methods.
We also must say that exist a way to reduce these equations. Using the small-angle approximation, same masses, and same length. We can transform the nonlinear equation a couple of:


# Numerical methods
In this chapter, I decide not to use others numerical methods (Remember that in this serie we see RK4, Verlet, Crank-Nicolson, Euler...), instead of doing that, we are going to see the numerical pythonic methods are.
The methods are (from the scipy.integrate solve_ivp library): RK45, RADAU (implicit), DOP853, and BDF.
With a great touch we will see how this methods performs:

- RK45
- RADAU
- DOP853
- BDF
- How to use it?
# Errors and analysis
Now we have introduced the numerical methods, we can talk about: which one is better than other one, or which one require the least runtime... We discover that and much more than this through the graphics.

- Runtime
- Drift energy
    - Through different dt
    - Through different method

- Convergence and stability

If you want to discover other issues that have not been discussed here, you can check the code and change a few lines.

# The fractal revolution

# Bibliography
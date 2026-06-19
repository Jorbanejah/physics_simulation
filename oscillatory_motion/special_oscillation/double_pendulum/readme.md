## Double pendulum: A gateway to chaos

The double pendulum is a classic example of nonlinear, chaotic system. Despite its simple structure - two masses connected by rigid rods - it exhibits highly sensitive dependence on initial conditions, making long-term prediction nearly impossible without precise computation.

- **Lagrangian formulation:**
  
 As always, identifying the system's degrees of freedom (the angular coordinates $\theta_1$ and $\theta_2$), we derive the equations of motion using the Euler-Lagrange equations.

The results is a set of coupled, nonlinear, second-order differential equations that cannot be solved analytically in general.

- **For the first pendulum:**

$$
\ddot{\theta}_1 =
\frac{-g(2m_1 + m_2)\sin\theta_1 - m_2 g \sin(\theta_1 - 2\theta_2) - 2\sin(\theta_1 - \theta_2)\,m_2\left(\dot{\theta}_2^{\,2} l_2 + \dot{\theta}_1^{\,2} l_1 \cos(\theta_1 - \theta_2)\right)}{ l_1\left(2m_1 + m_2 - m_2\cos(2\theta_1 - 2\theta_2)\right)}
$$

- **For the second pendulum:**

$$
\ddot{\theta}_2 =
\frac{2\sin(\theta_1 - \theta_2)\left(\dot{\theta}_1^{\,2} l_1 (m_1 + m_2)+ g (m_1 + m_2)\cos\theta_1 + \dot{\theta}_2^{\,2} l_2 m_2 \cos(\theta_1 - \theta_2) \right)}{l_2\left(2m_1 + m_2 - m_2\cos(2\theta_1 - 2\theta_2)\right)}
$$

These equations describe the full nonlinear dynamics of the system, and are typically solved using numerical methods.
There also exists a way to simplify these equations. Using the small-angle approximation, same masses, and same length. We can transform the nonlinear equation to a couple of:

- **Small-angle approximation,  $m_1 = m_2 = m$; $l_1 = l_2 = l$**

$$
\ddot{\theta}_1 \approx -\frac{g}{l}\left(3\theta_1 - 2\theta_2\right)
$$

$$
\ddot{\theta}_2 \approx \frac{2g}{l}\left(\theta_1 - \theta_2\right)
$$

## How do we write it?
Firstly, we have to calculate the kinetic energy to get the mass matrix

$$
T = \tfrac12\ \dot{\boldsymbol{\theta}}^{\mathsf{T}} \mathbf{M}(\theta_1,\theta_2)\ \dot{\boldsymbol{\theta}}
$$

For the double pendulum:

$$
T = \tfrac12 (m_1+m_2) L_1^2 \dot{\theta}_1^2 + \tfrac12 m_2 L_2^2 \dot{\theta}_2^2 + m_2 L_1 L_2 \dot{\theta}_1 \dot{\theta}_2 \cos(\theta_1 - \theta_2)
$$

 This yields the mass matrix:
 
$$
\mathbf{M}(\theta_1,\theta_2) = \begin{pmatrix} (m_1+m_2)L_1^2 & m_2 L_1 L_2 \cos(\theta_1 - \theta_2) \\
m_2 L_1 L_2 \cos(\theta_1 - \theta_2) & m_2 L_2^2
\end{pmatrix}
$$

Secondly, we calculate, from the Lagrangian, the force vector - where each component is the generalized force conjugate to the generalized coordinate $\theta_i$:

$$
\mathbf{Q} = -\frac{\partial V}{\partial \boldsymbol{\theta}}                        + \text{Coriolis/Centrifugal terms from } \frac{d}{dt} \left( \frac{\partial T}{\partial \dot{\boldsymbol{\theta}}} \right)
$$

Potential energy:

$$
V = -(m_1+m_2) g L_1 \cos\theta_1 - m_2 g L_2 \cos\theta_2
$$

The vector Q is the generalized force vector. It contains all gravitacional terms, centrifugal, and coriolis terms that arise from the Langrangian. 
It represents all contributions to the dynamics that do not multiply the accelations. 

$$
\mathbf{Q}(\theta_1,\theta_2,\dot{\theta}_1,\dot{\theta}_2) = \begin{pmatrix}
-(m_1+m_2)gL_1 \sin\theta_1 - m_2 L_1 L_2 \dot{\theta}_2^{\,2} \sin(\theta_1 - \theta_2) \\
m_2 L_1 L_2 \dot{\theta}_1^{\,2} \sin(\theta_1 - \theta_2) - m_2 g L_2 \sin\theta_2 
\end{pmatrix}
$$

Full nonlinear dynamics:

$$
\mathbf{M}(\theta_1 \theta_2)\,\boldsymbol{\alpha} = \mathbf{Q}
$$

Therefore:

$$
\boldsymbol{\alpha} = \mathbf{M}^{-1} \mathbf{Q}
$$

In numerical simulations, this is solved at each timestep using a linear solver.

$$
\begin{pmatrix}
\ddot{\theta}_1 \[4pt] \ddot{\theta}_2\end{pmatrix} =
\mathbf{M}(\theta_1,\theta_2)^{-1}\
\mathbf{Q}(\theta_1,\theta_2,\dot{\theta}_1,\dot{\theta}_2)
$$


Note: do not worry if you do not get it at first. It took me a few weeks understand this nomenclature and why it works. As a comment, the standard manipulator equation is:

$$
Q(\theta) \ddot{\theta} + C(\theta, dot{\theta}) \dot{theta} + G(\theta)
$$

quite use in robotics and multibody dynamics.

## Numerical methods

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

If you want to discover other concerns that we have not been discussed here, you can check the code and change a few lines.

# The fractal revolution

# Bibliography

## Spherical pendulum: spinning around and around

Another classical example of special oscillaroty motion is the spherical pendulum. Why is it? Because we live in a three dimensional world. However, this motion still reprensents a two-degree-of-freedom system.

# Lagrangian formulation:

As always, identifying the system's degree of freedom (in this case angular coordenate $\theta$ (polar angle) and $\phi$ (axial angle)), we can derive the equation of motion using the Euler-Langrangian equations:

<p align="center">
  <img src="figures/Subplots_animation.gif" width="500"><br>
</p>

*Note*: the red unit vector is the $u_{\theta}$.

$$
\ddot{\phi} = -\frac{2  \dot{\theta}  \dot{\phi}}{tan{\theta}}
$$

this equation describes how the axial angle change through time. But, if we look at it with good ayes, we discover a critical point. When is the denominator equal to zero? When the angles is $-\frac{\pi}{2}$ or $\frac{\pi}{2}$. So, when it comes to happen, we have to jump it this issue with an interpolation.

$$
\ddot{\theta} = \frac{\dot{\phi}^2  sin(2  \theta)}{2} - \frac{g}{L} sin(\theta)
$$

These equations describe the motion of the spherical pendulum.

A curious description of this system is about angular momentum conservation. The first equation can be written through Euler-Lagrangian equation:

$$
\frac{\partial{mL^2 sin(\theta)^2 \dot{\phi}}}{\partial{t}} = 0
$$

through this equation we can describe the effective potential like:

$$
mL^2 sin(\theta) ^2 \dot{\phi} = J_z
$$

then:

$$
\dot{\phi} = \frac{J_z}{mL^2 sin(\theta)}
$$

and the second equation it reads as:

$$
\ddot{\theta} = \frac{J_z  cos(\theta)}{m^2 L^4 sin(\theta)^3} - \frac{g}{L} sin(\theta)
$$

These equations describe the full nonlinear dynamics of the system, and are typically solved using numerical methods.

There also exists a way to simplify these equations. Using the small-angle approximation, same masses, and same length. We can transform the nonlinear equation to a couple of:

$$
\ddot{\theta} = \frac{\dot{\phi}^2 2 \theta}{2} - \frac{g}{L} \theta
$$

$$
\ddot{\phi} = {-2 \dot{\theta} \dot{\phi}}{\theta}
$$

<p align="center">
  <img src="figures/Projections.gif" width="500"><br>
</p>

## Numerical methods:

In this chapter (as always) we have to use numerical methods (or a piece of paper and hundreds of hour) to solve the system. In this case, we use: RK45, DOP853 and, our fantastic and extraordinary method, RK4. 

I wanted to bring this RK4 method back to show how badly works in different angles and with different scales. In this way, I want to remind that the 4th Runge-Kutta method could be used in some cases with good results, but when physic goes hard the method becomes inadequate.

<div style="display: flex; flex-direction: column; align-items: center;"> 
  <img src="figures/convergence.png" width="220" style = "margin: 5px 0;">
  <img src="figures/stability.png" width="220" style = "margin: 5px 0;">
  <img src= "figures/time_compute.png" width = "220" style = "margin: 5px 0;">
</div>

## Errors and analysis:
Now we have introduced the numerical methods, and which one is better than the other. We can wonder how well it perfoms.

*Note*: The graphics have been computed with these parameters: $t_max = 150$, ratol = 1e-10, atol = 1e-12, $m = 1.0$ $L = 2.0$ and $dt =0.01$. However, in some error graphics we have to short the axes due to clarity.
For clarity in the right side we display the linearized equation, while in the left side display the normal equation.

<div style="display: flex; flex-direction: column; align-items: center;"> 
  <img src="figures/heat.png" width="330" style = "margin: 5px 0;">
  <img src="figures/heat_linearized.png" width="330" style = "margin: 5px 0;">
</div>

If we want a intersection plane:

<div style="display: flex; flex-direction: column; align-items: center;"> 
  <img src="figures/vertical.png" width="330" style = "margin: 5px 0;">
  <img src="figures/vertical_linearized.png" width="330" style = "margin: 5px 0;">
</div>


As we can see in the linearized graphics, the linearized equations *explode* when inital condition becomes bigger than twenty degrees or less. 

Once we discover when the method lose energy, and how it works through different initial conditions, we can talk about two underlying motion: nutation and precession. 

**Nutation:** rocking, swaying, or nodding motion in the axis of rotation of an object, such as gyroscope.

**Precession:** continuous change in the orientation of a rotating body's axis

These two motions probabibly sound familiar, due to the fact that these motions, with others, describes the Earth's movement. The energy is exchanged between *rocking* and *rotation*.

The spherical pendulum is a combination of: 
- Oscillatory motion in $\theta$ (go up and down, nutation)
- Rotation motion in $\phi$ (spinning aroung the z-axis, precession)

The conserved angular momentum creates the precession through this equation:

$$
mL^2 sin(\theta) ^2 \dot{\phi} = J_z
$$

We can see how if $\theta$ decreases then $\dot{\phi}$ increases;conversely, if $\theta$ increases then $\dot{\phi}$ decreases. So, the angular velocity depends on inclination of the pendulum. In the first figure we can clearly see how, as the angle $\theta$ decreases, the trajectory in the phace space becomes narrow and speeds up.

<div style="display: flex; flex-direction: column; align-items: center;"> 
  <img src="figures/nutation_preccesion.png" width="330" style = "margin: 5px 0;">
  <img src= "figures/phase_space.png" width = "330 " style = "margin: 5px 0;">
</div>

Futhermore, the second equation creates the nutation:

$$
\ddot{\theta} = \frac{J_z  cos(\theta)}{m^2 L^4 sin(\theta)^3} - \frac{g}{L} sin(\theta)
$$

This equation have two terms. The first one is called centrifugal - it pushes outward. While the second one is called gravitational - it pushes downward. In the second picture, the nutation can be appreciated when $\theta$ oscillates while the angle $\phi$ increases in a linear way with some variation.

The thrid figure indicates a continuous probabiliy distribution, where - at the same time as in second figure - the closed curves indicate nutation (modulating the amplitude) and the pattern repetition indicates the precession (modulating the phase). Also, in this figure we can appreciate diffetent densities - places where the pendulum takes a long time - that it coincide with maximus of nutation and slow preccesion.

  <img src="figures/kde_phase_space.png" width="400" style = "margin: 5px 0;">

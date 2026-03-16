<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Damping Vibration</title>
<style>
body{
    font-family: Arial, Helvetica, sans-serif;
    line-height:1.6;
    margin:40px;
}
h1,h2,h3{
    color:#333;
}
code{
    background:#f4f4f4;
    padding:3px 6px;
}
img{
    max-width:700px;
    display:block;
    margin:20px auto;
}
</style>
</head>

<body>

<h1>Damping Vibration</h1>

<p>Hello there,</p>

<p>
Let us talk about the next section of oscillatory systems: <b>Damping Vibration</b>.
We know that in real life motion is not infinite. There exists a concept called 
<b>energy dissipation</b>. Effects such as air friction, internal material friction, 
and many other mechanisms gradually remove energy from the system. As a result, 
every oscillatory motion eventually stops.
</p>

<p>
But how can we represent this kind of motion mathematically?
</p>

<p>
As usual, we assume that if the friction forces are relatively small, the motion
remains approximately periodic while its amplitude slowly decreases with time.
Under this assumption, we consider the same systems studied in the previous section:
</p>

<ul>
<li>Mass–spring system</li>
<li>Simple pendulum</li>
</ul>

<p>
The only additional ingredient is the presence of a <b>dissipative force</b>:
</p>

<p style="text-align:center;">
<b>F<sub>d</sub> = −γ v</b>
</p>

<p>
This force is non-conservative and acts opposite to the direction of motion,
therefore the work performed by this force is negative and removes energy
from the system.
</p>

<h2>Equations of Motion</h2>

<p>
Applying Newton's second law and establishing the dynamic equilibrium,
the equations of motion become:
</p>

<h3>Mass–Spring System</h3>

<p style="text-align:center;">
m x'' + γ x' + kx = 0
</p>

<h3>Pendulum System</h3>

<p style="text-align:center;">
θ'' + (γ/m) θ' + (g/L) sin(θ) = 0
</p>

<img src="figures/pendulum.gif">
<img src="figures/spring.gif">

<h2>Important Parameters</h2>

<p>
Before discussing how to solve these differential equations, it is useful
to define two important parameters:
</p>

<ul>
<li><b>Damping parameter</b>: β</li>
<li><b>Natural angular frequency</b>: ω<sub>0</sub></li>
</ul>

<p>
For the spring system:
</p>

<p style="text-align:center;">
β = γ / (2m)
</p>

<p style="text-align:center;">
ω<sub>0</sub>² = k / m
</p>

<p>
Using these parameters the equation becomes:
</p>

<p style="text-align:center;">
x'' + 2βx' + ω<sub>0</sub>² x = 0
</p>

<p>
These parameters determine the dynamical regime of the system:
</p>

<ul>
<li>If β &lt; ω<sub>0</sub> → <b>Underdamped motion</b></li>
<li>If β = ω<sub>0</sub> → <b>Critical damping</b></li>
<li>If β &gt; ω<sub>0</sub> → <b>Overdamped motion</b></li>
</ul>

<img src="figures/trajectory_spring.png">

<h2>Types of Motion</h2>

<h3>Underdamped Motion</h3>

<p>
In this regime the system still oscillates, but the amplitude decays
exponentially with time. The characteristic time that describes this
decay is called the <b>relaxation time</b>.
</p>

<p>
For the linear equation of the spring system the analytical solution is:
</p>

<p style="text-align:center;">
x(t) = A e<sup>-βt</sup> cos(ωt + φ)
</p>

<p>
This analytical solution exists only for the linear equation.
For nonlinear systems such as the pendulum, numerical methods must be used.
</p>

<h2>Numerical Solutions</h2>

<img src="figures/regime_summary_pendulum.png">

<p>
The previous figure shows the different regimes for the pendulum system,
each solved using different numerical methods.
Since the pendulum equation is nonlinear, we cannot use the analytical
solution.
</p>

<p>
We can observe that the Euler method performs poorly, especially for the
first regime. This occurs because Euler's method accumulates large
numerical errors and is not very stable for oscillatory systems.
</p>

<h2>How Do We Solve It?</h2>

<p>
This time we will not discuss the basic numerical methods in detail,
since they were already explained in the previous section
(Euler method and Runge–Kutta RK4).
</p>

<p>
However, another numerical method is introduced here:
<b>Crank–Nicolson method</b>.
</p>

<p>
The Crank–Nicolson method was developed by
John Crank and Phyllis Nicolson in 1947. It is an implicit
finite-difference method commonly used for solving differential
equations in numerical analysis.
</p>

<p>
It is based on averaging the explicit and implicit Euler methods,
which provides improved stability and second-order accuracy.
</p>

<h3>Crank–Nicolson Scheme</h3>

<p>
The following derivation corresponds to the notes used in the implementation:
</p>

<img src="figures/crank_nicolson_notes.png">

<h2>Energy Validation</h2>

<p>
How do we know if the numerical solution is correct?
</p>

<p>
One way is to analyze the system energy and apply the
<b>work–energy theorem</b>. The dissipative force performs negative
work, therefore the total mechanical energy must decrease with time.
</p>

<img src="figures/energy_initial_angle.png">

<p>
This figure also shows how poorly the Euler method performs.
Because it is not energy-consistent, it can artificially destroy
or generate energy.
</p>

<h2>Error Analysis</h2>

<p>
The first figure corresponds to the spring system. Since the analytical
solution is known, we can compute the absolute error of the numerical
solutions.
</p>

<img src="figures/numerical_errors_spring.png">

<p>
For the pendulum system we cannot compute analytical errors. Instead,
we analyze the <b>stability</b> and <b>convergence</b> of the methods.
</p>

<img src="figures/convergence_pendulum_rk4.png">
<img src="figures/convergence_pendulum_crank_nicolson.png">

<img src="figures/stability_pendulum_crank_nicolson.png">
<img src="figures/stability_pendulum_rk4.png">

<p>
From these results we observe that <b>RK4</b> generally provides higher
accuracy and stability. However, the <b>Crank–Nicolson method</b> also
performs very well and remains a strong alternative, especially in
problems where stability is essential.
</p>

</body>
</html>
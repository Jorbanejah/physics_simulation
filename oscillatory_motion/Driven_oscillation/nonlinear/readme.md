<h1>Nonlinear driven oscillation</h1>

<p>
In this section we study the <strong>nonlinear driven damped oscillator</strong> under a periodic external force. 
The governing equation is
</p>

<p>
$$
m L^2 \ddot{q} + \gamma \dot{q} + m g L sin(q) = F_0 \cos(\omega t).
$$
</p>

<p>
Each term corresponds to a real phyisical effect:
  <li> $m$ $L^2$ $\ddot{q}$ rotacional inertiaa x angular acceleration.</li>
  <li> $\gamma$ $\dot{q}$ viscous dampeing torque.</li>
  <li> $mgL$ $\sin(q)$- gravitational restoring torque.</li>
  <li> $F_0$ $\cos(\omega t)$ - external periodic driving torque.</li>

<p>
To simplify notation, we introduce the standard parameters:
</p>

<p>
$$
\alpha = \frac{F_0}{m L^2}, 
\qquad 
\beta = \frac{\gamma}{2m L^2}, 
\qquad
\omega_0^2 = \frac{g}{L},
$$
</p>

<p>
so the equation becomes
</p>

<p>
$$
\ddot{q} + 2\beta \dot{q} + \omega_0^2 sin(q) = \alpha \cos(\omega t).
$$
</p>

![Regime summary](../figures/regime_summary_nonlinear2.0_2.png)

<p>
  This equation have not analytical solution, so we have to dip into numerical methods: Crank Nicolson (CN), RK4 and Verlet integrator. As we can see, in the previous picture, our numerical method match energy system remarkable well. Also, we ilustrate how the drift energy change with different steps integration or 3-dimensional heatmaps graphics sweeping over normalizes values of $\beta$ and $\omega$.
</p>

![Drift energy](../figures/nonlinear_heatmap_Verlet_0.1.png)

![Steps](../figures/drift_energy_vs_dt.png)

<hr>

<h2>Resonance</h2>

<hr>

<h3>Chaotic motion</h3>

Take a look one more time our precious equation: 

<p>
  $$
  \ddot{q} + 2\beta \dot{q} + \omega_0^2 sin(q) = \alpha \cos(\omega t).
  $$
</p>

<p>
  In mathematics, a differential equation is called non-autonomous where the independent variable t does not appear explicitly. We can write the previous equation like that if we make a system of three first-order autonomous equation:
</p>

<p>
  $$
  \dot{q} = u
  \dot{u} = - 2\beta \dot{q} - \omega_0^2 sin(q) + \alpha \cos(\omega t)
  \dot{\phi} = \omega
  $$
</p>

<p>
  The necessary conditions for an autonomous system of differential equations to admit chaotic solutions are:
<ol>
<li>The system must have at least three independent dynamical variables — <i>condition satisfied</i>.</li>
<li>The system must contain at least one nonlinear coupling — <i>condition satisfied</i>.</li>
</ol>

But, what exactly is chaos? To be precise, what we are studying here is called deterministic chaos: chaotic behaviour arising from deterministic equations of motion. What we will named as chaos is: small changes in the initial conditions lead to exponentially diverging in trajectories. This sensitivity to initial condition is often called: <b>Butterfly effect</b>.

How can we illustrate chaos in a dynamical system?

<li>Bifurcation diagram: A plot showing how the long‑term behavior of the system changes as a parameter varies.
Period‑doubling cascades, windows of stability, and sudden transitions reveal the route to chaos</li>
<li>Lyapunov exponent: A quantitative measure of sensitivity to initial conditions.
A positive Lyapunov exponent indicates exponential divergence of nearby trajectories — the hallmark of chaos.</li>

[Bifurcation + Lyapunov](../figures/Bifurcation_and_lyapunov.png)

<li>Poincare sections:A stroboscopic “slice” of the phase space.
Regular motion produces smooth curves; chaotic motion fills regions irregularly, revealing the underlying structure of the attractor.</li>

[Poincare sections](..figures/figures/Poincare_sections_and_trajectories.png)
Note: these three graphics have an animation in the first readme.md. Moreover, a further explanation will be soon.
</p>


<hr>

<h4>Bibliography</h4>

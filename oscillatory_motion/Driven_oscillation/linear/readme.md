<h1>Linear Driven Oscillation</h1>

<p>
In this section we study the <strong>linear driven damped oscillator</strong> under a periodic external force. 
The governing equation is
</p>

<p>
$$
m \ddot{q} + \gamma \dot{q} + k q = F_0 \cos(\omega t).
$$
</p>

<p>
To simplify notation, we introduce the standard parameters:
</p>

<p>
$$
\alpha = \frac{F_0}{m}, 
\qquad 
\beta = \frac{\gamma}{2m}, 
\qquad
\omega_0^2 = \frac{k}{m},
$$
</p>

<p>
so the equation becomes
</p>

<p>
$$
\ddot{q} + 2\beta \dot{q} + \omega_0^2 q = \alpha \cos(\omega t).
$$
</p>

<hr>

<h2>General Solution</h2>

<p>
This is a second‑order linear non‑homogeneous ODE. Its solution is the sum of:
</p>

<p>
$$
q(t) = q_{\text{hom}}(t) + q_{\text{part}}(t).
$$
  
</p>

<h3>Homogeneous (transient-state) Solution</h3>

<p>
The homogeneous equation is
</p>

<p>
$$
\ddot{q} + 2\beta \dot{q} + \omega_0^2 q = 0.
$$
</p>

<p>
The characteristic equation
</p>

<p>
$$
r^2 + 2\beta r + \omega_0^2 = 0
$$
</p>

<p>
has roots
</p>

<p>
$$
r = -\beta \pm \sqrt{\beta^2 - \omega_0^2}.
$$
</p>

<p>
For the <strong>underdamped</strong> case, $\beta < \omega_0$, the solution becomes
</p>

<p>
$$
q_{\text{hom}}(t) = e^{-\beta t}
\left(A\cos(\omega_d t) + B\sin(\omega_d t)\right),
$$
</p>

<p>
where the <strong>damped natural frequency</strong> is
</p>

<p>
$$
\omega_d = \sqrt{\omega_0^2 - \beta^2}.
$$
</p>

<h3>Particular (Steady‑State) Solution</h3>

<p>
We assume a solution of the form
</p>

<p>
$$
q_{\text{part}}(t) = A \cos(\omega t - \delta).
$$
</p>

<p>
We substitute the trial solution into the differential equation. First compute de derivatives, then substitute into the equation of motion:
</p>

<p>
  $$
  - A \omega cos(\omega t -\delta) + 2 \beta A \omega sin(\omega t + \delta) + \omega_0^2 A cos(\omega  t - \delta) = \alpha cos(\omega t)
  $$
</p>

<p>
  We group cosine and sine terms and then, using the trigonometric identities, we can collecting coefficients of $cos(\omega t)$ and $sin(\omega t)$, obtaining the system.
</p>

<p>
  $$
  A [(\omega_0^2 - \omega^2)cos(\delta) - 2 \beta \omega sin(\delta)] = \alpha \\
  A [(\omega_0^2 - \omega^2)sin(\delta) +  2 \beta \omega cos(\delta)] = 0
  $$
</p>

<p>
  With the second equation, dividing by $cos(\delta)$, you can solve the phase lag:
</p>

<p>
$$
\tan\delta = \frac{2\beta\omega}{\omega_0^2 - \omega^2}.
$$
</p>

<p>
  To solve the amplitude, we have to eliminate the phase, we square and add both equation. Solving for the amplitude yields:
</p>

<p>
$$
A = \frac{\alpha}{\sqrt{(\omega_0^2 - \omega^2)^2 + (2\beta\omega)^2}},
$$
</p>

<p>
This steady‑state term dominates after transients decay.
</p>

<p>
  So the complete equation for a driven damped harmonic oscillator is:
</p>

<p>
  $$
  q (t) = e^{-\beta t} \left(A\cos(\omega_d t) + B\sin(\omega_d t)\right) + A \cos(\omega t - \delta)
  $$
</p>

<hr>

<h2>Verlet Method</h2>

<p>
As I promise when I start this oscillation folder, in every motion I teach you a new numerical method. To compare this time our known numerical methods, we introduce the <strong>velocity Verlet</strong> integrator, commonly used for conservative systems. Starting from Newton’s law
</p>

<p>
$$
m\ddot{q} = F(q),
$$
</p>

<p>
the update equations are
</p>

<p>
$$
q_{n+1} = q_n + v_n \Delta t + \frac{1}{2} a_n (\Delta t)^2,
$$
</p>

<p>
$$
a_{n+1} = \frac{F(q_{n+1})}{m},
$$
</p>

<p>
$$
v_{n+1} = v_n + \frac{1}{2}(a_n + a_{n+1})\Delta t.
$$
</p>

<p>
Although Verlet is designed for conservative systems, it performs remarkably well in this linear driven case, showing excellent energy behavior as well as RK4 and Crank–Nicolson.
</p>

![Regime summary γ = 2, F₀ = 1](../figures/regime_summary2.0_1.png)

<p>
  As we can see, the numerical methods match the analytical solution remarkably well. In particular, the total energy of the system (blue curve in the thrid row) remains essentially constant, demonstrating the excellent stability of the integrators.
  For a deeper analysis of numerical errors, several 3-dimensional heatmaps are included in figures folder. These explore how the energy drift varies when sweeping over different normalized values of $\omega$ and $\beta$, allowing to us to identify regions of stability (far from resonance) and instability (close to resonance) across the parameter space.

![Heatmap_rk4](../figures/heatmap_errorsVerlet_0.png)

  Before discussing resonance, we must understand how energy flows in a non-conservative system. 
  Of course, our mechanical energy is $E_mec = T + U = \frac{1}{2}m dq^2 + \frac{1}{2} k q^2$ - Note: the U's formula depends on potential. However, because of damping and external parameter forcing, the system is not conservative, so $\Delta E_mec \neq 0$.
  
How will we deal with? With the Work-Energy theorem states:
</p>

<p>
  $$
  \frac{d E_{mec}}{dt} = \sum P_{non-conservative} = P_{drive} - P_{diss},
  $$
</p>

<p>
As we know, $P = F v$ or $P = F dq$. While dissipative force is $F_{diss}= - \gamma dq$, drive force is $F_{ext} = F_0 cos(\omega t)$. So, $P_{drive} =  F_{ext} dq$ is the instantaneous power injected by the external force, and $P_{diss} = \gamma dq^2$ is the instantaneous power lost to damping.
</p>
<p>
Integrating over the time gives the work contributions. On the one hand, dissipative energy: $W_{diss} = \int \gamma dq^2\ dt$. On the othen hand, drive energy: $W_{drive}= \int F_{ext} dq\ dt$.
So, the total energy system is: 
</p>

<p>
  $$
  E = E_{mec} + W_{diss} - W{drive}
  $$
</p>
   
<h2>Resonance</h2>

<p>
Resonance occurs when the driving frequency ($\omega$) approaches the natural frequency ($\omega_0$) of the system as we can see in the following picture:
</p>

![Resonance](../figures/beta_vs_power.png)

<p>
In the long‑term, the motion is purely steady‑state because of exponential term:
</p>

<p>
$$
q_{\text{steady}}(t) = A \cos(\omega t - \delta),
$$
</p>

<p>
with amplitude:
</p>

<p>
$$
A = \frac{\alpha}{\sqrt{(\omega_0^2 - \omega^2)^2 + (2\beta\omega)^2}}.
$$
</p>


![Resonance](../figures/beta_vs_amplitude.png)

<hr>

<h3>Near Resonance</h3>

<p>
If
</p>

<p>
$$
4\beta^2 \omega^2 \ll (\omega_0^2 - \omega^2)^2,
$$
</p>

<p>
and $\omega \approx \omega_0$, the denominator becomes small and the amplitude grows large, producing the classical resonance peak as we have showed in the previous pictures.
</p>

<p>
Two equivalent ways to observe resonance:
   
  1_ Tune the natural frequency $\omega_0$ close to the driving frequency $\omega$:
</p>
<p>
  $$
    A \approx \frac{\alpha}{2\beta\omega_0}.
  $$
</p>

<p>
  2_  Tune the driving frequency $\omega$ close to $\omega_0$ .  
    To find the frequency that maximizes the steady‑state amplitude, we minimize the denominator:
</p>

<p>
  $$
   D(\omega) = (\omega_0^2 - \omega^2)^2 + (2\beta\omega)^2.
  $$
</p>

  <p>
    Differentiating and setting  $\frac{d D}{d \omega} = 0$ gives the resonance condition
  </p>
        
<p>
  $$
  \omega_{\text{res}} = \sqrt{\omega_0^2 - 2\beta^2}.
  $$
</p>

<p>
Thus the strongest response occurs when the driving frequency is slightly below the natural frequency, with the shift depending on the damping.
</p>    

<p>
Both approaches reveal the characteristic resonance behavior.
</p>

<hr>

<h3>Bibliography</h3>

<p>
The following resources were consulted during the preparation of this section:
</p>

<a href="https://ocw.mit.edu/courses/res-8-009-introduction-to-oscillations-and-waves-summer-2017/mitres_8_009su17_lec5.pdf">MIT OCW – Driven Oscillations (Lecture 5)</a>

Class notebooks and lecture notes

Physics notes from <a href="https://github.com/RubenCarrionCastro">https://github.com/RubenCarrionCastro</a>

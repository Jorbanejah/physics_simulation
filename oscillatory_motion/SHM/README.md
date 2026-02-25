<h1>Simple Harmonic Motion (SHM) — Numerical Exploration</h1>

<p>Hi, how is it going on?</p>

<p>
This repository is the <strong>first part of a three‑part series</strong> dedicated to understanding oscillatory motion.
Oscillations appear everywhere: from the periodic sway of a leaf to the vibrations inside machines.
But how can we <em>describe</em> this kind of motion?
</p>

<hr>

<h2>📌 1. The Differential Equation of SHM</h2>

<p>A system undergoing simple harmonic motion satisfies the linear second‑order differential equation:</p>

<p>
$$
\frac{d^2x}{dt^2} + \omega^2 x = 0
$$
</p>

<p>whose analytical solution is:</p>

<p>
$$
x(t) = A \cos(\omega t + \phi)
$$
</p>

<p>This equation describes, for example:</p>
<ul>
  <li>A <strong>1‑D mass–spring system</strong></li>
  <li>A <strong>pendulum under the small‑angle approximation</strong> $(\sin\theta \approx \theta)$</li>
</ul>

<hr>

<h2>📌 2. The Real Pendulum Equation</h2>
4
<p>The exact equation of motion for a pendulum of length $L$ is:</p>

<p>
$$
\frac{d^2\theta}{dt^2} + \frac{g}{L}\sin\theta = 0
$$
</p>

<p>This equation is <strong>non‑linear</strong>, and therefore <strong>has no closed‑form analytical solution</strong> in elementary functions.</p>

<hr>

<h2>📌 3. Numerical Methods: Euler vs. Runge–Kutta</h2>

<p>
Euler’s method initially seems to work, but when you compute the <strong>mechanical energy</strong>, you will notice that the motion is <em>not</em> conservative.
</p>

<p>This happens because:</p>
<ul>
  <li>Euler’s method accumulates error at every step.</li>
  <li>For oscillatory systems, this accumulated error grows without bound.</li>
</ul>

<p>A better alternative is the <strong>Runge–Kutta 4th order method (RK4)</strong>, which conserves energy much more accurately.</p>

<hr>

<h2>📌 4. Repository Structure</h2>

<p>This repository contains two main classes:</p>
<ul>
  <li><strong>Pendulum</strong></li>
  <li><strong>Spring</strong></li>
</ul>

<hr>

<h2>📌 5. Pendulum Class</h2>

<p>
The <code>Pendulum</code> class can generate a <strong>2‑dimensional animation</strong> when <code>animate=True</code>.
</p>

<ul>
  <li>If <code>approx=True</code>, it uses the <strong>small‑angle approximation</strong>.</li>
  <li>If <code>approx=False</code>, it solves the <strong>non‑linear equation</strong> using <strong>RK4</strong>.</li>
</ul>

<h3>🔹 Small‑angle approximation animation</h3>
<img src="./figures/pendulum_approx.gif" width="400">

<h3>🔹 Full non‑linear pendulum animation</h3>
<img src="./figures/pendulum.gif" width="400">

<hr>

<h2>📌 6. When Does the Approximation Fail?</h2>

<p>The following figure compares the <strong>energy error</strong> between the approximate and exact models:</p>
<img src="./figures/energy_error_pendulum.png" width="500">

<p>Comparison of the <strong>period</strong> as a function of the initial angle:</p>
<img src="./figures/periods_pendulum.png" width="500">

<hr>

<h2>📌 7. Phase Space Representation</h2>

<p>Phase space provides a powerful way to visualize the dynamics of oscillatory systems.</p>

<img src="./figures/phase_space_pendulum.png" width="500">

<p>A colormap of trajectories as a function of the initial angle:</p>

<img src="./figures/colormap_pendulum.png" width="500">

<p>
For small angles, the trajectory is nearly circular.
For large angles, the trajectory becomes elliptical due to the non‑linear term.
</p>

<hr>

<h2>📌 8. Spring Class</h2>

<p>
The <code>Spring</code> class models a <strong>1‑D mass–spring system</strong> located at \((x_0, y_0)\).
It computes:
</p>

<ul>
  <li>Position</li>
  <li>Velocity</li>
  <li>Mechanical energy</li>
  <li>Phase space trajectory</li>
</ul>

<p>and can animate the motion:</p>

<img src="./figures/spring.gif" width="400">

<hr>

<h2>📌 9. Influence of the Spring Constant \(k\)</h2>

<p>The period of a mass–spring system is:</p>

<p>
$$
T = 2\pi\sqrt{\frac{m}{k}}
$$
</p>

<p>Thus, increasing \(k\) <strong>reduces the period</strong>, as shown here:</p>

<img src="./figures/periods_spring.png" width="500">

<p>The phase space also changes with \(k\):</p>

<img src="./figures/phase_space_spring.png" width="500">

<p>
As \(k\) increases:
</p>

<ul>
  <li>The system oscillates faster.</li>
  <li>The phase‑space ellipse becomes <strong>narrower in position</strong> and <strong>wider in velocity</strong>.</li>
</ul>

<hr>

<h2>📌 10. What’s Next?</h2>

<p>
This is only the beginning.
In the next part, we will explore:
</p>

<ul>
  <li><strong>Damped oscillations</strong></li>
  <li><strong>Driven oscillations</strong></li>
  <li><strong>Resonance</strong></li>
  <li><strong>Energy dissipation</strong></li>
</ul>

<p>Stay tuned!</p>

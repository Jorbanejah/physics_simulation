<h2>Parabolic motion</h2>

<p>
Anyone who has ever studied basic kinematics has heard about <em>parabolic motion</em>. It is the result of throwing a projectile with initial position
$(x_0, y_0)$, initial speed $(v_0)$, and launch angle $(\alpha)$. The equations of motion along the x and y axes in the ideal case are:
</p>

<p>
Horizontal velocity (no horizontal acceleration):
</p>

$$
v_x(t) = v_{0x} = v_0 \cos\alpha
$$

<p>
Vertical velocity (affected by gravity):
</p>

$$
v_y(t) = v_{0y} - g t = v_0 \sin\alpha - g t
$$

<p>
Integrating these expressions, we obtain the position:
</p>

$$
x(t) = x_0 + v_0 \cos\alpha t
$$

$$
y(t) = y_0 + v_0 \sin\alpha t - \frac{1}{2} g t^2
$$

<p>
This produces the familiar parabolic trajectory, which can be represented as:
</p>

<table>
  <tr>
    <td align="center">
      <img src="projectile_motion_analytic.png" alt="Projectile motion without drag" width="300"/><br/>
      <strong>Projectile motion without air resistance</strong>
    </td>
  </tr>
</table>

<p>
Too boring, isn’t it?
</p>

<p>
Well, let’s make the system more interesting.
</p>

<hr>

<h2>🧠 Resistance coefficient (linear drag)</h2>

<p>
Imagine you are walking peacefully down the street, throwing a ball above your head, when suddenly a strong gust of wind bends its trajectory.
You still manage to catch it because you are faster, but you start wondering if there is a physical model that can describe this motion.
There is: we introduce a <strong>linear drag force</strong>.
</p>

<p>
In this case, the acceleration depends on the velocity:
</p>

$$
m a_x = -k m v_x
$$

$$
m a_y = -k m v_y - m g
$$

<p>
where k is the resistance coefficient. Dividing by m, we obtain:
</p>

$$
a_x = \frac{dv_x}{dt} = -k v_x
$$

$$
a_y = \frac{dv_y}{dt} = -k v_y - g
$$

<p>
We can solve these equations both analytically and numerically. We will do both.
</p>

<hr>

<h3>Analytical solution with linear drag</h3>

<h4>Horizontal motion</h4>

<p>
We solve
</p>

$$
\frac{dv_x}{dt} = -k v_x
$$

<p>
This is a first-order linear differential equation with solution:
</p>

$$
v_x(t) = v_{0x} e^{-k t}
$$

<p>
Integrating once more:
</p>

$$
x(t) = x_0 + \int_0^t v_{0x} e^{-k s} ds
     = x_0 + \frac{v_{0x}}{k} \left(1 - e^{-k t}\right)
$$

<h4>Vertical motion</h4>

<p>
Now we solve
</p>

$$
\frac{dv_y}{dt} + k v_y = -g
$$

<p>
The solution is:
</p>

$$
v_y(t) = \left(v_{0y} + \frac{g}{k}\right) e^{-k t} - \frac{g}{k}
$$

<p>
Integrating to obtain the vertical position:
</p>

$$
y(t) = y_0 + \int_0^t \left[\left(v_{0y} + \frac{g}{k}\right) e^{-k s} - \frac{g}{k}\right] ds
$$

$$
y(t) = y_0 + \frac{v_{0y} + \frac{g}{k}}{k} \left(1 - e^{-k t}\right) - \frac{g}{k} t
$$

<p>
These are the analytical equations of motion with linear air resistance.
</p>

<hr>

<h3>Maximum height</h3>

<p>
The maximum height is reached when the vertical velocity becomes zero:
</p>

$$
v_y(t_h) = 0
$$

$$
\left(v_{0y} + \frac{g}{k}\right) e^{-k t_h} - \frac{g}{k} = 0
$$

<p>
Solving for t_h:
</p>

$$
e^{-k t_h} = \frac{\frac{g}{k}}{v_{0y} + \frac{g}{k}} = \frac{g}{k v_{0y} + g}
$$

$$
t_h = \frac{1}{k} \ln\left(\frac{k v_{0y} + g}{g}\right)
$$

<p>
Substituting this time into y(t), we obtain the maximum height $H_{\max}$:
</p>

$$
H_{\max} = y(t_h)
= y_0 + \frac{v_{0y} + \frac{g}{k}}{k} \left(1 - e^{-k t_h}\right) - \frac{g}{k} t_h
$$

<p>
This expression can be simplified further, but in this form it already gives the exact maximum height with linear drag.
</p>

<hr>

<h3>Range and transcendental equation for the flight time</h3>

<p>
To compute the range, we need the total flight time T, i.e, the time when the projectile returns to the ground:
</p>

$$
y(T) = 0
$$

<p>
Assuming y_0 = 0, we have:
</p>

$$
0 = \frac{v_{0y} + \frac{g}{k}}{k} \left(1 - e^{-k T}\right) - \frac{g}{k} T
$$

<p>
Multiplying by k:
</p>

$$
\left(v_{0y} + \frac{g}{k}\right) \left(1 - e^{-k T}\right) - g T = 0
$$

<p>
Rewriting:
</p>

$$
T = \left(\frac{v_{0y}}{g} + \frac{1}{k}\right) \left(1 - e^{-k T}\right)
= \frac{k v_{0y} + g}{g k} \left(1 - e^{-k T}\right)
$$

<p>
This is a <strong>transcendental equation</strong> in T, so it cannot be solved analytically in closed form. However, we can approximate it.
</p>

<hr>

<h3>Approximation using a Taylor expansion</h3>

<p>
We expand the exponential in a Taylor series up to third order:
</p>

$$
e^{-k T} \approx 1 - k T + \frac{(k T)^2}{2} - \frac{(k T)^3}{6}
$$

<p>
Substituting into the transcendental equation:
</p>

$$
T = \left(\frac{v_{0y}}{g} + \frac{1}{k}\right)
\left[1 - \left(1 - k T + \frac{(k T)^2}{2} - \frac{(k T)^3}{6}\right)\right]
$$

$$
T = \left(\frac{v_{0y}}{g} + \frac{1}{k}\right)
\left(k T - \frac{(k T)^2}{2} + \frac{(k T)^3}{6}\right)
$$

<p>
Keeping terms up to first order in k, we obtain the approximate flight time:
</p>

$$
T \approx \frac{2 v_{0y}}{g} \left(1 - \frac{k v_{0y}}{3 g}\right)
$$

<p>
Once we have the time, we can approximate the range R. The horizontal position with drag is:
</p>

$$
x(t) = x_0 + \frac{v_{0x}}{k} \left(1 - e^{-k t}\right)
$$

<p>
Expanding again for small k, we can write:
</p>

$$
R \approx v_{0x} \left(T - \frac{1}{2} k T^2\right)
$$

<p>
Substituting the approximate T into this expression gives an approximate range with drag.
</p>

<p>
On the other hand, without air resistance, the well-known range is:
</p>

$$
R_0 = \frac{v_0^2}{g} \sin(2\alpha)
$$

<p>
Comparing both, we can write an approximate correction:
</p>

$$
R' \approx R_0 \left(1 - \frac{4 k v_{0y}}{3 g}\right)
$$

<p>
So, the error up to first order in k is:
</p>

$$
\Delta R = R_0 - R' \approx R_0 \frac{4 k v_{0y}}{3 g}
$$

<hr>

<h3>Numerical solution using Euler’s method</h3>

<p>
We can also solve the problem numerically using Euler’s method. We rewrite the system as:
</p>

$$
\frac{dx}{dt} = v_x, \quad \frac{dy}{dt} = v_y
$$

$$
\frac{dv_x}{dt} = -k v_x
$$

$$
\frac{dv_y}{dt} = -k v_y - g
$$

<p>
Using a time step $(\Delta t)$, the Euler update equations are:
</p>

$$
v_x^{n+1} = v_x^n - k v_x^n \Delta t
$$

$$
v_y^{n+1} = v_y^n - (k v_y^n + g) \Delta t
$$

$$
x^{n+1} = x^n + v_x^n \Delta t
$$

$$
y^{n+1} = y^n + v_y^n \Delta t
$$

<p>
With these numerical equations, the problem becomes straightforward to implement in Python.
</p>

<table>
  <tr>
    <td align="center">
      <img src="energy.png" alt="Energy evolution" width="300"/><br/>
      <strong>Energy evolution with drag</strong>
    </td>
    <td align="center">
      <img src="range.png" alt="Range comparison" width="300"/><br/>
      <strong>Range: analytical vs approximate</strong>
    </td>
  </tr>
</table>

<hr>

<h3>Comparing analytical and numerical trajectories</h3>

<p>
Now that we have both the analytical and numerical equations of motion, we can compare how the trajectory changes for different values of \(k\),
and how much the numerical solution differs from the analytical one.
</p>

<div align="center">
  <img src="projectile_motion_analytic.png" width="400" alt="Projectile Motion Animation">
  <p><strong>Projectile motion with drag – animation</strong></p>
</div>

<div align="center">
  <img src="projectile_motion_comparison.png" width="400" alt="Projectile Motion Animation">
  <p><strong>Projectile motion with drag – animation</strong></p>
</div>



<p>
Why did we compute two different solution methods instead of just one? First, knowledge never takes up space. Second, because of simulation:
when you want a smooth animation, you want to avoid gaps, inconsistencies, and irregularities. We want a continuous and visually coherent motion.
Since the animation shows how the trajectory changes between different media and different values of k, we would have obtained irregular
trajectories if we had not used a numerical method.
</p>

<div align="center">
  <img src="projectile_motion.gif" width="600" alt="Projectile Motion Animation">
  <p><strong>Projectile motion with drag – animation</strong></p>
</div>


<h2>🧠 Bounce + Resistance coefficient</h2>

<p>
What if instead of catching the ball we had failed? We would have seen how the ball rebounds on the street, performing smaller and smaller
parabolic trajectories until it finally loses all its energy (the system is not perfect, sorry). How can we describe this behaviour?
</p>

<p>
Maintaining the same drag model as before, the only additional ingredient we need is the <strong>coefficient of restitution</strong> e.
When the projectile hits the ground, we compute the vertical velocity just before impact, $- v_y$. The rebound velocity is simply:
</p>

$$
v_y^{new} = - e v_y^{old}
$$

<p>
The horizontal velocity keeps evolving under drag, and the vertical velocity is reset according to the expression above. With each bounce,
the energy decreases, producing progressively smaller trajectories.
</p>

<div align="center" style="margin-top:20px;">
  <img src="projectile_motion_rebound_comparison.png" width="450" alt="Bounce with drag"><br>
  <strong>Rebound under linear drag</strong>
</div>

<div align="center" style="margin-top:25px;">
  <img src="bounce_animation.gif" width="550" alt="Bounce animation with drag"><br>
  <strong>Animation of multiple rebounds with drag</strong>
</div>



## Linear Driven Oscillation

First, we consider the **linear driven, damped oscillator**.

## Equation

<p>
$$
m\ddot{q} +  \gamma \dot{q} + k q = F0 \cos(\omega t)
$$
</p>

We introduce the standard parameters

<p>
$$
\alpha = \frac{F_0}{m}, \qquad 
\beta = \frac{\gamma}{2m}, \qquad
\omega_0^2 = \frac{k}{m},
$$
</p>

so the equation becomes

<p>
$$
\ddot{q} + 2\beta \dot{q} + \omega_0^2 q = \alpha \cos(\omega t).
$$
</p>

This is a **second‑order linear non‑homogeneous ODE**, so the general solution is

<p>
$$
q(t) = q_{\text{hom}}(t) + q_{\text{part}}(t).
$$
</p>

---

## Homogeneous Solution

The homogeneous equation is

<p>
$$
\ddot{q} + 2\beta \dot{q} + \omega_0^2 q = 0.
$$
</p>

The characteristic equation

<p>
$$
r^2 + 2\beta r + \omega_0^2 = 0
$$
</p>

has roots

<p>
$$
r = -\beta \pm \sqrt{\beta^2 - \omega_0^2}.
$$
</p>

For the underdamped case $\beta$ < $\omega_0$, the solution is

<p>
$$
q_{\text{hom}}(t) = e^{-\beta t}
\left(A\cos(\omega_d t) + B\sin(\omega_d t)\right),
$$
</p>

where the **damped natural frequency** is

<p>
$$
\omega_d = \sqrt{\omega_0^2 - \beta^2}.
$$
</p>

---

## Particular (Steady‑State) Solution

We assume a solution of the form

<p>
$$
q_{\text{part}}(t) = C \cos(\omega t - \delta).
$$
</p>

Solving for the amplitude C and phase lag $\delta$ gives

<p>
$$
C = \frac{\alpha}{\sqrt{(\omega_0^2 - \omega^2)^2 + (2\beta\omega)^2}},
$$
</p>

<p>
$$
\tan\delta = \frac{2\beta\omega}{\omega_0^2 - \omega^2}.
$$
</p>

---

## Verlet Method

After discussing resonance analytically, we introduce the **Verlet method**, a symplectic integrator used for conservative systems.

Starting from Newton’s second law

<p>
$$
m\ddot{q} = F(q),
$$
</p>

the **velocity Verlet** update equations are

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

As shown in the figure, Verlet performs extremely well compared with RK4 and Crank–Nicolson, both in energy conservation and in reproducing the analytical trajectory.

---

## Resonance

Now we are ready to discuss **resonance**.  
Resonance occurs when transient behavior dies out and only the steady‑state oscillation remains.

The long‑term solution is

<p>
$$
q_{\text{steady}}(t) = A \cos(\omega t - \delta),
$$
</p>

with amplitude

<p>
$$
A = \frac{\alpha}{\sqrt{(\omega_0^2 - \omega^2)^2 + (2\beta\omega)^2}}.
$$
</p>

### What happens near resonance?

If

<p>
$$
4\beta^2 \omega^2 \ll (\omega_0^2 - \omega^2)^2,
$$
</p>

and we choose \(\omega \approx \omega_0\),  
the denominator becomes small and the amplitude grows large.

This produces the classical **resonance peak**.

We can simulate two scenarios:

1. **Tune the natural frequency** $\omega_0$ close to the driving frequency $\omega$. Then:
<p>
$$
A = \frac{F0}{2 \beta \omega_0}
$$

<p>
2. **Tune the driving frequency** $\omega$ close to the natural frequency $\omega_0$.

Both produce the characteristic resonance behavior.

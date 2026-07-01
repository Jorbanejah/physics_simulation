"""
Special graphics:

- Frequency-maps analysis or Frequency Laskar maps.

    - 1_ Write the equation
    - 2_ Choose a grid of initial conditions: 
            1_ Choose an energy level
            2_ Scan over initial conditions angles
    - 3_Integrate the equation long enough times. (DOP853 with t_times = 250)
    - 4_Build a complex signal for frequency extraction -> E.g: Z (t) = theta(t) + j * phi(t)
    - 5_Apply Laskar-style frequency analysis (NAFF): use a Fourier to approximate Z(t), extract the fundamental omega_k and take a frequency vector (v_1, v_2)
    - 6_Construct the frequency map: associate the initial conditions with its frequencies, plot the frequencies vector v_1/ v_2vs the peak parameter with a fancy colormap


    A POWERFUL Laskar map could be: 
        
    - Two-window analysis: compute frequencies over two successive time windows for each trajectory.

    - Frequency drift: measure  Δv = v(1)-v(2); large drifts signal chaotic diffusion.
    
"""



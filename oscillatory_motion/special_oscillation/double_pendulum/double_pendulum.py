'''
The following code describe the motion of the double-pendulum system.
Description of the system: two masses join to a inextensable rod. The system will be need: L_1, L_2 -- string: m1, m2 -- masses
Structure of the code:

    Class Params: where we define the main parameters, on top of that a boolean parameter called: small-angle.

    def equation(params)
    def approx_equation(params)
    Class DoublePendulum: 
        def __init__(self, params, **kwargs):
        def run
        def Transform
        def energies
'''

import numpy as np
from scipy.integrate import solve_ivp
from dataclasses import dataclass

@dataclass
class Params():

    #Main parameter
    m1: float = 1
    m2: float = 1.5
    L1: float = 1
    L2: float = 2

    #Initial condition
    q0: float = np.deg2rad(30)
    dq0: float = 0

    #Times
    t: int = 15
    dt: float = 0.01

    #Extra parameters
    small_angle: bool = False
    Inextansable: bool = False

###
# ----------------Integrator -----------------------
###

def rk4(f, t, y, dt, params=None):
    """
    Generic RK4 integrator

    Parameters
    ----------
    f : function
        f(t, y, params) -> dy/dt
    t : float
    y : np.ndarray
    dt : float
    params : optional

    Returns
    -------
    y_new : np.ndarray
    """

    k1 = f(t, y, params)
    k2 = f(t + dt/2, y + dt/2 * k1, params)
    k3 = f(t + dt/2, y + dt/2 * k2, params)
    k4 = f(t + dt, y + dt * k3, params)

    return y + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)

###
#-----------------equations --------------
###

def equation_double_pendulum(theta, dtheta, p):
    """
    theta = [theta1, theta2]
    dtheta = [dtheta1, dtheta2]
    p = dict with m1, m2, L1, L2
    """
    g = 9.81
    m1, m2 = p["m1"], p["m2"]
    L1, L2 = p["L1"], p["L2"]

    th1, th2 = theta
    w1, w2 = dtheta

    delta = th1 - th2

    den = L1 * (2*m1 + m2 - m2 * np.cos(2*delta))

    num1 = (
        -g*(2*m1 + m2)*np.sin(th1)
        - m2*g*np.sin(th1 - 2*th2)
        - 2*np.sin(delta)*m2*(w2**2 * L2 + w1**2 * L1 * np.cos(delta))
    )

    num2 = (
        2*np.sin(delta) * (
            w1**2 * L1 * (m1 + m2)
            + g*(m1 + m2)*np.cos(th1)
            + w2**2 * L2 * m2 * np.cos(delta)
        )
    )

    dd1 = num1 / den
    dd2 = num2 / (L2 * (2*m1 + m2 - m2*np.cos(2*delta)))

    return dd1, dd2

def equation_approx(theta, dtheta, p):
    """
    theta = [theta1, theta2]
    dtheta = [dtheta1, dtheta2]
    p = dict with m1, m2, L1, L2

    Here we keep de velocity coupling but linearize only the angles

    sin(theta1 - theta2) = theta1 -theta2
    cos(theta1- theta2) = 1
    """
    g = 9.81
    m1, m2 = p["m1"], p["m2"]
    L1, L2 = p["L1"], p["L2"]

    th1, th2 = theta
    w1, w2 = dtheta

    delta = th1 - th2

    dem_theta1 = L1*(2*m2 + m1)
    dem_theta2  = L2 * (2*m2 + m1)

    num_theta1 = -g * (2*m1 + m2)* th1 - m2 * g *(th1 - 2 * th2) - 2*m2(delta) * (w2 ** 2 * L2 + w1 ** 2 * L1)
    num_theta2 = 2 * delta * ((w1 ** 2 *L1 *(m1 + m2) + g *(m1 + m2)  + w2 **2 *L2 *m2))

    return num_theta1 / dem_theta1, num_theta2 / dem_theta2

###
# ------------ Main class ------------
###

class DoublePendulum():
    '''
    Double Pendulum system.

    ------
    Parameters:
        m1 : float
        first mass pendulum.
        m2 : float
        second mass pendulum.
        L1 : float
        length first rod pendulum. 
        L2 : float
        length second rod pendulum.
        q0 : float
        Initial generalize coordenate
        dq0 : float
        Initial generalize velocity
        t:
        total system time (default: 15)
        dt:
        Time step for numerical integration (default: 0.01)

    ------
    Other parameters:

        Inextansable : bool
        a parameter that simulates either inextensable or non-inextensable rod for L1 paramater(default: False)

        small_angle : bool
        parameter that decides which equation use either non-approximation or approximation (default: False)

    ------
    Note:
    This class implements a double pendulum model that can operate in 
    either small-angle approximation or without approximation, depending on
    the selected configuration of small_angle paramter
    '''

    def __init__(self, m1, m2, L1, L2, q0, dq0 , t, dt, **kwargs):

        positive = [m1, m2, L1, L2]

        invalid = [i for i in positive if i < 0]

        if invalid:
            raise ValueError(f"The following parameters must be positive: {invalid}")

        self.m1 = m1
        self.m2 = m2
        self.L1 = L1
        self.L2 = L2

        self.q0 = q0
        self.dq0 = dq0

        self.t = t
        self.dt = dt

        self.params = kwargs
    
    def run(self):

        self.inextansable = self.params.get("Inextansable")
        self.small_angle = self.params.get("small_angle")

        if self.inextansable is False:

            if self.small_angle is False:

                sol = solve_ivp(equation_double_pendulum, t_span= self.t, y0 = ([self.q0, self.dq0]), method = "RK45", dense_output= True)
            
            elif self.small_angle is True:

                sol = solve_ivp(equation_approx, t_span= self.t, y0 = ([self.q0, self.dq0]), method = "RK45", dense_output=True)

            else:

                raise ValueError("The small-angle parameter must be True or False")








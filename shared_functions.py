import numpy as np

def lagrangean_derivative_coefs(dx):
    """
    Returns the coeficients for the Lagrangeand derivative of the differences
    array 'dx'. The first point has a right derivative, last point a left
    derivative and central difference is for the mid-points.
    """
    ldc1 = np.concatenate(([-(2*dx[0]+dx[1])/(dx[0]*(dx[0]+dx[1]))],
                          -dx[1:]/(dx[:-1]*(dx[:-1]+dx[1:])),
                          [dx[-1]/(dx[-2]*(dx[-2]+dx[-1]))]))
    ldc2 = np.concatenate(([dx[0]+dx[1]/(dx[1]*dx[0])],
                          (dx[1:] - dx[:-1])/dx[:-1]/dx[1:],
                          [(dx[-1]+dx[-2])/(dx[-2]*dx[-1])]))
    ldc3 = np.concatenate(([-dx[0]/(dx[1]*(dx[1]+dx[0]))],
                           dx[:-1]/(dx[1:]*(dx[:-1]+dx[1:])),
                           [(2*dx[-1]+dx[-2])/(dx[-1]*(dx[-2]+dx[-1]))]))

    return ldc1, ldc2, ldc3

def h2Kh(h, n, m, gamma, Ks, Kh = None, tmp1 = None, tmp2 = None):
    if dudh and tmp1 and tmp2:
        tmp1[:] = np.power( gamma * h, n-1.)
        tmp2[:] = np.power(1 + gamma * h * tmp1, m/2.)
        Kh[:]   = Ks/tmp2 * np.power(1-tmp1/np.power(tmp2, 2), 2)
    else:
        tmp1  = np.power( gamma * h, n-1.)
        tmp2 = np.power(1 + gamma * h * tmp1, m/2.)
        Kh   = Ks/tmp2 * np.power(1-tmp1/np.power(tmp2, 2), 2)

    return Kh

def dudh(h, n, m, gamma, Ks, dudh = None, tmp1 = None, tmp2 = None):
    if dudh and tmp1 and tmp2:
        tmp1[:] = np.power(gamma*h, n-1)
        tmp2[:] = np.power(1 + gamma*h * tmp1, m+1)
        dudh[:] = - gamma*(n-1) * tmp1 / tmp2
    else:
        tmp1 = np.power(gamma*h, n-1)
        tmp2 = np.power(1 + gamma*h * tmp1, m+1)
        dudh = - gamma*(n-1) * tmp1 / tmp2

    return dudh

def h2u(h, n, m, gamma, u = None):
    if u:
        u[:] = 1/np.power(1+ np.power(gamma * h, n), m)
    else:
        u    = 1/np.power(1+ np.power(gamma * h, n), m)

    return u

def u2h(u, n, m, gamma, h = None):
    if h:
        h[:] =  np.power(np.power(u, -1/m) - 1, 1/n) / gamma
    else:
        h    =  np.power(np.power(u, -1/m) - 1, 1/n) / gamma

    return h

def characteristics(t, u, mass_in, mass_out, s1, s2, model):
    porosity = model.porosity
    y  = model.y
    dy = model.dy
    r0 = model.r0
    L  = model.l
    l0_out = L + model.l0_out
    l_out  = L + model.l0_out - mass_out

    ds = s2 - s1

    omega2g = (np.power(model.omega_start  + (model.omega - model.omega_start)
                                  * (1 - np.exp(-model.omega_gamma*t)), 2)
                / model.g)

    GC = np.empty(t.shape, float)
    RM = np.empty(t.shape, float)

    P = np.pi * model.d / 4

    # Water mass
    wm_sat = (ds/2  * (dy[0]* u[:, 0] + dy[-1]*u[:, -1]
                + np.sum((dy[:-1] + dy[1:])*u[:, 1:-1], 1)))

    WM = model.density * P * (mass_in + model.porosity*wm_sat + mass_out)

    for i in range(len(t)):
        # Gravitational center
        gc_unsat = (porosity * 1/2 * model.density * ds[i]
                    * ((r0 + s1[i])*dy[0]*u[i, 0] + (r0 + s2[i])*dy[-1]*u[i, -1]
                        + np.sum((dy[:-1]+dy[1:])*(r0 + s1[i] + ds[i]*y[1:-1])*u[i, 1:-1])))
        gc_sat   = (1/2 * model.density
                    * (porosity * (np.power(r0 + s1[i], 2) - np.power(r0, 2))
                       + (np.power(r0, 2) - np.power(r0 - mass_in[i], 2))
                       + (np.power(r0 + l0_out, 2) - np.power(r0 + l_out[i], 2))))
        GC[i] =  P * (gc_unsat + gc_sat) / WM[i]

        # Rotational momentum
        rm_unsat = (porosity * 1/4 * model.density * ds[i]
                    * (np.power(r0 + s1[i], 2)*dy[0]*u[i, 0]
                       + np.power(r0 + s2[i], 2)*dy[-1]*u[i, -1]
                       + np.sum((dy[:-1]+dy[1:]) * u[i, 1:-1]
                                 * np.power(r0 + s1[i] + ds[i]*y[1:-1], 2))))
        rm_sat   = (1/6 * model.density
                    * (porosity * (np.power(r0 + s1[i], 3) - np.power(r0, 3))
                        + (np.power(r0, 3) - np.power(r0 - mass_in[i], 3))
                        + (np.power(r0 + l0_out, 3) - np.power(r0 + l_out[i], 3))))

        RM[i] = omega2g[i] * P * (rm_unsat + rm_sat)

    return GC, RM, WM

def f1(t):
    return 1.7032046506 * np.power(t, 1.233644749)
def f2(t):
    return 0.630314472 * np.log(t) + 8.4248850255
def f3(t):
    return 0.1332308098 * np.log(t) + 9.5952480661

def find_omega2g(t, model):
    if model.omega == 10:
        if t > 22:
            omega = model.omega
        elif t > 21:
            omega = (22-t)*f3(t) + (t - 21)*model.omega
        elif t > 12:
            omega = f3(t)
        elif t > 10:
            omega = (12-t)/2*f2(t) + (t - 10)/2*f3(t)
        elif t > 5:
            omega =  f2(t)
        elif t > 4:
            omega = (5-t)*f1(t) + (t-4)*f2(t)
        else:
            omega = f1(t)

        print('t = ', t, ', omega = ', omega*60)
            return np.power(omega, 2)/model.g
    else:
        return (np.power(model.omega_start  + (model.omega - model.omega_start)
                         *(1 - np.exp(-model.omega_gamma*t)), 2)
                / model.g)

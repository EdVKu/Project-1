import numpy as np
def euler(f,f0,n,h):
    t = 0
    b = n*h
    s0 = f0
    s, T = [],[]
    while t <= b:
        ssq = s0 + h*f(t, s0)
        s.append(ssq)
        T.append(t)
        s0 = ssq
        t += h
    return s, T

def rkuta2(f,f0,n,h):
    t = 0
    b = n*h
    s0 = f0
    s,T = [],[]
    while t<=b:
        k1n = h*f(t, s0)
        k2n = h*f(t, s0 + 0.5*k1n)
        s.append(s0 + k2n)
        T.append(t)
        s0 += k2n
        t += h
    return s, T

def rkuta4(f,f0,n,h):
    t = 0
    b = n*h
    s0 = f0
    s,T = [],[]

    while t<=b:
        k1n = h*f(t, s0)
        k2n = h*f(t, s0 + 0.5*k1n)
        k3n = h*f(t, s0 + 0.5*k2n)
        k4n = h*f(t, s0 + k3n)
        si = s0 + k1n/6 + k2n/3 + k3n/3 + k4n/6
        s.append(si)
        T.append(t)
        s0 = si

        t += h
    return s, T

def SEuler(f,g,f0,n,h):
    p, q = [], []
    return p, q

def StormerV(f,g,f0,n,h):
    p, q = [], []

    return p, q

def dynamics_solve(f, f0 = 1, h = 0.1, n = 100,d = 1, t0 = 0, method = 0):
    method = input("please enter a number to choose your method (0 = Euler, 1 = Runge Kutta 2, 2 = Runge Kutta 4)")
    if int(method) == 0:
        vsol = []
        for i in range(D):
            vsol.append(euler(f[i], f0, n, h))
        return np.array(vsol)
    elif int(method) == 1:
        vsol = []
        for i in range(D):
            vsol.append(rkuta2(f[i], f0, n, h))
        return np.array(vsol)
    elif int(method) == 2:
        vsol = []
        for i in range(D):
            vsol.append(rkuta4(f[i], f0, n, h))
        return np.array(vsol)
    else:
        return None
def hamiltonian_solve(f,g, f0 = 1, h = 0.1, n = 100,D = 1, t0 = 0, method = 0):
    if type(f)!= "list":
        f = list(f)
    method = input("please enter a number to choose your method (0 = Euler, 1 = Runge Kutta 2, 2 = Runge Kutta 4, 3 = Symplectic Euler, 4 = Stormer-Verlet)")
    if 0==int(method)<=2:
        dynamics_solve(f, f0, h, n, D, t0, method)
    elif int(method) == 3:
        vsol = []
        for i in range(D):
            vsol.append(SEuler(f[i], f0, n, h))
        return np.array(vsol)
    elif int(method) == 4:
        vsol = []
        for i in range(D):
            vsol.append(StormerV(f[i], f0, n, h))
        return np.array(vsol)
    else:
        return "None"
    
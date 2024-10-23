import numpy as np
def euler(f,f0,n,h):
    t = 0
    b = n*h
    s0 = f0
    s = []
    while t <= b:
        ssq = s0 + h*f(t, s0)
        s.append(ssq)
        s0 = ssq
        t += h
    return s

def rkuta2(f,f0,n,h):
    t = 0
    b = n*h
    s0 = f0
    s = []
    while t<=b:
        k1n = h*f(t, s0)
        k2n = h*f(t, s0 + 0.5*k1n)
        s.append(s0 + k2n)
        s0 += k2n
        t += h
    return s

def rkuta4(f,f0,n,h):
    t = 0
    b = n*h
    s0 = f0
    s = []
    while t<=b:
        k1n = h*f(t, s0)
        k2n = h*f(t, s0 + 0.5*k1n)
        k3n = h*f(t, s0 + 0.5*k2n)
        k4n = h*f(t, s0 + k3n)
        si = s0 + k1n/6 + k2n/3 + k3n/3 + k4n/6
        s.append(si)
        s0 = si
        t += h
    return s

def SEuler(f,f0,n,h):
    return 0

def StormerV(f,f0,n,h):
    return 0

def dynamics_solve(f, f0, h, n):
    method = input("please enter a number to choose your method (0 = Euler, 1 = Runge Kutta 2, 2 = Runge Kutta 4)")
    if int(method) == 0:
        return euler(f, f0, n, h)
    elif int(method) == 1:
        return rkuta2(f, f0, n, h)
    elif int(method) == 2:
        return rkuta4(f, f0, n, h)
    else:
        return None
def hamiltonian_solve(f, f0 ,n ,h):
    method = input("please enter a number to choose your method (0 = Euler, 1 = Runge Kutta 2, 2 = Runge Kutta 4, 3 = Symplectic Euler, 4 = Stormer-Verlet)")
    if int(method) == 0:
        return euler(f, f0, n, h)
    elif int(method) == 1:
        return rkuta2(f, f0, n, h)
    elif int(method) == 2:
        return rkuta4(f, f0, n, h)
    elif int(method) == 3:
        return SEuler(f, f0, n, h)
    elif int(method) == 4:
        return StormerV(f, f0, n, h)
    else:
        return None
a = dynamics_solve(lambda t,x : np.exp(-x*t), 1,5e-5,100)

print(a)
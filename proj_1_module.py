import numpy as np
def euler(f,f0,n,h,t_o):
    t = t_o
    b = h*(n+1)
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
    b = h*(n+1)
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
def T(n,h,t_o = 0):
    dt = h*(n)-t_o
    return [i*dt for i in range(n+1)]
def rkuta4(f,f0,n,h):
    t = 0
    b = h*(n+1)
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
    p, q, T = np.zeros(n), np.zeros(n), []


    return p, q, T

def StormerV(f,g,f0,n,h):
    p, q, T = [], [], []

    return p, q, T

def dynamics_solve(f, f0 = 1, h = 0.1, n = 100, D = 1, t0 = 0, method = "Euler"):
    method_dict = {
        "Euler": euler,
        "RK2": rkuta2,
        "RK4": rkuta4
    }
    if method not in method_dict:
        return None
    
    selected_method = method_dict[method]
    vsol = T(n,h,t0)
    vsol.extend(selected_method(f[i], f0, n, h)[0] for i in range(D))
    return np.array(vsol)
def hamiltonian_solve(f,g, f0 = 1, h = 0.1, n = 100, D = 1, t0 = 0, method = "Euler"):
    if method != "SE" or method != "SV":
        for i in range(D):
            vsol = [dynamics_solve(f[i],)]
    elif method == "SE":
        vsol = [SEuler(f[0], g[0], f0, n, h)[1]]
        for i in range(D):
            vsol.append(SEuler(f[i], g[i], f0, n, h)[0])
        return np.array(vsol)
    elif method == "SV":
        vsol = [StormerV(f[0], g[0], f0, n, h)[1]]
        for i in range(D):
            vsol.append(StormerV(f[i], g[i], f0, n, h)[0])
        return np.array(vsol)
    else:
        return "None"


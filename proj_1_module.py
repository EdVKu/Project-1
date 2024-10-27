import numpy as np
def euler(f, f0, n, h, t_o = 0):
    t = t_o
    b = h*(n)-t_o
    s0 = f0
    s, T = [],[]
    while t <= b:
        ssq = s0 + h*f(t, s0)
        s.append(ssq)
        T.append(t)
        s0 = ssq
        t += h
    return s, T

def rkuta2(f, f0, n, h, t_o = 0):
    t = t_o
    b = h*(n)-t_o
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

    
def rkuta4(f, f0, n, h, t_o = 0):
    t = 0
    b = h*(n)-t_o
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
def SEuler(f,g,f0,g0,h,n,t_o):
    p, q= np.zeros(n+1), np.zeros(n+1)
    T = [t_o]
    p[0], q[0] = f0, g0
    t = t_o
    for i in range(n):
        p[i+1] = p[i] + f(q[i]) * h
        q[i+1] = q[i] + g(p[i+1]) * h
        t += h
        T.append(t)

    return p, q, T

def StormerV(f,g,f0,g0,h,n,t_o):
    p, q, T = [], [], []

    return p, q, T


def dynamics_solve(f, f0 = 1, h = 0.1, n = 100, D = 1, t_o = 0, method = "Euler"):
    method_dict = {
        "Euler": euler,
        "RK2": rkuta2,
        "RK4": rkuta4
    }
    if method not in method_dict:
        return None
    
    selected_method = method_dict[method]
    vsol = [selected_method(f[0], f0, n, h, t_o)[-1]]
    for i in range(D):
        vsol.append(selected_method(f[i], f0, n, h, t_o)[0])

    return np.array(vsol)
def hamiltonian_solve(f,g, f0 = 1,g0=1, h = 0.1, n = 100, D = 1, t_o = 0, method = "Euler"):
    if method == "Euler":
      vsolp, vsolq = [euler(f[0],f0,h,n,t_o)[-1]],[euler(g[0],g0,h,n,t_o)[-1]]
      for i in range(D):
        vsolp.append(euler(f[i],f0,h,n,t_o))
        vsolq.append(euler(g[i],g0,h,n,t_o))

      return np.array(vsolp),np.array(vsolq)
    elif method == "RK2":
      vsolp, vsolq = [rkuta2(f[0],f0,h,n,t_o)[-1]],[rkuta2(g[0],g0,h,n,t_o)[-1]]
      for i in range(D):
        vsolp.append(rkuta2(f[i],f0,h,n,t_o))
        vsolq.append(rkuta2(g[i],g0,h,n,t_o))

      return np.array(vsolp),np.array(vsolq)
    elif method == "RK4":
      vsolp, vsolq = [rkuta4(f[0],f0,h,n,t_o)[-1]],[rkuta4(g[0],g0,h,n,t_o)[-1]]
      for i in range(D):
        vsolp.append(rkuta2(f[i],f0,h,n,t_o))
        vsolq.append(rkuta2(g[i],g0,h,n,t_o))

      return np.array(vsolp),np.array(vsolq)
    elif method == "SV":
      stor = StormerV(f[0],g[0],f0,g0,h,n,t_o)
      solp, solq = [],[]
      for i in range(D):
        solp.append(StormerV(f[i],g[i],f0,g0,h,n,t_o)[0])
        solq.append(StormerV(f[i],g[i],f0,g0,h,n,t_o)[1])
      
      return np.array(solp), np.array(solq), stor[-1]
    elif method == "SE":
      seu = SEuler(f[0],g[0],f0,g0,h,n,t_o)
      solp, solq = [],[]
      for i in range(D):
        psol,qsol = SEuler(f[i],g[i],f0,g0,h,n,t_o)
        solp.append(psol)
        solq.append(qsol)
      
      return np.array(solp), np.array(solq), seu[-1]
    else:
      return "None"
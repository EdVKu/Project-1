import numpy as np
def euler(f, f0, n, hf, t_o = 0):
    t = t_o
    h = ((hf)-t_o)/n
    s0 = f0
    s, T = [],[]
    while t <= hf:
        ssq = s0 + h*f(t, s0)
        s.append(ssq)
        T.append(t)
        s0 = ssq
        t += h
    return s, T
def euler2(f, f0,g0, n, hf, t_o = 0):
  h = ((hf)-t_o)/n
  s, v, T = [f0], [g0], [t_o]
  for i in range(n):
    snext = s[i] + h*v[i]
    vnext = v[i] + h*f(T[i],s[i])
    s.append(snext)
    v.append(vnext)
    T.append(T[i]+h)
  return v,s, T

def rkuta2(f, f0, n, hf, t_o = 0):
    t = t_o
    h = ((hf)-t_o)/n
    s0 = f0
    s,T = [],[]
    while t<=hf:
        k1n = h*f(t, s0)
        k2n = h*f(t, s0 + 0.5*k1n)
        s.append(s0 + k2n)
        T.append(t)
        s0 += k2n
        t += h
    return s, T

def rk2_2(f, f0,g0, n, hf, t_o = 0):
  s = [f0]
  v = [g0]
  h = ((hf)-t_o)/n
  T = [t_o]
  for i in range(n):
    k1s = h*s[i]
    k1v = h*f(T[i],s[i])
    k2s = h * (v[i] + 0.5 * k1v)
    k2v = h * f(T[i] + 0.5 * h, s[i] + 0.5 * k1s)
    sp = s[i] + k2s
    vp = v[i] + k2v
    s.append(sp)
    v.append(vp)
    T.append(T[i] + h)
  return v,s, T

def rkuta4(f, f0, n, hf, t_o = 0):
    t = 0
    h = ((hf)-t_o)/n
    s0 = f0
    s,T = [],[]

    while t<=hf:
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
def rk4_2(f, f0,g0, n, hf, t_o = 0):
  s = [f0]
  v = [g0]
  T = [t_o]
  h = ((hf)-t_o)/n
  for i in range(n):
    k1s = h*s[i]
    k1v = h*f(T[i],s[i])
    k2s = h * (v[i] + 0.5 * k1v)
    k2v = h * f(T[i] + 0.5 * h, s[i] + 0.5 * k1s)
    k3s = h * (v[i] + 0.5 * k2v)
    k3v = h * f(T[i] + 0.5 * h, s[i] + 0.5 * k2s)
    k4s = h * (v[i] + k3v)
    k4v = h * f(T[i] + 0.5 * h, s[i] + k3s)
    sp = s[i] + (k1s + 2 * k2s + 2 * k3s + k4s) / 6
    vp = v[i] + (k1v + 2 * k2v + 2 * k3v + k4v) / 6
    s.append(sp)
    v.append(vp)
    T.append(T[i] + h)
  return  v,s, T


def StormerV(f, f0,g0, n, h, t_o = 0):
  p, q = np.zeros(n), np.zeros(n)
  T = [t_o, t_o + h]
  p[0], q[0] = f0, g0
  p[1] = p[0] + q[0]*h + h**2*f(T[0],p[0])*0.5
  for i in range(1,n-1):
    p[i+1] = p[i]+q[i]*h +f(T[i],p[i])*h**2
    q[i+1] = q[i] + f(T[i],p[i])*h
    T.append(T[i] + h)

  return list(p), list(q), T

def SEuler(f,g, f0,g0, n, h, t_o = 0):
  p, q= np.zeros(n+1), np.zeros(n+1)
  T = [t_o]
  p[0], q[0] = f0, g0
  for i in range(n):
    p[i+1] = p[i] + f(T[i],q[i]) * h
    q[i+1] = q[i] + g(T[i],p[i+1]) * h
    T.append(T[i] + h)

  return list(p), list(q), T

def dynamics_solve(f, f0 = [1],  n = 100,hf = 100, D = 1, t_o = 0, method = "Euler"):
    method_dict = {
        "Euler": euler,
        "RK2": rkuta2,
        "RK4": rkuta4
    }
    if method not in method_dict:
        return None

    selected_method = method_dict[method]
    vsol = [selected_method(f[0], f0[0], n, hf, t_o)[-1]]
    for i in range(D):
        vsol.append(selected_method(f[i], f0[i], n, hf, t_o)[0])

    return vsol
def hamiltonian_solve(f,g, f0 = [1],g0=[1],  n = 100, h = 0.1,hf = 100, D = 1, t_o = 0, method = "Euler"):
  if method == "Euler":
    if D==1:
      vsolp, vsolq, T = euler2(f[0], f0[0], g0[0], n, hf, t_o)
      return (vsolp), (vsolq), (T)
    else:
      vsolp, vsolq, T = euler2(f[0], f0[0], g0[0], n, hf, t_o)
      for i in range(D):
        eul = euler2(f[i], f0[i], g0[i], n, hf, t_o)
        vsolp.append(eul[0])
        vsolq.append(eul[1])

      return vsolp, vsolq, T
  elif method == "RK2":
    if D==1:
      vsolp, vsolq, T = rk2_2(f[0], f0[0], g0[0], n, hf, t_o)
      return vsolp, vsolq, T
    else:
      vsolp, vsolq, T = rk2_2(f[0], f0[0], g0[0], n, hf, t_o)
      for i in range(D):
        eul = rk2_2(f[i], f0[i], g0[i], n, hf, t_o)
        vsolp.append(eul[0])
        vsolq.append(eul[1])

      return vsolp, vsolq, T
  elif method == "RK4":
    if D==1:
      vsolp, vsolq, T = rk4_2(f[0], f0[0], g0[0], n, hf, t_o)
      return vsolp, vsolq, T
    else:
      vsolp, vsolq, T = rk4_2(f[0], f0[0], g0[0], n, hf, t_o)
      for i in range(D):
        eul = rk4_2(f[i], f0[i], g0[i], n, hf, t_o)
        vsolp.append(eul[0])
        vsolq.append(eul[1])

      return vsolp, vsolq, T
  elif method == "SV":
    solp, solq, T = StormerV(f[0], f0[0], g0[0], n, h, t_o)
    if D == 1:
      return solp, solq, T
    else:
      for i in range(D):
        sverle = StormerV(f[i], f0[i], g0[i], n, h, t_o)
        solp.append(sverle[0])
        solq.append(sverle[1])
      return solp, solq, T
  elif method == "SE":
    solp, solq, T = SEuler(f[0], g[0], f0[0], g0[0], n, h, t_o)
    if D ==1:
      return solp, solq, T
    else:
      for i in range(D):
        seu = SEuler(f[i], g[i], f0[i], g0[i], n, h, t_o)
        solp.append(seu[0])
        solq.append(seu[1])

    return solp, solq, T
  else:
    return "None"
import numpy as np

def dynamics_solve(f, f0, h, n):
    method = input("please enter a number to choose your method (0 = Euler, 1 = Runge Kutta 2, 2 = Runge Kutta 4)")
    if int(method) == 0:
        t = 0
        b = n*h
        s0 = f0
        s = []
        while t <= b:
            ssq = s0 + h*f(t,s0)
            s.append(ssq)
            s0 = ssq
            t += h
        return s
    elif int(method) == 1:
        return "World"
    elif int(method) == 2:
        return "!:D"
    else:
        return None
def hamiltonian_solve():
    return 0

dynamics_solve(lambda x : np.exp(-x), 1,5e-5,100)


import numpy as np
from scipy.special import comb, factorial

dt = 0.001
t = 20
eps = 0.0
delta = 0.5
c = np.sqrt(2) / 1
beta = 16.0
omega = 1.0 # - 0.1j  # the imaginary part characterizes the loss. (Brownian oscillator model)

initstate = 0 # hard-coded here... should be flexible
tier = 51

def RK4(F, y0, h):

    k1 = h * F
    k2 = h * F + h * k1 / 2
    k3 = h * F + h * k2 / 2
    k4 = h * F + h * k3
    return y0 + (k1 + 2 * k2 + 2 * k3 + k4) / 6

def init(M = tier, N = 1):

    omega_ = np.real(omega)
    factor = 1 - np.exp( - beta * omega_)
    A = np.zeros((M, 4 * (M + 1) * N**M), dtype=complex)
    A[0, 0] = 0.0
    A[0, 1] = 0.0
    A[0, 2] = 1.0
    A[0, 3] = 1.0
    for j in range(1, int(M/2)):
        for n in range(100):
            A[2 * j, 5 * j + 2] += factor * np.exp(- j * beta * omega_) * factorial(j) * comb(n,j) * np.exp(- n * beta * omega_)
        A[2 * j, 7 * j + 3] = A[2 * j, 5 * j + 2]

    return(A)

def prop(A):

    M = tier
    # tier = 0
    A[0, 0] = RK4(- 2 * eps * A[0, 1] - np.sqrt(2) * c * (A[1, 2] + A[1, 3]), A[0, 0], dt)
    A[0, 1] = RK4(2 * eps * A[0, 0] - 2 * delta * A[0, 2] + np.sqrt(2) * c * (A[1, 0] + A[1, 1]), A[0, 1], dt)
    A[0, 2] = RK4(2 * delta * A[0, 1], A[0, 2], dt)
    A[0, 3] = 1.0

    # tier = 1, 2, ..., M - 1
    for j in range(1, M - 1):

        for k in range(j + 1):

            # omega -> Re[omega] + (-1)^(k) i Im[omega]. EOM needs to be proved... imaginary part is correct. How to prove the real part? why the absolute value?

            H_x = (- 2 * eps * A[j, j + k + 1] + ( np.abs(j - 2 * k) * np.imag(omega) - 1.0j * (j - 2 * k) * np.real(omega) ) * A[j, k] 
                - np.sqrt(2) * c * (A[j + 1, j + k + 2] + A[j + 1, j + k + 3] + 0.5 * (k * A[j - 1, j + k - 1] + (j - k) * A[j - 1, j + k])))
            H_y = (2 * eps * A[j, k] - 2 * delta * A[j, 2 * j + 2 + k] + ( np.abs(j - 2 * k) * np.imag(omega) - 1.0j * (j - 2 * k) * np.real(omega) ) * A[j, j + 1 + k] 
                + np.sqrt(2) * c * (A[j + 1, k] + A[j + 1, k + 1] + 0.5 * (k * A[j - 1, k - 1] + (j - k) * A[j - 1, k])))
            H_z = (2 * delta * A[j, j + 1 + k] + ( np.abs(j - 2 * k) * np.imag(omega) - 1.0j * (j - 2 * k) * np.real(omega) ) * A[j, 2 * j + 2 + k] 
                + 1.0j * np.sqrt(2) * c * 0.5 * (k * A[j - 1, 3 * j + k - 1] - (j - k) * A[j - 1, 3 * j + k]))
            I_H = (( np.abs(j - 2 * k) * np.imag(omega) - 1.0j * (j - 2 * k) * np.real(omega) ) * A[j, 3 * j + 3 + k] 
                + 1.0j * np.sqrt(2) * c * 0.5 * (k * A[j - 1, 2 * j + k - 1] - (j - k) * A[j - 1, 2 * j + k]))

            A[j, k] = RK4(H_x, A[j, k], dt)
            A[j, j + 1 + k] = RK4(H_y, A[j, j + 1 + k], dt)
            A[j, 2 * j + 2 + k] = RK4(H_z, A[j, 2 * j + 2 + k], dt)
            A[j, 3 * j + 3 + k] = RK4(I_H, A[j, 3 * j + 3 + k], dt)

    return A

def dynamics(M, N): # here N = 1 is hard-coded here, need to generate to N > 1 cases.

    A = init(M, N)
    PiiFile = open("Pt.txt","w") 
    for k in range(int(t / dt)):
        A = prop(A)
        if (k % int(0.1 / dt) == 0):
            PiiFile.write(f"{round(k * dt, 3)} \t")
            PiiFile.write(str(np.real(A[0, 2])) + "\t")
            PiiFile.write("\n")
    PiiFile.close()

if __name__ == '__main__':

    dynamics(M = tier, N = 1)

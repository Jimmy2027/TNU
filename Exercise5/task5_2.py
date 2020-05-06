import numpy as np
from scipy.special import gammaln

mu = tau = 1
N = 10
eps = np.random.normal(mu, 1, N)
y = mu + eps

mu0 = 0
lambda0 = 3
a0 = 2
b0 = 2


def update_mu(a, b):
    s2 = 1 / ((N + lambda0) * a / b)
    m = (lambda0 * mu0 + N * np.mean(y)) / (lambda0 + N)
    return s2, m


def update_tau(a, b, m, s2):
    a = a + (N + 1) / 2
    b = b + 0.5 * np.sum((y - m) ** 2 + lambda0 * (mu0 - m) ** 2 + (N + lambda0) * s2)
    return a, b


def evaluate_free_energy(a, b, s2):
    return -a * np.log(b) + gammaln(a) - gammaln(a0) + a0 * np.log(b0) + 0.5 * np.log(lambda0) + np.log(
        np.sqrt(s2)) - 0.5 * N * np.log(2 * np.pi) + 0.5


"""
Starting values
"""
diff = free_energy = 1000
a, b = a0, b0
mu, s2 = mu0, b0 / (a0 * lambda0)

while diff > 0.001:
    s2_new, mu_new = update_mu(a, b)
    a_new, b_new = update_tau(a, b, mu, s2)
    new_free_energy = evaluate_free_energy(a, b, s2)
    diff = np.abs(free_energy - new_free_energy)
    print(free_energy, new_free_energy, diff)
    # print(a,b,m,s2)
    free_energy = new_free_energy
    a, b = a_new, b_new
    mu, s2 = mu_new, s2_new
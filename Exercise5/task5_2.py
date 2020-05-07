import numpy as np
from scipy.special import gammaln
import pandas as pd

df = pd.DataFrame()
mu = tau = 1
N = 10
eps = np.random.normal(mu, 1, N)
y = mu + eps

mu0 = 0
lambda0 = 3
a0 = 2
b0 = 2
s20 = b0 / (a0 * lambda0)


def update_mu(a, b):
    s2 = 1 / ((N + lambda0) * a / b)
    m = (lambda0 * mu0 + N * np.mean(y)) / (lambda0 + N)
    return s2, m


def update_tau(m, s2):
    a = a0 + (N + 1) / 2
    b = b0 + 0.5 * np.dot(np.transpose(y - np.ones(y.shape)), (y - np.ones(y.shape))) + lambda0 * (mu0 - m) ** 2 + (
            N + lambda0) * s2
    return a, b


def evaluate_free_energy(a, b, s2):
    return -a * np.log(b) + gammaln(a) - gammaln(a0) + a0 * np.log(b0) + 0.5 * np.log(lambda0) + np.log(
        np.sqrt(s2)) - 0.5 * N * np.log(2 * np.pi) + 0.5


"""
Starting values
"""
for i in range(10):
    diff = free_energy = 1000
    if not i == 0:
        a0 = np.random.randint(1, 10)
        b0 = np.random.randint(1, 10)
        mu0 = np.random.randint(0, 10)
        s20 = np.random.randint(0, 10)
    a, b = a0, b0
    mu, s2 = mu0, s20
    it = 0
    while diff > 0.0001:
        s2_new, mu_new = update_mu(a, b)
        a_new, b_new = update_tau(mu, s2)
        new_free_energy = evaluate_free_energy(a, b, s2)
        diff = np.abs(free_energy - new_free_energy)
        print(free_energy, new_free_energy, diff)
        print(a_new, b_new, mu_new, s2_new)
        free_energy = new_free_energy
        a, b = a_new, b_new
        mu, s2 = mu_new, s2_new
        it += 1

    df = df.append(
        pd.DataFrame([[a0, b0, mu0, s20, it, free_energy]], columns=['a0', 'b0', 'mu0', 's20', 'it', 'free_energy']))
df.to_csv('results.csv', index = False)
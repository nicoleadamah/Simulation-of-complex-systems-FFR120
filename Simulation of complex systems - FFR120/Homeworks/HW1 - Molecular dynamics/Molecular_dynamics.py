import numpy as np
from matplotlib import pyplot as plt
from tqdm import trange
from typing import Callable, Tuple
import argparse

SNAPSHOT = 0
ENERGIES = 1

def leapfrog(x0: np.ndarray, v0: np.ndarray, xh_f: Callable[[np.ndarray], np.ndarray], dt = 1.0) -> Tuple[np.ndarray, np.ndarray]:
    x0 = x0 + v0 * dt/2
    xh = xh_f(x0)
    v = v0 + xh * dt
    x = x0 + v * dt/2
    return (x, v)

def initialize_pos(size: Tuple[int, int], L, sigma):
    pos = np.zeros(size)

    for i in range(size[0]):
        point_placed = False

        while not point_placed:
            new_pos = np.random.rand(1, 2) * L
            distances = [ np.sqrt(np.sum((new_pos - pos[j,:])**2)) for j in range(i) ]
            instabilities = [ r > sigma for r in distances ]
            if all(instabilities):
                pos[i,:] = new_pos
                point_placed = True
    return pos

def initialize_v(size: Tuple[int, int], v0):
    angles = (np.random.rand(size[0]) * 2 * np.pi).reshape(size[0],1)
    directions = np.column_stack((np.cos(angles), np.sin(angles)))
    return directions * v0

def lennard_jones_f(positions, e, s):
    f = np.zeros_like(positions)
    for i in range(positions.shape[0]):
        for j in range(i+1, positions.shape[0]):
            r = np.sqrt(np.sum((positions[i,:] - positions[j,:])**2))
            magnitude = 4 * e * (12 * np.power(s, 12) * np.power(r, -13) - 6 * np.power(s, 6) * np.power(r, -7))
            direction = (positions[j,:] - positions[i,:]) / r
            f[i,:] -= magnitude*direction
            f[j,:] += magnitude*direction
    return f

def lennard_jones_p(positions:np.ndarray, e, s):
    p = np.zeros(positions.shape[0])
    for i in range(positions.shape[0]):
        for j in range(i+1, positions.shape[0]):
            r = (np.sum((positions[i,:] - positions[j,:])**2))
            magnitude = 4 * e * (np.power(s ** 2 / r, 6) - np.power(s ** 2 / r, 3))
            p[i] += magnitude
            p[j] += magnitude
    return p

def potential_energy(x, e, s):
    return 0.5 * np.sum(lennard_jones_p(x, e, s))

def kinetic_energy(v, m, v_scale):
    return 0.5 * m * np.sum((v/v_scale)**2)

def plot_snap(N, position_history):
    for i in range(N):
        plt.plot(position_history[:, i, 0].squeeze(), position_history[:, i, 1].squeeze())
    plt.gca().set_xlim([0, L])
    plt.gca().set_ylim([0, L])
    plt.ylabel('$y$')
    plt.xlabel('$x$')
    plt.title('Snapchot of Lennard-Jones gas')
    plt.show()

def plot_energies(e,steps,freq,t_scale,Ek,Ep,dt):
    plt.subplot(3, 1, 1)
    plt.plot(np.arange(steps // freq) * dt / t_scale, Ek / e, '', label="Kinetic Energy",markersize=1)
    plt.ylabel("$E_k$")
    plt.subplot(3, 1, 3)
    plt.plot(np.arange(steps // freq) * dt / t_scale, Ep / e, '', label="Total Energy", markersize=1)
    plt.ylabel("$E_p$")
    plt.subplot(3, 1, 2)
    plt.plot(np.arange(steps // freq) * dt / t_scale, (Ek + Ep / e), '',label="Potential Energy", markersize=1)
    plt.ylabel("$E$")
    plt.show()

def main(args):
    e = 1.0 #epsilon
    s = 1.0 #sigma
    m = 1 #mass
    v_scale = np.sqrt(2 * e / m)
    t_scale = s * np.sqrt(m / (2 * e))
    L = s * 100
    N = 10
    steps = 40000
    freq = 5
    dt = 0.001
    #dt = s / (2 * v_scale) * 0.02
    x = initialize_pos((N, 2), L, s) #positions
    v = initialize_v((N, 2), 2 * v_scale) #velocities
    x_history = np.empty((steps // freq, N, 2))
    Ek = np.empty(steps // freq)
    Ep = np.empty(steps // freq)

    plotting = ENERGIES
    for t in trange(steps):
        if plotting == SNAPSHOT:
            x_history[t // freq, :, :] = x
        if t % freq == 0 and plotting == ENERGIES:

            Ek[t // freq] = kinetic_energy(v, m, v_scale)
            Ep[t // freq] = potential_energy(x, e, s)

        (x, v) = leapfrog(x, v, lambda x: lennard_jones_f(x, e, s), dt=dt)

        for i in range(N):
            # Outside left bound
            if x[i, 0] < 0:
                x[i, 0] = x[i, 0] * -1
                v[i, 0] = v[i, 0] * -1

            # Outside right bound
            elif x[i, 0] > L:
                x[i, 0] = 2 * L - x[i, 0]
                v[i, 0] = v[i, 0] * -1

            # Outside lower bound
            elif x[i, 1] < 0:
                x[i, 1] = x[i, 1] * -1
                v[i, 1] = v[i, 1] * -1

            # Outside upper bound
            elif x[i, 1] > L:
                x[i, 1] = 2 * L - x[i, 1]
                v[i, 1] = v[i, 1] * -1

    if plotting == SNAPSHOT:
        plot_snap(N, x_history)
    elif plotting == ENERGIES:
        plot_energies(e, steps, freq, t_scale, Ek, Ep, dt)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulate molecular dynamics")
    parser.add_argument("--outdir", "-o", type=str, default=".", help="Out directory")
    args = main(parser.parse_args())
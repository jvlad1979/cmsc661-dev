import numpy as np
from .rk import RungeKuttaSolver
import pytest
from numpy.testing import assert_allclose

# Define RungeKuttaSolver class (assuming it's already imported)

# Test 1: Exponential growth

def test_exponential_growth_erk4():
    k = 1.0
    func = lambda t, y: k * y
    y0 = np.array([1.0])
    t0 = 0.0
    h = 0.1
    N = 100

    # Explicit RK4 coefficients
    A = np.array([[0, 0, 0, 0],
                  [0.5, 0, 0, 0],
                  [0, 0.5, 0, 0],
                  [0, 0, 1, 0]])
    b = np.array([1/6, 1/3, 1/3, 1/6])
    c = np.array([0, 0.5, 0.5, 1])

    solver = RungeKuttaSolver(A, b, c)
    t, y = solver.solve(func, t0, y0, h, N)

    y_exact = y0 * np.exp(k * t)
    assert_allclose(y.flatten(), y_exact, rtol=1e-4)

# Test 2: Linear ODE

def test_linear_ode_erk4():
    func = lambda t, y: np.array([t])
    y0 = np.array([0.0])
    t0 = 0.0
    h = 0.1
    N = 100

    # Explicit RK4 coefficients
    A = np.array([[0, 0, 0, 0],
                  [0.5, 0, 0, 0],
                  [0, 0.5, 0, 0],
                  [0, 0, 1, 0]])
    b = np.array([1/6, 1/3, 1/3, 1/6])
    c = np.array([0, 0.5, 0.5, 1])

    solver = RungeKuttaSolver(A, b, c)
    t, y = solver.solve(func, t0, y0, h, N)

    y_exact = 0.5 * t**2
    assert_allclose(y.flatten(), y_exact, rtol=1e-4)

# Test 3: Simple Harmonic Oscillator

def test_simple_harmonic_oscillator_erk4():
    omega = 1.0

    def func(t, y):
        return np.array([y[1], -omega**2 * y[0]])

    y0 = np.array([1.0, 0.0])  # Initial position 1, velocity 0
    t0 = 0.0
    h = 0.1
    N = 100

    # Explicit RK4 coefficients
    A = np.array([[0, 0, 0, 0],
                  [0.5, 0, 0, 0],
                  [0, 0.5, 0, 0],
                  [0, 0, 1, 0]])
    b = np.array([1/6, 1/3, 1/3, 1/6])
    c = np.array([0, 0.5, 0.5, 1])

    solver = RungeKuttaSolver(A, b, c)
    t, y = solver.solve(func, t0, y0, h, N)

    y_exact = np.cos(omega * t)
    assert_allclose(y[:, 0], y_exact, rtol=1e-3)

def test_different_implicit_methods():
    k = 1.0
    func = lambda t, y: k * y
    jac = lambda t, y: np.array([[k]])
    y0 = np.array([1.0])
    t0 = 0.0
    h = 0.001
    N = 100

    methods = {
        "Implicit Euler": (np.array([[1]]), np.array([1]), np.array([1])),
        "Trapezoidal Rule": (np.array([[0.5, 0.5], [0.5, 0.5]]), np.array([0.5, 0.5]), np.array([0.5, 0.5])),
        "Radau IIA (Order 3)": (np.array([[5/12, -1/12], [3/4, 1/4]]), np.array([3/4, 1/4]), np.array([1/3, 1]))
    }

    for method_name, (A, b, c) in methods.items():
        solver = RungeKuttaSolver(A, b, c)
        t, y = solver.solve(func, t0, y0, h, N, jac)
        y_exact = y0 * np.exp(k * t)
        assert_allclose(y.flatten(), y_exact, rtol=1e-4, err_msg=f"Failed for {method_name}")
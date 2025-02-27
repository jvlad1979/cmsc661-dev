import numpy as np
from scipy.linalg import solve

class BDFSolver:
    # BDF coefficients for N = 1 to 6
    BDF_COEFFICIENTS = {
        1: (np.array([1, -1]), 1),
        2: (np.array([3/2, -2, 1/2]), 2/3),
        3: (np.array([11/6, -3, 3/2, -1/3]), 6/11),
        4: (np.array([25/12, -4, 3, -4/3, 1/4]), 12/25),
        5: (np.array([137/60, -5, 5, -10/3, 5/4, -1/5]), 60/137),
        6: (np.array([147/60, -6, 15/2, -20/3, 15/4, -6/5, 1/6]), 60/147)
    }

    def __init__(self, order, max_iter=100, tol=1e-8, initial_solver=None):
        if order < 1 or order > 6:
            raise ValueError("BDF order must be between 1 and 6.")
        self.order = order
        self.alpha, self.gamma = self.BDF_COEFFICIENTS[order]
        self.max_iter=max_iter
        self.tol=tol
        if initial_solver is None:
            initial_solver = BDFSolver(order - 1)
        else:
            self.initial_solver = initial_solver

    def solve(self, func, jac, t0, y0, h, N):
        t = np.zeros(N + 1)
        y = np.zeros((N + 1, len(y0)))
        t[:self.order] = t0 + h * np.arange(self.order)

        # Use initial solver for the first few steps
        t_init, y_init = self.initial_solver.solve(func, t0, y0, h, self.order - 1, jac)
        y[:self.order] = y_init

        for i in range(self.order - 1, N):
            y_guess = y[i]
            for iteration in range(self.max_iter):  # Newton-Raphson iterations
                F = self._bdf_residual(func, t[i+1], y_guess, y[i - self.order + 1:i+1], h)
                J = self._bdf_jacobian(jac, t[i+1], y_guess, h)
                delta_y = solve(J, -F)
                y_guess += delta_y

                if np.linalg.norm(delta_y) < self.tol:
                    break
            else:
                raise RuntimeError("Newton-Raphson did not converge in BDF method.")

            y[i+1] = y_guess
            t[i+1] = t[i] + h

        return t, y

    def _bdf_residual(self, func, t_new, y_new, y_prev, h):
        return self.alpha[0] * y_new - np.sum(self.alpha[1:, np.newaxis] * y_prev[::-1], axis=0) - h * self.gamma * func(t_new, y_new)

    def _bdf_jacobian(self, jac, t_new, y_new, h):
        n = len(y_new)
        return self.alpha[0] * np.eye(n) - h * self.gamma * jac(t_new, y_new)
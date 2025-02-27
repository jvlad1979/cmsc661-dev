import numpy as np
from scipy.linalg import solve

class RungeKuttaSolver:
    def __init__(self,A,b,c,max_iter=100,tol=1e-6):
        self.A = A
        self.b = b
        self.c = c
        self.max_iter = max_iter
        self.tol = tol
    
    def solve_explicit_step(self,func,t_i,y_i,h,):
        s = len(self.b)
        n = len(y_i)
        k = np.zeros((s, n))

        for j in range(s):
            t_stage = t_i + self.c[j] * h
            y_stage = y_i + h * np.sum(self.A[j, :j][:, np.newaxis] * k[:j], axis=0)
            k[j] = func(t_stage, y_stage)

        return k

    # def solve_implicit_step(self, func, jac, t_i, y_i, h):
    #     """Solves implicit step using Newton-Raphson with precomputations and np.linalg.solve()."""
    #     s = len(self.b)  # Number of stages
    #     n = len(y_i)  # Dimension of the system
    #     k = np.zeros((s, n))

    #     # Precompute h * A to avoid redundant multiplications

    #     for j in range(s):
    #         t_stage = t_i + self.c[j] * h
    #         k_j = np.zeros(n)  # Initial guess for k_j

    #         # Precompute Jacobian matrix once per Newton-Raphson iteration
    #         J_func = jac(t_stage, y_i)  # Evaluate Jacobian at base point
    #         J = np.eye(n) - h * self.A[j, j] * J_func  # Newton matrix

    #         for _ in range(self.max_iter):
    #             # Compute y_stage using precomputed `hA`
    #             y_stage = y_i + np.dot(h * self.A[j, :], k)

    #             F = k_j - func(t_stage, y_stage)  # Residual function
    #             delta_k = np.linalg.solve(J, -F)  # Solve linear system instead of inverting

    #             k_j += delta_k  # Update k_j

    #             if np.linalg.norm(delta_k) < self.tol:
    #                 break  # Converged

    #         k[j] = k_j  # Store computed k_j

    #     return k

    def solve_implicit_step(self, func, jac, t_i, y_i, h):
        """Optimized Newton-Raphson solver with adaptive Jacobian updates and better initial guesses."""
        s = len(self.b)  # Number of stages
        n = len(y_i)  # Dimension of the system
        k = np.zeros((s, n))
        # Precompute h * A
        k[0] = func(t_i,y_i)
        for j in range(s):
            t_stage = t_i + self.c[j] * h
            # k_j = k[j]
            for iter in range(self.max_iter):
                y_stage = y_i + np.dot(h * self.A[j, :j], k[:j])
                F = k[j] - func(t_stage, y_stage)  
                # Update Jacobian adaptively (mimicking DIRK2)
                DF = np.eye(n) - h * self.A[j, j] * jac(t_stage, y_stage)  
                delta_k = np.linalg.solve(DF, F)  
                k[j] -= delta_k  # Update k_j
                if np.linalg.norm(delta_k) < self.tol:
                    break  # Converged
            print('monke')
            # if iter%10 == 0:
            #     print(iter)
            #     print(k)
            # k[j] = k_j  # Store computed k_j

        return k


    def solve_k(self,func,jac,t_i,y_i,h):
        if np.allclose(self.A, np.tril(self.A)) and np.all(np.diag(self.A) == 0):
            k = self.solve_explicit_step(func,t_i,y_i,h)
        else:
            if jac is None:
                raise NotImplementedError("Jacobian is required for implicit methods... For now...")
            k = self.solve_implicit_step(func,jac,t_i,y_i,h)
        return k

    def solve(self,func,t0,y0,h,N,jac=None):
        t = np.zeros(N+1)
        y = np.zeros((N+1,len(y0)))
        t[0] = t0
        y[0] = y0
        for i in range(N):
            k = self.solve_k(func,jac,t[i],y[i],h)
            y[i + 1] = y[i] + h * np.sum(self.b[:, np.newaxis] * k, axis=0)
            t[i + 1] = t[i] + h
        return t, y

class DIRK2(RungeKuttaSolver):
    def __init__(self,max_iter=100,tol=1e-8):
        gamma = 1-1/np.sqrt(2)
        A = np.array([[gamma,0],[1-gamma,gamma]])
        b = np.array([1-gamma,gamma])
        c = np.array([gamma,1])
        super().__init__(A,b,c,max_iter=max_iter,tol=tol)

class DIRK3(RungeKuttaSolver):
    def __init__(self,max_iter=100,tol=1e-8):
        gamma = 1/2 + np.sqrt(3)/6
        A = np.array([[gamma,0],[1-2*gamma,gamma]])
        b = np.array([1/2,1/2])
        c = np.array([gamma,1-gamma])
        super().__init__(A,b,c,max_iter=max_iter,tol=tol)

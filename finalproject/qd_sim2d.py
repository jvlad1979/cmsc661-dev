import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
from scipy.sparse import diags, kron, identity, csc_matrix
from scipy.sparse.linalg import eigsh, spsolve

# Constants
q = const.e
eps0 = const.epsilon_0
hbar = const.hbar
m0 = const.m_e
kB = const.k

def create_laplacian_2d(Nx, Ny, dx, dy):
    """Constructs 2D Laplacian with Dirichlet BCs using Kronecker products."""
    main_x = -2 * np.ones(Nx)
    off_x = np.ones(Nx - 1)
    Dx = diags([off_x, main_x, off_x], [-1, 0, 1], shape=(Nx, Nx)) / dx**2

    main_y = -2 * np.ones(Ny)
    off_y = np.ones(Ny - 1)
    Dy = diags([off_y, main_y, off_y], [-1, 0, 1], shape=(Ny, Ny)) / dy**2

    L = kron(identity(Ny), Dx) + kron(Dy, identity(Nx))
    return csc_matrix(L)

def triple_dot_potential(X, Y, V0=0.2, w=30e-9, sep=60e-9):
    """Creates a triple quantum dot confinement potential (sum of Gaussians)."""
    centers = [-sep, 0.0, sep]
    V = np.zeros_like(X)
    for xc in centers:
        V += V0 * np.exp(-((X - xc)**2 + Y**2) / (2 * w**2))
    return V

def schrodinger_solver_2d(V, dx, dy, m_eff, num_states=3):
    """Solves 2D Schrödinger equation using finite differences."""
    Nx, Ny = V.shape
    N = Nx * Ny
    V_flat = V.ravel()
    L = create_laplacian_2d(Nx, Ny, dx, dy)
    prefactor = - (hbar**2) / (2 * m_eff * m0 * q)
    H = prefactor * L + diags(V_flat, 0)
    E, psi = eigsh(H, k=num_states, which='SA')
    psi = psi.reshape((Nx, Ny, num_states))
    return E, psi

def poisson_solver_2d(rho, dx, dy, eps_r):
    """Solves 2D Poisson equation with Dirichlet BCs (phi=0 at edges)."""
    Nx, Ny = rho.shape
    L = create_laplacian_2d(Nx, Ny, dx, dy)
    rhs = -rho.ravel() / (eps_r * eps0)
    phi = spsolve(L, rhs)
    return phi.reshape((Nx, Ny))

def fermi_density(psi, E, Ef, T, dx, dy):
    """Computes 2D electron density from wavefunctions and eigenvalues."""
    rho = np.zeros_like(psi[:, :, 0])
    for i in range(len(E)):
        f_occ = 2 / (np.exp((E[i] - Ef) * q / (kB * T)) + 1)
        rho += f_occ * np.abs(psi[:, :, i])**2
    return rho / (dx * dy)

def run_self_consistent_2d(
    Lx=300e-9, Ly=150e-9, Nx=100, Ny=50, m_eff=0.067, eps_r=12.9,
    Ef=0.01, T=4.0, V0=0.2, w=30e-9, sep=60e-9,
    max_iter=100, tol=1e-4, num_states=3
):
    dx, dy = Lx / Nx, Ly / Ny
    x = np.linspace(-Lx/2, Lx/2, Nx)
    y = np.linspace(-Ly/2, Ly/2, Ny)
    X, Y = np.meshgrid(x, y, indexing='ij')

    V_gate = triple_dot_potential(X, Y, V0=V0, w=w, sep=sep)
    V_old = V_gate.copy()

    for it in range(max_iter):
        E, psi = schrodinger_solver_2d(V_old, dx, dy, m_eff, num_states=num_states)
        rho = fermi_density(psi, E, Ef, T, dx, dy)
        phi = poisson_solver_2d(-q * rho, dx, dy, eps_r)
        V_new = V_gate + q * phi

        dV = np.linalg.norm(V_new - V_old)
        print(f"Iteration {it+1}: ΔV = {dV:.2e} V")

        if dV < tol:
            break
        V_old = V_new

    # Visualization
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    im0 = axes[0].imshow(V_new.T, extent=[x[0]*1e9,x[-1]*1e9,y[0]*1e9,y[-1]*1e9],
                         origin='lower', aspect='auto')
    axes[0].set_title("Self-Consistent Potential [eV]")
    fig.colorbar(im0, ax=axes[0])

    im1 = axes[1].imshow(rho.T / 1e6, extent=[x[0]*1e9,x[-1]*1e9,y[0]*1e9,y[-1]*1e9],
                         origin='lower', aspect='auto')
    axes[1].set_title("Electron Density [cm⁻³]")
    fig.colorbar(im1, ax=axes[1])
    plt.tight_layout()
    plt.show()

    return X, Y, V_new, rho, E, psi
if __name__ == "__main__":
    run_self_consistent_2d()
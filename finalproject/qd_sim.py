# import numpy as np
# from numpy.fft import fft2, ifft2
# from scipy.fft import dst, idst

# import scipy.sparse
# import scipy.sparse.linalg
# import matplotlib.pyplot as plt
# from hilbertcurve.hilbertcurve import HilbertCurve

# # Constants
# e = 1.602e-19
# epsilon_0 = 8.854e-12
# kB = 1.38e-23
# T = 4.2
# energy_offset1 = 0
# energy_offset2 = 0
# mu_lead = 0
# U = 100e-3 * e  # Interdot Coulomb interaction energy (example value)

# # Simulation grid (must be power of 2)
# Nx = 32
# Ny = 32
# L = 100e-9
# dx = L / (Nx - 1)
# dy = L / (Ny - 1)

# # Dot parameters
# dot_radius1 = 7e-9  # Increased radius for dot 1
# dot_radius2 = 7e-9  # Increased radius for dot 2
# dot1_center = (L / 4, L / 2)
# dot2_center = (3 * L / 4, L / 2)

# # Gate voltage grid (also power of 2)
# Nvg = 32
# V_range = np.linspace(-2, 2, Nvg)

# # Hilbert curve in gate voltage space
# p = int(np.log2(Nvg))
# vg_curve = HilbertCurve(p, 2)
# num_points = 2 ** (p * 2)
# grid_coords = [vg_curve.point_from_distance(d) for d in range(num_points)]

# grid_coords = [(i, j, V_range[i], V_range[j]) for i, j in grid_coords]

# # Storage
# n1_results = np.zeros((Nvg, Nvg))
# n2_results = np.zeros((Nvg, Nvg))

# # Initialize phi and rho for warm-starting
# phi = np.zeros((Nx, Ny))
# rho = np.zeros((Nx, Ny))

# # Simulation helpers
# def solve_poisson(phi, rho, max_iter=500, tol=1e-5):
#     for _ in range(max_iter):
#         phi_old = phi.copy()
#         phi[1:-1,1:-1] = 0.25 * (
#             phi[2:,1:-1] + phi[:-2,1:-1] + phi[1:-1,2:] + phi[1:-1,:-2]
#             + (dx**2 / epsilon_0) * rho[1:-1,1:-1]
#         )
#         if np.max(np.abs(phi - phi_old)) < tol:
#             break
#     return phi

# def solve_poisson_spectral(phi, rho, dx, Vg1, Vg2):
#     """
#     Solves Poisson's equation using the Discrete Sine Transform (DST)
#     for Dirichlet boundary conditions (adapted for non-zero boundaries).

#     Args:
#         phi (np.ndarray): Current electrostatic potential (Nx x Ny) - for warm-starting.
#         rho (np.ndarray): Charge density (Nx x Ny).
#         dx (float): Grid spacing.
#         Vg1 (float): Gate voltage on the left boundary (x=0).
#         Vg2 (float): Gate voltage on the right boundary (x=L).

#     Returns:
#         np.ndarray: Updated electrostatic potential (Nx x Ny).
#     """
#     Nx, Ny = rho.shape
#     Lx = (Nx - 1) * dx
#     Ly = (Ny - 1) * dx

#     # 1. Create x as a column vector for boundary potential
#     x = np.linspace(0, Lx, Nx).reshape(-1, 1)
#     boundary_potential_x_col = Vg1 + (Vg2 - Vg1) * x / Lx

#     # 2. Expand it to have the shape (Nx, Ny)
#     boundary_potential_x = np.tile(boundary_potential_x_col, (1, Ny))

#     # 3. Work with the deviation from the boundary potential
#     phi_prime = phi - boundary_potential_x
#     rho_prime = rho # We're solving for phi' such that laplacian(phi') = -rho/epsilon_0

#     # Apply DST to phi_prime and rho_prime along x
#     phi_prime_dst_x = dst(phi_prime[1:-1, :], axis=0, type=1)
#     rho_prime_dst_x = dst(rho_prime[1:-1, :], axis=0, type=1)

#     kx = np.pi * np.arange(1, Nx - 1) / Lx
#     KX, _ = np.meshgrid(kx, np.ones(Ny - 2)) # Shape (Nx-2, Ny-2) after second DST

#     # Apply DST along y
#     phi_prime_dst_xy = dst(phi_prime_dst_x[:, 1:-1], axis=1, type=1)
#     rho_prime_dst_xy = dst(rho_prime_dst_x[:, 1:-1], axis=1, type=1)

#     ky = np.pi * np.arange(1, Ny - 1) / Ly
#     _, KY = np.meshgrid(np.ones(Nx - 2), ky, indexing='ij') # Shape (Nx-2, Ny-2)

#     # 4. Solve in the transformed space
#     K_sq = KX**2 + KY**2
#     K_sq[K_sq == 0] = 1e-12

#     phi_prime_dst_xy = -rho_prime_dst_xy / (epsilon_0 * K_sq)

#     # 5. Perform inverse DST
#     phi_prime_interior = idst(idst(phi_prime_dst_xy, axis=1, type=1), axis=0, type=1)

#     # 6. Add back the boundary potential to get the full potential
#     phi_new = np.zeros((Nx, Ny))
#     phi_new[0, :] = Vg1
#     phi_new[-1, :] = Vg2
#     phi_new[1:-1, 1:-1] = boundary_potential_x[1:-1, 1:-1] + phi_prime_interior

#     return phi_new

# def fermi_dirac(E, mu):
#     return 1.0 if T == 0 and E <= mu else 1 / (np.exp((E - mu) / (kB * T)) + 1)

# def get_dot_potentials(phi):
#     i1, j1 = int(dot1_center[0] / dx), int(dot1_center[1] / dy)
#     i2, j2 = int(dot2_center[0] / dx), int(dot2_center[1] / dy)
#     return phi[i1, j1], phi[i2, j2]

# def compute_occupancy(phi1, phi2, n1_old, n2_old):
#     # Include the effect of the charge on the other dot in the chemical potential
#     mu1 = -e * phi1 + energy_offset1 + U * n2_old
#     mu2 = -e * phi2 + energy_offset2 + U * n1_old
#     return fermi_dirac(mu1, mu_lead), fermi_dirac(mu2, mu_lead)

# def update_charge_density(n1, n2):
#     rho = np.zeros((Nx, Ny))
#     dot_area1 = np.pi * dot_radius1**2
#     dot_area2 = np.pi * dot_radius2**2
#     charge1 = n1 * e / dot_area1
#     charge2 = n2 * e / dot_area2
#     for i in range(Nx):
#         for j in range(Ny):
#             x, y = i * dx, j * dy
#             if (x - dot1_center[0])**2 + (y - dot1_center[1])**2 <= dot_radius1**2:
#                 rho[i,j] = charge1
#             elif (x - dot2_center[0])**2 + (y - dot2_center[1])**2 <= dot_radius2**2:
#                 rho[i,j] = charge2
#     return rho

# # Main loop with warm starting
# for k, (i, j, Vg1, Vg2) in enumerate(grid_coords):
#     print(f"Simulating point {k+1}/{num_points}: Vg1={Vg1:.2f}, Vg2={Vg2:.2f}")
#     # Boundary conditions
#     phi[0,:] = Vg1
#     phi[-1,:] = Vg2

#     # Initialize occupancies for the self-consistent loop
#     n1 = 0.0
#     n2 = 0.0

#     # Self-consistent loop
#     for iter_count in range(20):  # Increased iterations for convergence
#         phi = solve_poisson_spectral(phi, rho, dx, Vg1, Vg2)
#         # phi = solve_poisson(phi, rho)

#         # phi = solve_poisson_cg(phi, rho, dx, Nx, Ny, boundary_values, max_iter=500, tol=1e-5)
#         phi1, phi2 = get_dot_potentials(phi)
#         n1_new, n2_new = compute_occupancy(phi1, phi2, n1, n2)
#         new_rho = update_charge_density(n1_new, n2_new)

#         if np.allclose(rho, new_rho, atol=1e-3):
#             n1 = n1_new
#             n2 = n2_new
#             print(f"  Converged in {iter_count+1} iterations: n1={n1:.2f}, n2={n2:.2f}")
#             break
#         rho = new_rho
#         n1 = n1_new
#         n2 = n2_new
#     else:
#         print(f"  Warning: Self-consistent loop did not converge for Vg1={Vg1:.2f}, Vg2={Vg2:.2f}")

#     n1_results[i, j] = n1
#     n2_results[i, j] = n2

# if __name__ == "__main__":
#     # Combine dot occupancies
#     n_total = n1_results + n2_results

#     # Plot the combined charge stability diagram
#     plt.figure(figsize=(6, 5))
#     plt.imshow(
#         n_total.T,  # Transpose for correct orientation
#         origin='lower',
#         extent=[V_range.min(), V_range.max(), V_range.min(), V_range.max()],
#         cmap='viridis',
#         aspect='auto'
#     )
#     plt.colorbar(label='Total Occupancy (n1 + n2)')
#     plt.xlabel('V_gate1 (V)')
#     plt.ylabel('V_gate2 (V)')
#     plt.title(f'Combined Dot Occupancy with U = {U/e:.3f} eV, Larger Dots')
#     plt.tight_layout()
#     plt.show()
import numpy as np
from numpy.fft import fft2, ifft2
from scipy.fft import dst, idst
import matplotlib.pyplot as plt
from hilbertcurve.hilbertcurve import HilbertCurve

# --- Constants ---
e = 1.602e-19
epsilon_0 = 8.854e-12
kB = 1.38e-23
T = 4.2
mu_lead = 0

# --- Simulation Grid ---
Nx = 32
Ny = 32
L = 100e-9
dx = L / (Nx - 1)
dy = L / (Ny - 1)

# --- Dot Parameters (for 2 dots as in the initial problem) ---
num_dots = 2
dot_radii = [7e-9, 7e-9]
dot_centers = [(L / 4, L / 2), (3 * L / 4, L / 2)]
energy_offsets = [0, 0]

# Inter-dot Coulomb interaction energy (example value)
U = 100e-3 * e
U_matrix = np.array([[0, U], [U, 0]])

# --- Gate Voltage Grid ---
Nvg = 32
V_range = np.linspace(-2, 2, Nvg)

# Hilbert curve in gate voltage space
p = int(np.log2(Nvg))
vg_curve = HilbertCurve(p, 2)
num_points = 2 ** (p * 2)
grid_coords = [vg_curve.point_from_distance(d) for d in range(num_points)]
grid_coords_vg = [(V_range[i], V_range[j]) for i, j in grid_coords]

# Storage
n_results = np.zeros((Nvg, Nvg, num_dots))

# Initialize phi and rho for warm-starting
phi = np.zeros((Nx, Ny))
rho = np.zeros((Nx, Ny))

# --- Simulation Helpers ---
def solve_poisson_spectral(phi, rho, dx, boundary_values):
    """
    Solves Poisson's equation using the Discrete Sine Transform (DST)
    for Dirichlet boundary conditions (adapted for non-zero boundaries).

    Args:
        phi (np.ndarray): Current electrostatic potential (Nx x Ny) - for warm-starting.
        rho (np.ndarray): Charge density (Nx x Ny).
        dx (float): Grid spacing.
        boundary_values (list): [Vg1, Vg2] - Gate voltages on the left and right boundaries.

    Returns:
        np.ndarray: Updated electrostatic potential (Nx x Ny).
    """
    Nx, Ny = rho.shape
    Lx = (Nx - 1) * dx
    Ly = (Ny - 1) * dx

    phi_prime = phi.copy()
    phi_prime[0, :] -= boundary_values[0]
    phi_prime[-1, :] -= boundary_values[1]

    rho_prime = rho

    phi_prime_dst_x = dst(phi_prime[1:-1, :], axis=0, type=1)
    rho_prime_dst_x = dst(rho_prime[1:-1, :], axis=0, type=1)

    kx = np.pi * np.arange(1, Nx - 1) / Lx
    KX, _ = np.meshgrid(kx, np.ones(Ny - 2))

    phi_prime_dst_xy = dst(phi_prime_dst_x[:, 1:-1], axis=1, type=1)
    rho_prime_dst_xy = dst(rho_prime_dst_x[:, 1:-1], axis=1, type=1)

    ky = np.pi * np.arange(1, Ny - 1) / Ly
    _, KY = np.meshgrid(np.ones(Nx - 2), ky, indexing='ij')

    K_sq = KX**2 + KY**2
    K_sq[K_sq == 0] = 1e-12

    phi_prime_dst_xy = -rho_prime_dst_xy / (epsilon_0 * K_sq)

    phi_prime_interior = idst(idst(phi_prime_dst_xy, axis=1, type=1), axis=0, type=1)

    phi_new = np.zeros((Nx, Ny))
    phi_new[0, :] = boundary_values[0]
    phi_new[-1, :] = boundary_values[1]
    phi_new[1:-1, 1:-1] = boundary_values[0] + (boundary_values[1] - boundary_values[0]) * np.linspace(0, 1, Nx - 2).reshape(-1, 1) + phi_prime_interior

    return phi_new

def fermi_dirac(E, mu):
    return 1.0 if T == 0 and E <= mu else 1 / (np.exp((E - mu) / (kB * T)) + 1)

def get_dot_potentials(phi, dot_centers):
    potentials = []
    for center in dot_centers:
        i, j = int(center[0] / dx), int(center[1] / dy)
        potentials.append(phi[i, j])
    return potentials

def compute_occupancies_n_dots(phi_dots, n_old, energy_offsets, U_matrix, mu_lead):
    """Computes the occupancy of N quantum dots."""
    N = len(phi_dots)
    n_new = np.zeros(N)
    for i in range(N):
        mu_i = -e * phi_dots[i] + energy_offsets[i]
        for j in range(N):
            if i != j:
                mu_i += U_matrix[i, j] * n_old[j]
        n_new[i] = fermi_dirac(mu_i, mu_lead)
    return n_new

def update_charge_density_n_dots(n, dot_centers, dot_radii):
    rho = np.zeros((Nx, Ny))
    for k in range(len(dot_centers)):
        center = dot_centers[k]
        radius = dot_radii[k]
        dot_area = np.pi * radius**2
        charge = n[k] * e / dot_area
        for i in range(Nx):
            for j in range(Ny):
                x, y = i * dx, j * dy
                if (x - center[0])**2 + (y - center[1])**2 <= radius**2:
                    rho[i, j] = charge
    return rho

# --- Main Loop ---
n_old = np.zeros(num_dots)  # Initialize occupancies

for k, (Vg1, Vg2) in enumerate(grid_coords_vg):
    print(f"Simulating point {k+1}/{num_points}: Vg1={Vg1:.2f}, Vg2={Vg2:.2f}")
    boundary_values = [Vg1, Vg2]

    # Initialize occupancies for the self-consistent loop
    n = np.zeros(num_dots)
    rho = np.zeros((Nx, Ny)) # Reset charge density

    # Self-consistent loop
    for iter_count in range(20):
        phi = solve_poisson_spectral(phi, rho, dx, boundary_values)
        phi_dots = get_dot_potentials(phi, dot_centers)
        n_new = compute_occupancies_n_dots(phi_dots, n, energy_offsets, U_matrix, mu_lead)
        new_rho = update_charge_density_n_dots(n_new, dot_centers, dot_radii)

        if np.allclose(n, n_new, atol=1e-3):
            n = n_new
            print(f"  Converged in {iter_count+1} iterations: n = {n}")
            break
        rho = new_rho
        n = n_new
    else:
        print(f"  Warning: Self-consistent loop did not converge for Vg1={Vg1:.2f}, Vg2={Vg2:.2f}")

    i, j = grid_coords[k]
    n_results[i, j, :] = n

# --- Visualization for 2 Dots ---
if __name__ == "__main__":
    n_total = np.sum(n_results, axis=-1)

    plt.figure(figsize=(6, 5))
    plt.imshow(
        n_total.T,  # Transpose for correct orientation
        origin='lower',
        extent=[V_range.min(), V_range.max(), V_range.min(), V_range.max()],
        cmap='viridis',
        aspect='auto'
    )
    plt.colorbar(label='Total Occupancy (n1 + n2)')
    plt.xlabel('V_gate1 (V)')
    plt.ylabel('V_gate2 (V)')
    plt.title(f'{num_dots} Dots - Total Occupancy with U = {U/e:.3f} eV, Larger Dots')
    plt.tight_layout()
    plt.show()

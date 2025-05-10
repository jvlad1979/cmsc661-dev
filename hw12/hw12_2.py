import numpy as np
import matplotlib.pyplot as plt

# --- Plotting Settings ---
plt.rcParams.update(
    {
        "text.usetex": True,  # Use LaTeX for text rendering
        "font.family": "serif",  # Set the font family to serif
        "font.serif": ["Times New Roman"],  # Specify the serif font
        "font.size": 12,  # Set the default font size - Adjusted for better readability in subplots
    }
)


# --- Problem Definition & Parameters ---
# Flux function for Burgers' equation: f(u) = 0.5 * u^2
def flux_burgers(u):
    return 0.5 * u**2


# Spatial domain
x_min = -2.0
x_max = 8.0
nx = 400  # Number of spatial points
dx = (x_max - x_min) / nx
# Using cell centers for numerical methods
x_coords = np.linspace(x_min + dx / 2, x_max - dx / 2, nx)

# Time domain and parameters
t_final = 6.0
cfl = 0.5  # Courant-Friedrichs-Lewy number
# Max characteristic speed for Burgers' with initial u in [0,1] is max(|u|) = 1
# dt = cfl * dx / max_speed
dt = cfl * dx / 1.0  # Assuming max speed is 1 initially

plot_times = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]

# Initial condition: u0(x) = 1 for x in [0,1], 0 otherwise
u_initial = np.zeros(nx)
u_initial[(x_coords >= 0) & (x_coords <= 1)] = 1.0


# --- Exact Solution ---
def exact_solution_burgers(x, t):
    """
    Computes the exact solution for the Burgers' equation with the given
    initial condition: u0(x) = 1 for x in [0,1], 0 otherwise.
    Based on Hyperbolic.pdf, Section 4.3.
    """
    u_exact = np.zeros_like(x)
    if t == 0:
        u_exact[(x >= 0) & (x <= 1)] = 1.0
        return u_exact

    # Shock position before interaction with rarefaction (t <= 2)
    shock_pos_early = 1 + 0.5 * t
    # Shock position after interaction (t > 2)
    shock_pos_late = np.sqrt(2 * t)

    for i, xi in enumerate(x):
        if t <= 2.0:
            if xi < 0:
                u_exact[i] = 0.0
            elif xi <= t:  # Rarefaction fan
                u_exact[i] = xi / t
            elif xi < shock_pos_early:  # Plateau
                u_exact[i] = 1.0
            else:  # Right of shock
                u_exact[i] = 0.0
        else:  # t > 2.0
            if xi < 0:
                u_exact[i] = 0.0
            elif xi <= shock_pos_late:  # Rarefaction fan up to the shock
                # The state in the rarefaction fan is x/t
                u_exact[i] = xi / t
            else:  # Right of shock
                u_exact[i] = 0.0
    return u_exact


# --- Numerical Methods ---


# Helper for Godunov: Determine u* based on uL and uR
def godunov_u_star(uL, uR):
    """
    Determines u* for the Godunov scheme for Burgers' equation.
    f(u) = 0.5*u^2, f'(u) = u. Sonic point u_s = 0.
    Based on Hyperbolic.pdf, Section 5.6, Cases 1-4.
    """
    # Case 1: f'(uL) >= 0 and f'(uR) >= 0 (i.e., uL >= 0 and uR >= 0)
    if uL >= 0 and uR >= 0:
        return uL
    # Case 2: f'(uL) <= 0 and f'(uR) <= 0 (i.e., uL <= 0 and uR <= 0)
    elif uL <= 0 and uR <= 0:
        return uR
    # Case 3: f'(uL) >= 0 >= f'(uR) (i.e., uL >= 0 and uR <= 0)
    elif uL >= 0 and uR <= 0:
        # Shock speed s = 0.5 * (uL + uR) - Note: This is the Rankine-Hugoniot speed.
        # u* is uL if shock speed > 0, uR if shock speed < 0.
        s = 0.5 * (uL + uR)
        if s > 1e-9:  # Using a small tolerance for comparison with 0
            return uL
        elif s < -1e-9:  # Using a small tolerance for comparison with 0
            return uR
        else:  # s is close to 0
            # This happens when uL is close to -uR.
            # For Burgers, if uL > 0 and uR < 0 and uL+uR=0, the state at x/t=0 is u=0 (sonic point)
            return 0.0  # Sonic point value when shock is stationary
    # Case 4: f'(uL) < 0 < f'(uR) (i.e., uL < 0 and uR > 0) -> Transonic rarefaction
    elif uL < 0 and uR > 0:
        return 0.0  # u_s (sonic point)

    # Fallback for any unhandled theoretical edge cases
    # print(f"Warning: Godunov u* unhandled case: uL={uL}, uR={uR}")
    # Using a simplified Oleinik-like condition based on the sign of u_L + u_R
    if uL + uR > 0:
        return uL
    else:
        return uR


# Generic solver function
def solve_method(u_init, dx, dt, t_final_sim, method_func, current_plot_times):
    u = u_init.copy()
    nx_sim = len(u)
    t = 0.0

    results = {}
    next_plot_idx = 0

    # Store initial condition
    if next_plot_idx < len(current_plot_times) and np.isclose(
        t, current_plot_times[next_plot_idx], atol=dt / 10
    ):
        results[current_plot_times[next_plot_idx]] = u.copy()
        next_plot_idx += 1

    # Time evolution
    while t < t_final_sim:
        u = method_func(u, dt, dx)
        t += dt

        # Store solution if close to a plot time
        if (
            next_plot_idx < len(current_plot_times)
            and t >= current_plot_times[next_plot_idx] - dt / 10
        ):
            # Find the correct plot time index in case we stepped over one or more
            while (
                next_plot_idx < len(current_plot_times)
                and t >= current_plot_times[next_plot_idx] - dt / 10
            ):
                # Only store if the time hasn't been stored or is significantly closer
                if current_plot_times[next_plot_idx] not in results or abs(
                    t - current_plot_times[next_plot_idx]
                ) < abs(
                    t
                    - (
                        list(results.keys())[-1]
                        if results
                        else -np.inf  # Compare with the last stored time
                    )
                ):
                    results[current_plot_times[next_plot_idx]] = u.copy()

                next_plot_idx += 1

    # Ensure the final time is captured if it's a plot time
    if t_final_sim in current_plot_times and t_final_sim not in results:
        results[t_final_sim] = u.copy()

    return results


# --- Implementations of specific methods using cell averages ---


# Lax-Friedrichs Method
def method_lax_friedrichs(u_n, dt, dx):
    nx_sim = len(u_n)
    u_np1 = np.zeros(nx_sim)
    # Use ghost cells for boundary conditions (assuming u=0 outside domain)
    u_padded = np.pad(u_n, 1, mode="constant", constant_values=0.0)
    f_padded = flux_burgers(u_padded)

    for j in range(0, nx_sim):
        # Corresponds to index j in the original array, which is j+1 in the padded array
        u_np1[j] = 0.5 * (u_padded[j] + u_padded[j + 2]) - 0.5 * dt / dx * (
            f_padded[j + 2] - f_padded[j]
        )

    return u_np1


# Richtmyer Method (Two-step Lax-Wendroff)
def method_richtmyer(u_n, dt, dx):
    nx_sim = len(u_n)
    u_np1 = np.zeros(nx_sim)
    u_half = np.zeros(nx_sim - 1)  # Stores U_{j+1/2}^{n+1/2}

    # Predictor step at half-steps (staggered grid)
    # Using ghost cells assuming u=0 outside domain
    u_padded = np.pad(u_n, 1, mode="constant", constant_values=0.0)
    f_padded = flux_burgers(u_padded)
    for j in range(nx_sim - 1):
        # Corresponds to U_{j+1/2}^{n+1/2} using u_j and u_{j+1} (indices j+1 and j+2 in padded)
        u_half[j] = 0.5 * (u_padded[j + 1] + u_padded[j + 2]) - 0.5 * dt / dx * (
            f_padded[j + 2] - f_padded[j + 1]
        )

    # Corrector step at full steps
    # Using ghost cells for u_half (assuming u_half=0 outside its domain)
    u_half_padded = np.pad(u_half, 1, mode="constant", constant_values=0.0)
    f_half_padded = flux_burgers(u_half_padded)
    for j in range(0, nx_sim):
        # Corresponds to U_j^{n+1} using u_{j-1/2} and u_{j+1/2} (indices j and j+1 in padded u_half)
        u_np1[j] = u_n[j] - dt / dx * (f_half_padded[j + 1] - f_half_padded[j])

    return u_np1


# MacCormack Method
def method_maccormack(u_n, dt, dx):
    nx_sim = len(u_n)
    u_np1 = np.zeros(nx_sim)
    u_star = np.zeros(nx_sim)

    # Predictor step (forward difference)
    # Using ghost cells for forward difference (assuming u=0 to the right)
    u_padded_fwd = np.pad(u_n, (0, 1), mode="constant", constant_values=0.0)
    f_padded_fwd = flux_burgers(u_padded_fwd)
    for j in range(0, nx_sim):
        u_star[j] = u_n[j] - dt / dx * (f_padded_fwd[j + 1] - f_padded_fwd[j])

    # Corrector step (backward difference)
    # Using ghost cells for backward difference (assuming u_star=0 to the left)
    u_star_padded_bwd = np.pad(u_star, (1, 0), mode="constant", constant_values=0.0)
    f_star_padded_bwd = flux_burgers(u_star_padded_bwd)
    for j in range(0, nx_sim):
        # Corresponds to index j in original array, which is j+1 in the padded array
        u_np1[j] = 0.5 * (u_n[j] + u_star[j]) - 0.5 * dt / dx * (
            f_star_padded_bwd[j + 1] - f_star_padded_bwd[j]
        )

    return u_np1


# Godunov Method
def method_godunov(u_n, dt, dx):
    nx_sim = len(u_n)
    u_np1 = np.zeros(nx_sim)
    # Numerical flux at cell interfaces (j+1/2)
    fluxes_interface = np.zeros(nx_sim + 1)

    # Using ghost cells (assuming u=0 outside domain) to calculate fluxes at boundaries
    u_padded = np.pad(u_n, 1, mode="constant", constant_values=0.0)

    # Calculate intercell fluxes F_{j+1/2} for j from -1 to nx-1 (using padded indices)
    # F_j+1/2 corresponds to fluxes_interface[j+1]
    for j in range(nx_sim + 1):  # j from 0 to nx
        uL = u_padded[j]
        uR = u_padded[j + 1]
        u_star_val = godunov_u_star(uL, uR)
        fluxes_interface[j] = flux_burgers(u_star_val)

    # Update cell averages (from Hyperbolic.pdf, equation 62)
    # U_j^{n+1} = U_j^n - dt/dx * (F_{j+1/2} - F_{j-1/2})
    for j in range(0, nx_sim):
        u_np1[j] = u_n[j] - dt / dx * (fluxes_interface[j + 1] - fluxes_interface[j])

    return u_np1


if __name__ == "__main__":
    # --- Run Solvers ---
    print(f"dx = {dx}, dt = {dt}, CFL = {dt / dx}")

    # Pass the correct time points and use the general solver
    results_lf = solve_method(
        u_initial, dx, dt, t_final, method_lax_friedrichs, plot_times
    )
    results_richtmyer = solve_method(
        u_initial, dx, dt, t_final, method_richtmyer, plot_times
    )
    results_maccormack = solve_method(
        u_initial, dx, dt, t_final, method_maccormack, plot_times
    )
    results_godunov = solve_method(
        u_initial, dx, dt, t_final, method_godunov, plot_times
    )

    # --- Plotting ---
    num_plots = len(plot_times)
    cols = 4  # Number of columns for subplots
    rows = 2  # Number of rows for subplots

    # Calculate figure height to accommodate 2 rows and width for 4 columns
    fig, axes = plt.subplots(
        rows, cols, figsize=(18, 9), squeeze=False
    )  # Adjusted figsize
    axes_flat = axes.flatten()  # Flatten for easy iteration

    for i, t_plot in enumerate(plot_times):
        # Ensure we don't exceed the number of axes
        if i >= rows * cols:
            print(
                f"Warning: More plot times ({num_plots}) than available subplots ({rows * cols}). Skipping remaining plot times."
            )
            break

        ax = axes_flat[i]

        # Exact solution
        u_exact_plot = exact_solution_burgers(x_coords, t_plot)
        ax.plot(x_coords, u_exact_plot, "k-", label="Exact", linewidth=2)

        # Numerical solutions - Plot if the time was captured by the solver
        if t_plot in results_lf:
            ax.plot(
                x_coords,
                results_lf[t_plot],
                "r.--",
                label="Lax-Friedrichs",
                markersize=3,
                alpha=0.7,
            )
        if t_plot in results_richtmyer:
            ax.plot(
                x_coords,
                results_richtmyer[t_plot],
                "b.--",
                label="Richtmyer",
                markersize=3,
                alpha=0.7,
            )
        if t_plot in results_maccormack:
            ax.plot(
                x_coords,
                results_maccormack[t_plot],
                "g.--",
                label="MacCormack",
                markersize=3,
                alpha=0.7,
            )
        if t_plot in results_godunov:
            ax.plot(
                x_coords,
                results_godunov[t_plot],
                "m.--",
                label="Godunov",
                markersize=3,
                alpha=0.7,
            )

        ax.set_title(f"t = {t_plot:.1f}")  # Simplified title
        ax.set_xlabel("$x$")
        ax.set_ylabel("$u$")  # Simplified label
        ax.set_xlim([x_min, x_max])
        ax.set_ylim([-0.2, 1.2])  # Adjusted y-limit for better visualization
        ax.grid(True, linestyle=":", alpha=0.6)
        if i == 0:  # Add legend to the first plot
            ax.legend(loc="upper right")

    # Hide any unused subplots in the 2x4 grid (there should be one unused axis)
    for j in range(num_plots, rows * cols):
        fig.delaxes(axes_flat[j])

    # Adjust layout to prevent overlapping titles/labels and center the subplots
    plt.tight_layout()
    plt.show()

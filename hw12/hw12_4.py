import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft, fftfreq

plt.rcParams.update(
    {
        "text.usetex": True,  # Use LaTeX for text rendering
        "font.family": "serif",  # Set the font family to serif
        "font.serif": ["Times New Roman"],  # Specify the serif font
        "font.size": 22,  # Set the default font size
    }
)

# --- Problem Parameters ---
N_POINTS = 128  # Number of spatial grid points
L_DOMAIN = 32 * np.pi  # Length of the physical spatial domain for x
L0_SCALE = 16.0  # Scaling factor for x, so X = x/L0 is in [0, 2*pi]

# Time parameters
T_FINAL = 200.0
DT_INITIAL = 0.1  # Initial and maximum time step
DT_MIN_ALLOWED = 1e-6  # Minimum allowed time step to prevent infinite loops
DT_REDUCTION_FACTOR = 2.0  # Factor to reduce dt by on failure
DT_INCREASE_FACTOR = 1.1  # Factor to increase dt by on success (cautiously)
MAX_STEP_RETRIES = 10  # Max retries for a single time point before aborting

# --- Function Definitions ---


def initialize_system(n_pts, l_domain, l0_scale):
    """Initializes spatial grids, wavenumbers, and linear operator symbol."""
    x_X_grid = np.linspace(0, 2 * np.pi, n_pts, endpoint=False)
    x_physical_grid = np.linspace(0, l_domain, n_pts, endpoint=False)
    # Wavenumbers in standard FFT order: 0, 1, ..., N/2-1, -N/2, ..., -1
    k_integers_fft_order = fftfreq(n_pts) * n_pts

    K_physical_wavenumbers = k_integers_fft_order / l0_scale  # Still in FFT order
    linear_op_symbol = (K_physical_wavenumbers**2) - (K_physical_wavenumbers**4)
    return (
        x_X_grid,
        x_physical_grid,
        K_physical_wavenumbers,
        linear_op_symbol,
        k_integers_fft_order,
    )


def get_initial_condition(x_X_grid):
    """Computes the initial condition in X-space and its FFT."""
    u0_X_physical_vals = np.cos(x_X_grid) * (1 + np.sin(x_X_grid))
    u0_hat_vals = fft(u0_X_physical_vals)
    return u0_X_physical_vals, u0_hat_vals


def rhs_ks_etdrk4(
    u_hat_stage_vals,
    dt_for_stage_factors,
    L_op_symb_stage,
    K_phys_stage,
    k_int_fft_order_stage,
    n_pts_stage,
):
    """
    Computes the transformed nonlinear term for a stage of the ETDRK4 scheme,
    including 2/3 dealiasing rule.
    k_int_fft_order_stage: integer wavenumbers in standard FFT order [0, 1, ..., -N/2, ...]
    n_pts_stage: N_POINTS
    """
    exp_L_dt_stage = np.exp(L_op_symb_stage * dt_for_stage_factors)
    exp_neg_L_dt_stage = np.exp(-L_op_symb_stage * dt_for_stage_factors)

    u_prime_hat = exp_L_dt_stage * u_hat_stage_vals
    u_prime_physical = ifft(u_prime_hat)

    u_prime_sq_physical = u_prime_physical**2
    u_prime_sq_hat = fft(u_prime_sq_physical)

    # Apply 2/3 dealiasing rule
    # Zero out modes |k| >= N/3.
    # k_int_fft_order_stage contains wavenumbers like [0, 1, ..., N/2-1, -N/2, ..., -1]
    dealias_cutoff_k = n_pts_stage / 3.0
    dealias_mask = np.abs(k_int_fft_order_stage) < dealias_cutoff_k
    u_prime_sq_hat_dealiased = u_prime_sq_hat * dealias_mask

    N_u_prime_hat = (
        -0.5 * 1j * K_phys_stage * u_prime_sq_hat_dealiased
    )  # Use dealiased version
    return exp_neg_L_dt_stage * N_u_prime_hat


def run_simulation_adaptive(
    u_hat_initial,
    t_final_sim,
    dt_initial_sim,
    L_op_symbol_sim,
    K_physical_sim,
    k_int_fft_order_sim,
    n_pts_sim,
):
    """Runs the ETDRK4 time-stepping simulation with adaptive dt and dealiasing."""
    u_hat = u_hat_initial.copy()

    solution_data_u_physical_sim = []
    time_points_sim = []

    u0_physical_for_storage = ifft(u_hat_initial)
    solution_data_u_physical_sim.append(np.real(u0_physical_for_storage).copy())
    time_points_sim.append(0.0)

    t_current_sim = 0.0
    dt_current = dt_initial_sim

    plot_interval_sim = 1.0
    next_plot_time_sim = plot_interval_sim
    if next_plot_time_sim > t_final_sim + dt_current / 2:
        next_plot_time_sim = t_final_sim

    print(
        f"Starting adaptive simulation with dealiasing for Kuramoto-Sivashinsky equation (ETDRK4)."
    )
    print(f"N = {n_pts_sim}, t_final = {t_final_sim}, Initial dt = {dt_initial_sim}")

    steps_taken = 0
    consecutive_successes = 0
    total_failed_attempts = 0

    while t_current_sim < t_final_sim:
        if t_current_sim + dt_current > t_final_sim:
            dt_current = t_final_sim - t_current_sim
            if dt_current <= 1e-12:
                break

        exp_L_dt_full_step = np.exp(L_op_symbol_sim * dt_current)
        u_hat_before_step = u_hat.copy()

        step_successful = False
        retries_for_this_logical_step = 0

        while not step_successful and retries_for_this_logical_step < MAX_STEP_RETRIES:
            try:
                k1 = rhs_ks_etdrk4(
                    u_hat_before_step,
                    0.0,
                    L_op_symbol_sim,
                    K_physical_sim,
                    k_int_fft_order_sim,
                    n_pts_sim,
                )

                k2_input_u_hat = u_hat_before_step + 0.5 * dt_current * k1
                k2 = rhs_ks_etdrk4(
                    k2_input_u_hat,
                    0.5 * dt_current,
                    L_op_symbol_sim,
                    K_physical_sim,
                    k_int_fft_order_sim,
                    n_pts_sim,
                )

                k3_input_u_hat = u_hat_before_step + 0.5 * dt_current * k2
                k3 = rhs_ks_etdrk4(
                    k3_input_u_hat,
                    0.5 * dt_current,
                    L_op_symbol_sim,
                    K_physical_sim,
                    k_int_fft_order_sim,
                    n_pts_sim,
                )

                k4_input_u_hat = u_hat_before_step + dt_current * k3
                k4 = rhs_ks_etdrk4(
                    k4_input_u_hat,
                    dt_current,
                    L_op_symbol_sim,
                    K_physical_sim,
                    k_int_fft_order_sim,
                    n_pts_sim,
                )

                u_hat_intermediate = u_hat_before_step + (dt_current / 6.0) * (
                    k1 + 2 * k2 + 2 * k3 + k4
                )
                u_hat_new_attempt = exp_L_dt_full_step * u_hat_intermediate

                if np.any(np.isnan(u_hat_new_attempt)) or np.any(
                    np.isinf(u_hat_new_attempt)
                ):
                    # This specific error message helps pinpoint where NaN/Inf originates
                    nan_in_k1 = np.any(np.isnan(k1)) or np.any(np.isinf(k1))
                    nan_in_k2 = np.any(np.isnan(k2)) or np.any(np.isinf(k2))
                    nan_in_k3 = np.any(np.isnan(k3)) or np.any(np.isinf(k3))
                    nan_in_k4 = np.any(np.isnan(k4)) or np.any(np.isinf(k4))
                    nan_in_inter = np.any(np.isnan(u_hat_intermediate)) or np.any(
                        np.isinf(u_hat_intermediate)
                    )
                    raise ValueError(
                        f"NaN/Inf in u_hat_new. k1:{nan_in_k1}, k2:{nan_in_k2}, k3:{nan_in_k3}, k4:{nan_in_k4}, inter:{nan_in_inter}"
                    )

                u_hat = u_hat_new_attempt
                step_successful = True

            except (OverflowError, ValueError) as e:
                retries_for_this_logical_step += 1
                total_failed_attempts += 1
                dt_old = dt_current
                dt_current /= DT_REDUCTION_FACTOR
                # Recalculate propagator for the new, smaller dt_current for the retry
                exp_L_dt_full_step = np.exp(L_op_symbol_sim * dt_current)
                print(
                    f"--- Step failed at t={t_current_sim:.3f} with dt={dt_old:.2e}. Error: {type(e).__name__}. Retrying with dt={dt_current:.2e} (Retry {retries_for_this_logical_step}/{MAX_STEP_RETRIES}) ---"
                )

                if dt_current < DT_MIN_ALLOWED:
                    print(
                        f"!!! FATAL: dt_current ({dt_current:.2e}) below minimum allowed ({DT_MIN_ALLOWED:.2e}). Aborting. !!!"
                    )
                    return np.array(solution_data_u_physical_sim), np.array(
                        time_points_sim
                    )

        if not step_successful:
            print(
                f"!!! FATAL: Failed to complete step starting at t={t_current_sim:.3f} after {MAX_STEP_RETRIES} retries. Aborting. !!!"
            )
            return np.array(solution_data_u_physical_sim), np.array(time_points_sim)

        t_current_sim += dt_current
        steps_taken += 1

        if retries_for_this_logical_step == 0:
            consecutive_successes += 1
            if (
                consecutive_successes >= 10
            ):  # Increase dt after 10 consecutive successes
                dt_new = min(dt_current * DT_INCREASE_FACTOR, dt_initial_sim)
                if dt_new > dt_current + 1e-9:  # Ensure meaningful increase
                    dt_current = dt_new
                consecutive_successes = 0
        else:
            consecutive_successes = 0

        if t_current_sim >= next_plot_time_sim - dt_current / 2.0 or np.isclose(
            t_current_sim, t_final_sim
        ):
            if not time_points_sim or not np.isclose(
                time_points_sim[-1], t_current_sim, atol=dt_current / 2
            ):
                u_physical_plot = ifft(u_hat)
                solution_data_u_physical_sim.append(np.real(u_physical_plot).copy())
                time_points_sim.append(t_current_sim)

                if steps_taken % (
                    200 if dt_initial_sim > 0.05 else 800
                ) == 0 or np.isclose(t_current_sim, t_final_sim):
                    print(
                        f"Stored data at t = {t_current_sim:.3f}. Current dt = {dt_current:.2e}. Total steps: {steps_taken}. Failed attempts: {total_failed_attempts}."
                    )

            if (
                t_current_sim >= next_plot_time_sim - dt_current / 2.0
                and t_current_sim < t_final_sim - dt_current / 1.9
            ):
                next_plot_time_sim += plot_interval_sim
                if next_plot_time_sim > t_final_sim:
                    next_plot_time_sim = t_final_sim

    print(
        f"Simulation finished. Last computed time: {t_current_sim:.3f}. Total steps taken: {steps_taken}. Total failed attempts: {total_failed_attempts}."
    )
    return np.array(solution_data_u_physical_sim), np.array(time_points_sim)


def plot_results(
    solution_matrix_plot,
    time_points_plot,
    x_physical_plot,
    n_pts_plot,
    dt_initial_plot,
    t_final_plot,
):
    """Plots the simulation results."""
    if solution_matrix_plot.shape[0] <= 1 or len(time_points_plot) <= 1:
        print("Not enough data points available for plotting.")
        if len(time_points_plot) > 0:
            print(f"Collected time points: {time_points_plot}")
        return

    plt.figure(figsize=(10, 7))
    custom_vmin = -3.0
    custom_vmax = 3.0

    plot_t_min = time_points_plot.min()
    plot_t_max = time_points_plot.max()

    print(f"Plotting data from t={plot_t_min:.3f} to t={plot_t_max:.3f}")
    print(f"Solution matrix shape: {solution_matrix_plot.shape}")
    print(f"Number of time points for plot: {len(time_points_plot)}")

    plt.imshow(
        solution_matrix_plot,
        aspect="auto",
        origin="lower",
        extent=[x_physical_plot.min(), x_physical_plot.max(), plot_t_min, plot_t_max],
        cmap="viridis",
        vmin=custom_vmin,
        vmax=custom_vmax,
    )

    plt.colorbar(label="u(x,t)")
    plt.xlabel("x (Physical domain)")
    plt.ylabel("t (Time)")
    title_latex = r"$u_t + u_{xxxx} + u_{xx} + \frac{1}{2}(u^2)_x = 0$"
    plt.title(
        f"Kuramoto-Sivashinsky Equation \n{title_latex}\n(N={n_pts_plot}, Initial dt={dt_initial_plot:.2e})"
    )
    plt.ylim(0, t_final_plot)
    plt.show()
    plt.savefig("kuramoto.png", dpi=300)


# --- Main Execution Block ---
if __name__ == "__main__":
    (
        x_X_main,
        x_physical_main,
        K_physical_main,
        L_op_symbol_main,
        k_int_fft_order_main,
    ) = initialize_system(N_POINTS, L_DOMAIN, L0_SCALE)

    _, u_hat_initial_main = get_initial_condition(x_X_main)

    solution_matrix_main, time_points_main = run_simulation_adaptive(
        u_hat_initial_main,
        T_FINAL,
        DT_INITIAL,
        L_op_symbol_main,
        K_physical_main,
        k_int_fft_order_main,  # Pass the integer wavenumbers for dealiasing
        N_POINTS,
    )

    plot_results(
        solution_matrix_main,
        time_points_main,
        x_physical_main,
        N_POINTS,
        DT_INITIAL,
        T_FINAL,
    )

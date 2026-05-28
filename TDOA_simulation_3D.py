"""
Basic 3D TDOA localization spacing simulation.

This script analyzes how the distance between sensors affects localization
error for a 7-sensor 3D TDOA array. Sensor 1 is used as the reference sensor.
"""

from pathlib import Path

import matplotlib

# Prefer Tk windows for interactive chart/table display and avoid Qt backend
# issues. If Tk is not available, the script still saves CSV/PNG normally.
try:
    import tkinter as tk
    from tkinter import ttk

    matplotlib.use("TkAgg")
    GUI_AVAILABLE = True
except Exception:
    matplotlib.use("Agg")
    GUI_AVAILABLE = False

from matplotlib.figure import Figure
import numpy as np
import pandas as pd
from scipy.optimize import least_squares


# ---- System parameters ----
SOUND_SPEED_MPS = 343.0  # m/s, speed of sound in air
PULSE_WIDTH_S = 100e-9  # s, 100 ns pulse width. Kept for traceability.
SAMPLING_RATE_HZ = 100_000.0  # Hz
TIMING_NOISE_STD_S = 1.0 / (np.sqrt(12.0) * SAMPLING_RATE_HZ)  # s
MONTE_CARLO_TRIALS = 200

CSV_FILENAME = "tdoa_spacing_result.csv"
PLOT_FILENAME = "tdoa_error_vs_distance.png"


def create_sensor_array(spacing_m):
    """
    Create the fixed 7-sensor 3D array.

    Parameters
    ----------
    spacing_m : float
        Sensor spacing in meters.

    Returns
    -------
    np.ndarray
        Sensor coordinates with shape (7, 3), unit: meter.
    """
    return np.array(
        [
            [0.0, 0.0, 0.0],
            [spacing_m, 0.0, 0.0],
            [2.0 * spacing_m, 0.0, 0.0],
            [0.0, spacing_m, 0.0],
            [0.0, 2.0 * spacing_m, 0.0],
            [0.0, 0.0, spacing_m],
            [0.0, 0.0, 2.0 * spacing_m],
        ],
        dtype=float,
    )


def generate_source_position(distance_m):
    """
    Generate source position along the (1, 1, 1) direction.

    Parameters
    ----------
    distance_m : float
        Source distance from the origin in meters.

    Returns
    -------
    np.ndarray
        Source coordinate [x, y, z], unit: meter.
    """
    unit_direction = np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0)
    return distance_m * unit_direction


def compute_toa(source_position, sensors, sound_speed):
    """
    Compute time of arrival for each sensor.

    Parameters
    ----------
    source_position : np.ndarray
        Source coordinate [x, y, z], unit: meter.
    sensors : np.ndarray
        Sensor coordinates with shape (N, 3), unit: meter.
    sound_speed : float
        Sound speed, unit: m/s.

    Returns
    -------
    np.ndarray
        TOA for each sensor, unit: second.
    """
    distances = np.linalg.norm(source_position - sensors, axis=1)
    return distances / sound_speed


def compute_tdoa(toa):
    """
    Compute TDOA values using sensor 1 as the reference.

    Parameters
    ----------
    toa : np.ndarray
        TOA for all sensors, unit: second.

    Returns
    -------
    np.ndarray
        TDOA for sensors 2 to 7 relative to sensor 1, unit: second.
    """
    return toa[1:] - toa[0]


def add_tdoa_noise(tdoa, timing_noise_std):
    """
    Add independent Gaussian timing noise to each TDOA measurement.

    Parameters
    ----------
    tdoa : np.ndarray
        Clean TDOA values, unit: second.
    timing_noise_std : float
        Timing noise standard deviation, unit: second.

    Returns
    -------
    np.ndarray
        Noisy TDOA values, unit: second.
    """
    noise = np.random.normal(loc=0.0, scale=timing_noise_std, size=tdoa.shape)
    return tdoa + noise


def tdoa_residual(position, sensors, measured_tdoa, sound_speed):
    """
    Residual function for nonlinear least squares.

    residual_i = predicted_tdoa_i1(position) - measured_tdoa_i1

    Parameters
    ----------
    position : np.ndarray
        Candidate source coordinate [x, y, z], unit: meter.
    sensors : np.ndarray
        Sensor coordinates with shape (7, 3), unit: meter.
    measured_tdoa : np.ndarray
        Measured TDOA values for sensors 2 to 7, unit: second.
    sound_speed : float
        Sound speed, unit: m/s.

    Returns
    -------
    np.ndarray
        Residual vector, unit: second.
    """
    ranges = np.linalg.norm(position - sensors, axis=1)
    predicted_tdoa = (ranges[1:] - ranges[0]) / sound_speed
    return predicted_tdoa - measured_tdoa


def estimate_position(sensors, measured_tdoa, sound_speed, initial_guess):
    """
    Estimate source position with scipy.optimize.least_squares.

    Parameters
    ----------
    sensors : np.ndarray
        Sensor coordinates with shape (7, 3), unit: meter.
    measured_tdoa : np.ndarray
        Measured TDOA values for sensors 2 to 7, unit: second.
    sound_speed : float
        Sound speed, unit: m/s.
    initial_guess : np.ndarray
        Initial source position guess [x, y, z], unit: meter.

    Returns
    -------
    np.ndarray
        Estimated source coordinate [x, y, z], unit: meter.

    Raises
    ------
    RuntimeError
        If least_squares does not converge to a finite solution.
    """
    result = least_squares(
        tdoa_residual,
        x0=initial_guess,
        args=(sensors, measured_tdoa, sound_speed),
        method="lm",
        max_nfev=1000,
    )

    if not result.success or not np.all(np.isfinite(result.x)):
        raise RuntimeError(result.message)

    return result.x


def run_simulation(spacing_cm):
    """
    Run the full Monte Carlo simulation.

    Parameters
    ----------
    spacing_cm : float
        Sensor spacing in centimeters.

    Returns
    -------
    pd.DataFrame
        Summary table for each source distance.
    """
    spacing_m = spacing_cm / 100.0
    sensors = create_sensor_array(spacing_m)
    rows = []

    print("\nSensor coordinates (meter):")
    for idx, sensor in enumerate(sensors, start=1):
        print(
            f"sensor {idx}: "
            f"({sensor[0]:.6f}, {sensor[1]:.6f}, {sensor[2]:.6f})"
        )

    for distance_m in range(10, 101, 10):
        true_position = generate_source_position(float(distance_m))
        clean_toa = compute_toa(true_position, sensors, SOUND_SPEED_MPS)
        clean_tdoa = compute_tdoa(clean_toa)

        errors_m = []
        failure_count = 0

        for _ in range(MONTE_CARLO_TRIALS):
            measured_tdoa = add_tdoa_noise(clean_tdoa, TIMING_NOISE_STD_S)

            # This intentionally starts close to the truth so this basic model
            # focuses on TDOA noise and sensor geometry, not global search.
            initial_guess = true_position + np.random.normal(0.0, 1.0, size=3)

            try:
                estimated_position = estimate_position(
                    sensors=sensors,
                    measured_tdoa=measured_tdoa,
                    sound_speed=SOUND_SPEED_MPS,
                    initial_guess=initial_guess,
                )
            except RuntimeError:
                failure_count += 1
                continue

            error_m = np.linalg.norm(estimated_position - true_position)
            if np.isfinite(error_m):
                errors_m.append(error_m)
            else:
                failure_count += 1

        errors_m = np.array(errors_m, dtype=float)

        if errors_m.size == 0:
            mean_error_m = np.nan
            rmse_m = np.nan
            median_error_m = np.nan
            p95_error_m = np.nan
            relative_rmse_percent = np.nan
        else:
            mean_error_m = float(np.mean(errors_m))
            rmse_m = float(np.sqrt(np.mean(errors_m**2)))
            median_error_m = float(np.median(errors_m))
            p95_error_m = float(np.percentile(errors_m, 95))
            relative_rmse_percent = float((rmse_m / distance_m) * 100.0)

        rows.append(
            {
                "distance_m": float(distance_m),
                "true_x": float(true_position[0]),
                "true_y": float(true_position[1]),
                "true_z": float(true_position[2]),
                "mean_error_m": mean_error_m,
                "rmse_m": rmse_m,
                "median_error_m": median_error_m,
                "p95_error_m": p95_error_m,
                "relative_rmse_percent": relative_rmse_percent,
                "failure_count": failure_count,
            }
        )

    return pd.DataFrame(rows)


def create_error_figure(df, spacing_m, figsize=(9, 6), dpi=100):
    """
    Create the RMSE and mean localization error figure.

    Parameters
    ----------
    df : pd.DataFrame
        Simulation result table.
    spacing_m : float
        Sensor spacing in meters.
    figsize : tuple
        Figure size in inches.
    dpi : int
        Figure DPI.

    Returns
    -------
    matplotlib.figure.Figure
        Figure object containing the error curves.
    """
    fig = Figure(figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111)
    ax.plot(
        df["distance_m"],
        df["rmse_m"],
        marker="o",
        linewidth=2,
        label="RMSE error",
    )
    ax.plot(
        df["distance_m"],
        df["mean_error_m"],
        marker="s",
        linewidth=2,
        label="Mean error",
    )
    ax.set_xlabel("Source distance (m)")
    ax.set_ylabel("Localization error (m)")
    ax.set_title(f"TDOA Localization Error vs Distance, spacing = {spacing_m:g} m")
    ax.legend()
    ax.grid(True, alpha=0.35)
    fig.tight_layout()
    return fig


def plot_results(df, spacing_m):
    """
    Plot RMSE and mean localization error versus source distance, then save PNG.

    Parameters
    ----------
    df : pd.DataFrame
        Simulation result table.
    spacing_m : float
        Sensor spacing in meters.
    """
    fig = create_error_figure(df, spacing_m)

    output_path = Path(__file__).resolve().parent / PLOT_FILENAME
    fig.savefig(output_path, dpi=300)


def _format_table_value(value):
    """Format table values for display in the Tkinter result window."""
    if pd.isna(value):
        return "nan"
    if isinstance(value, (int, np.integer)):
        return str(value)
    if isinstance(value, (float, np.floating)):
        return f"{value:.6g}"
    return str(value)


def show_results_windows(df, spacing_m):
    """
    Show separate chart and table windows.

    The implementation uses Tkinter/TkAgg to avoid Qt platform plugin errors.
    If a GUI backend cannot be initialized, the saved CSV and PNG are still
    available.
    """
    if not GUI_AVAILABLE:
        print("\n視窗顯示不可用：目前 Python 環境沒有可用的 Tkinter/TkAgg GUI backend。")
        return

    try:
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

        table_window = tk.Tk()
    except Exception as exc:
        print(f"\n視窗顯示不可用，已改為只儲存 CSV/PNG。原因：{exc}")
        return

    table_window.title(f"TDOA Results Table, spacing = {spacing_m:g} m")
    table_window.geometry("1180x420")

    table_frame = ttk.Frame(table_window, padding=8)
    table_frame.pack(fill=tk.BOTH, expand=True)

    columns = list(df.columns)
    tree = ttk.Treeview(table_frame, columns=columns, show="headings")

    for column in columns:
        tree.heading(column, text=column)
        width = max(95, min(170, len(column) * 9 + 28))
        tree.column(column, width=width, anchor=tk.CENTER, stretch=True)

    for _, row in df.iterrows():
        values = [_format_table_value(row[column]) for column in columns]
        tree.insert("", tk.END, values=values)

    y_scrollbar = ttk.Scrollbar(table_frame, orient=tk.VERTICAL, command=tree.yview)
    x_scrollbar = ttk.Scrollbar(table_frame, orient=tk.HORIZONTAL, command=tree.xview)
    tree.configure(yscrollcommand=y_scrollbar.set, xscrollcommand=x_scrollbar.set)

    tree.grid(row=0, column=0, sticky="nsew")
    y_scrollbar.grid(row=0, column=1, sticky="ns")
    x_scrollbar.grid(row=1, column=0, sticky="ew")
    table_frame.rowconfigure(0, weight=1)
    table_frame.columnconfigure(0, weight=1)

    chart_window = tk.Toplevel(table_window)
    chart_window.title(f"TDOA Error vs Distance, spacing = {spacing_m:g} m")
    chart_window.geometry("920x680")

    fig = create_error_figure(df, spacing_m)
    canvas = FigureCanvasTkAgg(fig, master=chart_window)
    canvas.draw()

    toolbar = NavigationToolbar2Tk(canvas, chart_window, pack_toolbar=False)
    toolbar.update()
    toolbar.pack(side=tk.BOTTOM, fill=tk.X)
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    table_window.mainloop()


def main():
    """Read user input, run the simulation, and save outputs."""
    user_input = input("請輸入探頭間距 spacing，單位 cm，例如 50 代表 50 cm: ").strip()
    spacing_cm = float(user_input)
    if spacing_cm <= 0:
        raise ValueError("spacing_cm must be positive.")

    spacing_m = spacing_cm / 100.0

    print("\nSimulation parameters:")
    print(f"sound_speed = {SOUND_SPEED_MPS} m/s")
    print(f"pulse_width = {PULSE_WIDTH_S:.3e} s")
    print(f"sampling_rate = {SAMPLING_RATE_HZ:.0f} Hz")
    print(f"timing_noise_std = {TIMING_NOISE_STD_S:.3e} s")
    print(f"monte_carlo_trials = {MONTE_CARLO_TRIALS}")
    print(f"spacing = {spacing_m:g} m")

    df = run_simulation(spacing_cm)

    script_dir = Path(__file__).resolve().parent
    csv_path = script_dir / CSV_FILENAME
    plot_path = script_dir / PLOT_FILENAME

    df.to_csv(csv_path, index=False)
    plot_results(df, spacing_m)

    print("\nCSV 儲存位置:")
    print(csv_path)
    print("\n圖片儲存位置:")
    print(plot_path)
    print("\n結果表格:")
    print(df.to_string(index=False))

    print("\n正在開啟圖表視窗與表格視窗...")
    show_results_windows(df, spacing_m)


if __name__ == "__main__":
    main()

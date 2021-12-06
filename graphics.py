import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from math import pi
import os


def plot_feature_importances():
    """ """
    pass


def polar_plot(data_set: pd.DataFrame, name=None, save_path=None):

    if name is None:
        name = name
    data = pd.DataFrame()

    # calculate normalized lipinski

    data["MolWt"] = [i / 500 for i in data_set["MolWt"]]
    data["LogP"] = [i / 5 for i in data_set["LogP"]]
    data["HBA"] = [i / 10 for i in data_set["NumHAcceptors"]]
    data["HBD"] = [i / 5 for i in data_set["NumHDonors"]]
    data["RotB"] = [i / 10 for i in data_set["NumRotatableBonds"]]
    data["TPSA"] = [i / 140 for i in data_set["TPSA"]]
    categories = list(
        data.columns
    )  # This will set up the parameters for the angles of the radar plot.
    N = len(categories)
    values = data[categories].values[0]
    values = np.append(values, values[:1])
    angles = [n / float(N) * 2 * pi for n in range(N)]
    angles += angles[:1]

    Ro5_up = [1, 1, 1, 1, 1, 1, 1]  # The upper limit for bRo5
    Ro5_low = [0.5, 0.1, 0, 0.25, 0.1, 0.5, 0.5]  # The lower limit for bRo5
    # fig=plt.figure()

    fig, ax = plt.subplots(subplot_kw={"projection": "polar"})

    plt.xticks(
        angles[:-1],
        categories,
        color="k",
        size=20,
        ha="center",
        va="top",
        fontweight="book",
    )

    plt.tick_params(axis="y", width=4, labelsize=12, grid_alpha=0.05)

    ax.set_rlabel_position(0)

    # ax.fill(angles, Ro5_up, 'red', alpha=0.2)
    ax.fill(angles, Ro5_low, "orangered", alpha=0.2)
    mean_vals = data.mean().to_list()
    mean_vals = np.append(mean_vals, mean_vals[:1])
    for i in data.index[:]:
        values = data[categories].values[i]
        values = np.append(values, values[:1])
        ax.plot(angles, values, linewidth=0.7, color="steelblue", alpha=0.5)
        # ax.fill(angles, values, 'C2', alpha=0.025)
    ax.plot(angles, mean_vals, linewidth=2, linestyle="--", color="red")
    ax.plot(angles, Ro5_up, linewidth=2, linestyle="-", color="black")
    ax.plot(angles, Ro5_low, linewidth=2, linestyle="-", color="black")

    ax.grid(axis="y", linewidth=1.5, linestyle="dotted", alpha=0.8)
    ax.grid(axis="x", linewidth=2, linestyle="-", alpha=1)
    save_fig(ax, name=name, save_path=save_path)


def save_fig(ax, name: str, save_path=None, dpi=300):
    """Frequently used snippet to save a plot

    Args:
        ax ([type]): Seaborn axis
        name (str): Name of the plot
        save_path ([type], optional): Where to save the plot. Defaults to None:str.
        dpi (int, optional): Resolution of the plot. Defaults to 300.

    Returns:
        (None)
    """
    if save_path is None:
        save_path = name
    else:
        save_path.joinpath(name)
    fig = ax.get_figure()
    fig.savefig(save_path, dpi=dpi)
    return None


if __name__ == "__main__":
    from processing import load_data

    data = pd.read_csv("converted_data.csv")
    biodeg = data[data["ReadyBiodegradability"] == 1]
    n_biodeg = data[data["ReadyBiodegradability"] == 0]
    polar_plot(biodeg, name="BiodegPolar.png")
    polar_plot(n_biodeg, name="NBiodegPolar.png")

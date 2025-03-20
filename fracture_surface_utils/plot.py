import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns

from fracture_surface_utils.data_processing import Data_Processing


class Plotter:

    def __init__(self, Result: Data_Processing):

        self.analysis = Result

        self.sliced_contour = Result.sliced_results
        self.shearlip_analysis = Result.shearlips_results
        self.piecewise_linear_fit = Result.pwlf_results
        self.modi_analysis = Result.modi_results

    def plot_fitted_contour(self):
        """Plotter - self explaining ."""

        analysis = self.analysis
        sliced_contour = self.sliced_contour
        piecewise_linear_fit = self.piecewise_linear_fit

        x_coordinate = analysis.x_coordinate
        specimen_name = analysis.specimen_name
        side = analysis.side
        output_path = analysis.output_path
        plotpath = os.path.join(output_path, "01_Contour_Plots")
        if not os.path.exists(plotpath):
            os.makedirs(plotpath)

        fig, axs = plt.subplots(
            1,
            1,
            gridspec_kw={"hspace": 0, "wspace": 0, "height_ratios": [1]},
            figsize=(4, 4),
        )

        sns.scatterplot(
            x=sliced_contour[x_coordinate]["Z_Coordinates_Raw"],
            y=sliced_contour[x_coordinate]["Y_Coordinates_Raw"],
            ax=axs,
            marker="^",
            label=f"Scanned Contour",
        )

        for segment in [
            item for item in piecewise_linear_fit[x_coordinate]["Segment_Index"].keys()
        ]:
            sns.lineplot(
                x=piecewise_linear_fit[x_coordinate]["Segment_Index"][segment][
                    "Z_Coordinates_Fit"
                ],
                y=piecewise_linear_fit[x_coordinate]["Segment_Index"][segment][
                    "Y_Coordinates_Fit"
                ],
                ax=axs,
                label="Linear Fit",
                color="k",
            )

        handles, labels = axs.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        axs.legend(
            by_label.values(),
            by_label.keys(),
            bbox_to_anchor=(0, 1.02, 1, 0.2),
            loc="lower left",
            mode="expand",
            borderaxespad=0,
            ncol=2,
        )

        axs.set_xlabel(r"$\it z$ [mm]")
        axs.set_ylabel(r"$\it y$ [mm]")

        for bp in piecewise_linear_fit[x_coordinate]["Breakpoints"]:
            axs.axvline(bp, color="grey", linestyle="--")

        axs.axis("equal")
        axs.set_xlim(-1.5, 1.5)
        plt.tight_layout()

        output_name = f"{x_coordinate}_{specimen_name}_{side}.png"
        save = os.path.join(plotpath, output_name)
        plt.savefig(save, dpi=300, bbox_inches="tight")
        plt.clf()
        plt.close()

    def plot_data_evaluation(self, legend_off: bool = True):
        """Plotter - self explaining .

        Parameters
        ----------

        legend_off : bool, default = True
                Turns the legend at plot 2,3 and 4 off.

        """

        analysis = self.analysis

        sliced_contour = self.sliced_contour
        shearlip_analysis = self.shearlip_analysis
        piecewise_linear_fit = self.piecewise_linear_fit
        modi_analysis = self.modi_analysis

        x_coordinate = analysis.x_coordinate
        specimen_name = analysis.specimen_name
        side = analysis.side
        output_path = analysis.output_path
        plotpath = os.path.join(output_path, "02_Data_Evaluation_Plots")
        if not os.path.exists(plotpath):
            os.makedirs(plotpath)

        sns.color_palette("colorblind")
        colorpalette = sns.color_palette("colorblind")
        color_map = [colorpalette[1], colorpalette[0], colorpalette[2], colorpalette[3]]

        fig, axs = plt.subplots(
            4,
            1,
            sharex="col",
            gridspec_kw={"hspace": 0, "wspace": 0, "height_ratios": [1, 0.6, 0.6, 0.6]},
            figsize=(4, 6),
        )

        sns.scatterplot(
            x=sliced_contour[x_coordinate]["Z_Coordinates_Contour"],
            y=sliced_contour[x_coordinate]["Y_Coordinates_Contour"],
            ax=axs[0],
            marker="^",
            label=f"Contour",
        )

        for segment in [
            item for item in piecewise_linear_fit[x_coordinate]["Segment_Index"].keys()
        ]:
            segment_idx = segment - 1

            sns.lineplot(
                x=piecewise_linear_fit[x_coordinate]["Segment_Index"][segment][
                    "Z_Coordinates_Fit"
                ],
                y=piecewise_linear_fit[x_coordinate]["Segment_Index"][segment][
                    "Y_Coordinates_Fit"
                ],
                ax=axs[0],
                label="Linear Fit",
                color="k",
            )

            array_lenght = len(
                piecewise_linear_fit[x_coordinate]["Segment_Index"][segment][
                    "Z_Coordinates_Fit"
                ]
            )
            sns.lineplot(
                x=piecewise_linear_fit[x_coordinate]["Segment_Index"][segment][
                    "Z_Coordinates_Fit"
                ],
                y=[
                    piecewise_linear_fit[x_coordinate]["Segment_Index"][segment][
                        "R_Squared"
                    ]
                ]
                * array_lenght,
                ax=axs[1],
                color="k",
                label="Correlation-\n coeffizient",
            )

            sns.lineplot(
                x=piecewise_linear_fit[x_coordinate]["Segment_Index"][segment][
                    "Z_Coordinates_Fit"
                ],
                y=[modi_analysis[x_coordinate]["Fracture_Mode_Classification_Category"]]
                * array_lenght,
                ax=axs[3],
                color=color_map[
                    modi_analysis[x_coordinate]["Fracture_Mode_Classification_Category"]
                ],
                label="Fracture Mode",
            )

            sns.lineplot(
                x=piecewise_linear_fit[x_coordinate]["Segment_Index"][segment][
                    "Z_Coordinates_Fit"
                ],
                y=[
                    shearlip_analysis[x_coordinate]["Shearlip_Classification_Category"][
                        segment_idx
                    ]
                ]
                * array_lenght,
                color="k",
                ax=axs[2],
                label="Shearlip",
            )

        axs[0].axis("equal")
        axs[0].set_xlim(-1.5, 1.5)

        axs[0].set_xlabel(r"$\it z$ [mm]")
        axs[0].set_ylabel(r"$\it y$ [mm]")

        axs[1].set_ylim(-0.15, 1.15)
        axs[1].axhspan(-0.15, analysis.corr_coeff_min, facecolor="red", alpha=0.2)
        axs[1].axhspan(analysis.corr_coeff_min, 1.15, facecolor="green", alpha=0.2)

        axs[1].set_ylabel(
            f"Correlation \n coefficient",
            rotation=90,
        )

        axs[2].set_ylim(-0.5, 2.5)
        axs[2].set_yticks([-1, 0, 1, 2, 3])
        axs[2].set_yticklabels(["", "Undef.", "No", "Yes", ""])
        axs[2].set_ylabel(
            f"Shearlip \n mode",
            rotation=90,
        )

        axs[3].set_ylim(-0.5, 3.5)
        axs[3].set_yticks(np.arange(0, 4, 1))
        axs[3].set_yticklabels(["Flat", "Slant", "V-Mode", "Transition"])
        axs[3].set_ylabel(
            f"Fracture \n mode",
            rotation=90,
        )
        axs[3].set_xlabel(r"$\it z$ [mm]")

        for i in range(4):
            handles, labels = axs[i].get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            if i == 0:
                axs[i].legend(
                    by_label.values(),
                    by_label.keys(),
                    bbox_to_anchor=(0, 1.02, 1, 0.2),
                    loc="lower left",
                    mode="expand",
                    borderaxespad=0,
                    ncol=1,
                )
            else:
                axs[i].legend(by_label.values(), by_label.keys(), loc="lower center")
                if np.any(legend_off):
                    axs[i].get_legend().remove()
            for bp in piecewise_linear_fit[x_coordinate]["Breakpoints"]:
                axs[i].axvline(bp, color="grey", linestyle="--")

        plt.tight_layout()
        output_name = f"{x_coordinate}_{specimen_name}_{side}.png"
        save = os.path.join(plotpath, output_name)
        plt.savefig(save, dpi=300, bbox_inches="tight")
        plt.clf()
        plt.close()
        print(f"Plotted coordinate x={x_coordinate} mm")

    def plot_xy_contour(
        self,
        contour_to_plot: str = "Contour",
        limits_x_plot: list = None,
        limits_angle_plot: list = None,
    ):
        """
        Plotter - self explaining.

        Parameters
        ----------

        plot_limits : bool, default = True
                Adjusts the plotted contour towards the analyzed limits of the sliced xy contour.

        contour_to_plot: str, default = 'Contour'
                Plot raw contour 'Raw' or interpolated contour 'Contour'.


        limits_x_plot: list = [float, float]
                Adjust x limits for plotting.

        limits_angle_plot: list = [float, float]
                Adjust angle limits for plotting.
        """

        analysis = self.analysis

        sliced_contour_xy = self.analysis.sliced_contour_xy

        contour_id = contour_to_plot

        specimen_name = analysis.specimen_name
        side = analysis.side
        output_path = analysis.output_path
        plotpath = os.path.join(output_path, "04_Crack_Path_Plots")
        if not os.path.exists(plotpath):
            os.makedirs(plotpath)

        fig, axs = plt.subplots(
            2,
            1,
            sharex="col",
            gridspec_kw={"hspace": 0, "wspace": 0, "height_ratios": [2, 1]},
            figsize=(8, 6),
        )

        ax1 = axs[0]
        ax2 = axs[1]

        ax1.grid(True, which="both")
        ax2.grid(True, which="both")
        sns.lineplot(
            x=sliced_contour_xy[self.analysis.z_coordinate_slice][
                f"X_Coordinates_{contour_id}"
            ],
            y=sliced_contour_xy[self.analysis.z_coordinate_slice][
                f"Y_Coordinates_{contour_id}"
            ],
            color="k",
            ax=ax1,
        )

        sns.lineplot(
            x=sliced_contour_xy[self.analysis.z_coordinate_slice][
                f"X_Coordinates_{contour_id}"
            ],
            y=sliced_contour_xy[self.analysis.z_coordinate_slice][
                f"Kinking_Angle_{contour_id}[°]"
            ],
            color="k",
            ax=ax2,
        )

        ax1.axis(
            "equal",
        )
        ax1.set_ylim(
            sliced_contour_xy[self.analysis.z_coordinate_slice][
                f"Y_Coordinates_{contour_id}"
            ].min()
            - 20,
            sliced_contour_xy[self.analysis.z_coordinate_slice][
                f"Y_Coordinates_{contour_id}"
            ].max()
            + 20,
        )

        if np.any(limits_x_plot):
            limits_x_plot = sorted([self.analysis.flip * els for els in limits_x_plot])
            ax1.set_xlim(min(limits_x_plot), max(limits_x_plot))
        if np.any(limits_angle_plot):
            ax2.set_ylim(min(limits_angle_plot), max(limits_angle_plot))

        axs[1].set_xlabel(r"$\it x$ [mm]")
        axs[0].set_ylabel(r"$\it y$ [mm]")
        axs[1].set_ylabel(r"$\phi_{0}$ [°]")

        plt.tight_layout()
        output_name = f"z_{self.analysis.z_coordinate_slice}_{specimen_name}_{side}.png"
        save = os.path.join(plotpath, output_name)
        plt.savefig(save, dpi=300, bbox_inches="tight")
        plt.clf()
        plt.close()
        print(f"Plotted slice z={self.analysis.z_coordinate_slice} mm")



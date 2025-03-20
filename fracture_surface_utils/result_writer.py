import os
import csv
import pandas as pd

from fracture_surface_utils.data_processing import Data_Processing


class Result_Writer:

    def __init__(self, Result: Data_Processing):

        self.analysis = Result

        self.sliced_contour = Result.sliced_results
        self.shearlip_analysis = Result.shearlips_results
        self.piecewise_linear_fit = Result.pwlf_results
        self.modi_analysis = Result.modi_results

    def write_to_csv(self):
        """
        Write analyzed data to .csv and analysis parameter to .txt .

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
        result_path = os.path.join(output_path, "03_Data_Evaluation")
        if not os.path.exists(result_path):
            os.makedirs(result_path)

        self.result_path = result_path

        filename = f"{x_coordinate}_{specimen_name}_{side}"

        user_input = {
            k: v for k, v in vars(analysis).items() if isinstance(v, (int, str, float))
        }
        user_input.pop("x_coordinate", None)

        with open(
            os.path.join(result_path, f"Parameter_{specimen_name}_{side}.txt"), "w"
        ) as txt_file:
            writer = csv.DictWriter(txt_file, user_input.keys())
            writer.writeheader()
            writer.writerow(user_input)

        writing_dict = {
            "X_Coordinate": x_coordinate,
            "Total_Width": piecewise_linear_fit[x_coordinate]["Total_Width"],
            "Total_Lenght": piecewise_linear_fit[x_coordinate]["Total_Lenght"],
            "Centerpoint": piecewise_linear_fit[x_coordinate]["Centerpoint"],
            "Central_Slope_Angle": piecewise_linear_fit[x_coordinate][
                "Central_Slope_Angle"
            ],
            "Fracture_Mode_Classification": modi_analysis[x_coordinate][
                "Fracture_Mode_Classification"
            ],
            "Number_of_Breakpoints": len(
                piecewise_linear_fit[x_coordinate]["Breakpoints"]
            ),
            "Number_of_Segments": analysis.number_of_segments,
        }

        for segment in [
            item for item in piecewise_linear_fit[x_coordinate]["Segment_Index"].keys()
        ]:
            breakpoint_idx = segment - 1

            writing_dict.update(
                {
                    f"Width_Segment_{segment}": piecewise_linear_fit[x_coordinate][
                        "Segment_Index"
                    ][segment]["Segment_Width"],
                    f"Lenght_Segment_{segment}": piecewise_linear_fit[x_coordinate][
                        "Segment_Index"
                    ][segment]["Segment_Lenght"],
                    f"Slope_Angle_{segment}": piecewise_linear_fit[x_coordinate][
                        "Segment_Index"
                    ][segment]["Slope_Angle"],
                    f"Shearlip_Classification_{segment}": str(
                        shearlip_analysis[x_coordinate]["Shearlip_Classification"][
                            segment - 1
                        ]
                    ),
                }
            )

        for breakpoint_idx in range(
            0, len(piecewise_linear_fit[x_coordinate]["Breakpoints"])
        ):
            writing_dict.update(
                {
                    f"Breakpoint_{breakpoint_idx}": piecewise_linear_fit[x_coordinate][
                        "Breakpoints"
                    ][breakpoint_idx]
                }
            )

        with open(os.path.join(result_path, f"{filename}.csv"), "w") as csv_file:
            writer = csv.DictWriter(csv_file, writing_dict.keys())
            writer.writeheader()
            writer.writerow(writing_dict)

        print(f"Wrote results for coordinate x={x_coordinate} mm")
        print(f"-------------------------------------------------")

    def write_to_csv_crackpath(self):
        """
        Write analyzed crack path coordinates to .csv

        """

        analysis = self.analysis
        sliced_contour_xy = self.analysis.sliced_contour_xy

        specimen_name = analysis.specimen_name
        side = analysis.side
        output_path = analysis.output_path
        result_path = os.path.join(output_path, "04_Crack_Path_Evaluation")
        if not os.path.exists(result_path):
            os.makedirs(result_path)

        self.result_path = result_path

        filename = f"z_{analysis.z_coordinate_slice}_{specimen_name}_{side}"

        writing_dict = {
            "X_Coordinate_Contour": sliced_contour_xy[analysis.z_coordinate_slice][
                "X_Coordinates_Contour"
            ],
            "Y_Coordinate_Contour": sliced_contour_xy[analysis.z_coordinate_slice][
                "Y_Coordinates_Contour"
            ],
            "Kinking_Angle_Contour[°]": sliced_contour_xy[analysis.z_coordinate_slice][
                "Kinking_Angle_Contour[°]"
            ],
        }

        df_contour = pd.DataFrame.from_dict(writing_dict)
        df_raw = pd.DataFrame(
            data={
                "X_Coordinate_Raw": sliced_contour_xy[analysis.z_coordinate_slice][
                    "X_Coordinates_Raw"
                ],
                "Y_Coordinate_Raw": sliced_contour_xy[analysis.z_coordinate_slice][
                    "Y_Coordinates_Raw"
                ],
            }
        ).sort_values(by="X_Coordinate_Raw")

        df_contour.to_csv(
            os.path.join(result_path, f"{filename}_contour.csv"), index=False
        )
        df_raw.to_csv(os.path.join(result_path, f"{filename}_raw.csv"), index=False)
        print(f"Wrote results for coordinate z={analysis.z_coordinate_slice} mm")
        print(f"-------------------------------------------------")

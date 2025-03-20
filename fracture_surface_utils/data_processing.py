import trimesh
import numpy as np
import pwlf
import os
from scipy.interpolate import interp1d
import math
from scipy.stats import linregress


class Data_Processing:
    def __init__(
        self,
        specimen_name: str = "not defined",
        side: str = None,
        data_path_input: str = None,
        filename: str = None,
        x_coordinate: float = 0,
        cutoff_z_perc: float = 0.1,
        cutoff_y_perc: float = 0.9,
        cutoff_y_abs: float = 2,
    ):
        """
        Parameter for slicing the processed .stl file for a given x coordinate or coordinate range.

        Parameters
        ----------
        specimen_name : str
                self-explaining
        side : str
                side, has to be "left" or "right"

        data_path_input : str or os.path.like
                datapath where the file is stored
        filename : str
                filename of .stl file to process

        x_coordinate: float
                slicing plane x coordinate
        cutoff_z_perc : int
                percentage of masked area within xy plane regarding the contour width
        cutoff_y_perc : int
                percentage of masked area within xy plane regarding the contour height
        cutoff_y_abs : int
                absolute value of masked area below maximum y value of the contour. Overrides the evaluation
                based on a given percentange



        Sometimes there are scanning artefacts within the .stl file that are not part of the fractured surface
        -> therefore we need to mask the contour.
        To ensure that we only analyze the fractured surface, a certain percentage of area at the outer edge of the
        contour within the xz-plane is masked.

        """

        global_path = os.getcwd()
        self.filename = filename
        self.data_path_input = data_path_input
        self.data_input = os.path.join(data_path_input, filename)

        self.specimen_name = specimen_name
        self.side = side
        if self.side == "left":
            self.flip = -1
        else:
            self.flip = 1

        self.output_path = os.path.join(
            global_path, "02_results", self.specimen_name, self.side
        )

        if not os.path.exists(self.output_path):
            os.makedirs(self.output_path)

        self.filename_stl = f"{specimen_name}_finemesh.stl"
        self.x_coordinate = x_coordinate

        self.cutoff_y_perc = cutoff_y_perc
        self.cutoff_z_perc = cutoff_z_perc
        self.cutoff_y_abs = cutoff_y_abs

        self.mesh = trimesh.load_mesh(self.data_input)

        self.number_of_segments = None

        self.segment_percentage_min = None
        self.segment_percentage_max = None
        self.corr_coeff_min = None
        self.angle_difference_central_min = None

        self.undefined_percentage_max = None
        self.shearlip_percentage_max = None
        self.reference_angle_max = None
        self.reference_angle_min = None

        self.sliced_results = None
        self.pwlf_results = None
        self.shearlips_results = None
        self.modi_results = None

    def slicer_yz_plane(self, interpolation_value: float, reposition: bool = False):
        """
        Function for slicing the .stl file for a given x coordinate at the zy plane.

        Parameters
        ----------

        interpolation_value : float
                number of datapoints to interpolate the given contour to
        repositon : bool, can be True or False
                if the specimen bends, the contour slice gets repositioned so the contour is centered and symmetrical
                around zero


        Returns:
        ----------
        x_slice_to_yz_plane: dict
                Dictionary with sliced x_coordinate as key and corresponding yz-plane coordinates as values.
                This includes both masked and interpolated ('_Contour') as well as raw ('_Raw') data.
        """

        # Read and process mesh

        x_slice_to_yz_plane = {}
        reposition = reposition

        x_coordinate = self.x_coordinate

        slice_mesh = self.mesh.section(
            plane_origin=[x_coordinate, 0, 0], plane_normal=[1, 0, 0]
        )
        slice_2D, to_3D = slice_mesh.to_planar()

        vertices = np.asarray(slice_2D.vertices)
        z_coordinates_raw = z_coordinates = vertices[:, 0]
        y_coordinates_raw = y_coordinates = vertices[:, 1]
        z_min, z_max = vertices[:, 0].min(), vertices[:, 0].max()
        y_max = vertices[:, 1].max()

        if reposition == True:

            self.thickness = abs(z_min) + abs(z_max)

            ###normalize
            col_y = vertices[:, 1].tolist()

            start = -1
            width = self.thickness
            end = start + width
            col_z_norm = (
                (vertices[:, 0] - vertices[:, 0].min())
                / (vertices[:, 0].max() - vertices[:, 0].min())
                * width
                + start
            ).tolist()

            vertices = np.stack((col_z_norm, col_y), axis=-1)

            z_coordinates_raw = vertices[:, 0]
            y_coordinates_raw = vertices[:, 1]

        # Mask data according to parameters defined in Data Processing

        z_lim_lower, z_lim_upper = (
            z_min - self.cutoff_z_perc * z_min,
            z_max - self.cutoff_z_perc * z_max,
        )

        if self.cutoff_y_abs is not None:
            y_lim = y_max - self.cutoff_y_abs
        else:
            y_lim = y_max - self.cutoff_y_perc * y_max

        mask_z = np.all(
            [z_coordinates >= z_lim_lower, z_coordinates <= z_lim_upper], axis=0
        )
        mask_y = np.all([y_coordinates > y_lim], axis=0)

        masked_vertices = vertices[mask_y & mask_z]
        masked_y_coordinates, masked_z_coordinates = (
            masked_vertices[:, 1],
            masked_vertices[:, 0],
        )

        # linear interpolation of unregular data

        self.z_int = np.linspace(
            masked_z_coordinates.min(),
            masked_z_coordinates.max(),
            num=interpolation_value,
            endpoint=True,
        )
        interpolate = interp1d(masked_z_coordinates, masked_y_coordinates)
        self.y_int = interpolate(self.z_int)

        z_int = self.z_int
        y_int = self.y_int

        x_slice_to_yz_plane[x_coordinate] = {
            "Z_Coordinates_Contour": z_int,
            "Y_Coordinates_Contour": y_int,
            "Z_Coordinates_Raw": z_coordinates_raw,
            "Y_Coordinates_Raw": y_coordinates_raw,
        }

        self.sliced_results = x_slice_to_yz_plane
        return x_slice_to_yz_plane

    def slicer_xy_plane(
        self,
        z_value: float = 0,
        side: str = None,
        limits_x: list = [10, 70],
        interpolation_value: float = 30,
        normal_vector: list = None,
    ):
        """Slicing along the xy plane to get the crack path.

        Parameters
        ----------
        z_value : float
                Position where place is sliced. Default = 0 equals specimen center.
        side : str
                side, has to be "left" or "right"
        limits_x: list, [int, int]
                absolute value of maximum and minimum x coordinate
        interpolation_value : float
                number of datapoints to interpolate the given contour to
        normal_vector : list, [float, float, float]
                manually overriding the preset normal vector, especially for bended specimen

        Returns:
        ----------
        z_slice_to_xy_plane: dict
                Dictionary with sliced z_coordinate as key and corresponding xy-plane coordinates as values.
                This includes both masked and interpolated ('_Contour') as well as raw ('_Raw') data as well as the
                kinking angle with regards to the interpolated contour..
        """

        z_slice_to_xy_plane = {}

        self.z_coordinate_slice = z_value
        self.limits_x = limits_x[0], limits_x[1]

        # normal vector depends on the sign of the z value

        # value_to_oondition = {0: [0, 0, 1], -1: [0, -0.0105, 0.99], 1:[0, 0.0153, -0.99]}

        value_to_oondition = {0: [0, 0, 1], -1: [0, -0.0105, 0.99], 1: [0, 0.0153, -1]}

        if self.z_coordinate_slice == 0:
            self.normal_vector = value_to_oondition.get(self.z_coordinate_slice)

        else:
            self.normal_vector = value_to_oondition.get(
                np.sign(self.z_coordinate_slice)
            )

        # overriding normal vector
        if np.any(normal_vector):
            self.normal_vector = normal_vector

        slice_mesh_xy = self.mesh.section(
            plane_origin=[0, 0, self.z_coordinate_slice],
            plane_normal=self.normal_vector,
        )
        slice_2D, to_3D = slice_mesh_xy.to_planar()

        vertices = np.asarray(slice_2D.vertices)
        x_coordinates_raw = x_coordinates = vertices[:, 0]
        y_coordinates_raw = y_coordinates = vertices[:, 1]

        # flip limits if testside == left
        self.limits_x = sorted([self.flip * els for els in self.limits_x])

        # mask according to x limits
        mask_x = np.all(
            [x_coordinates >= self.limits_x[0], x_coordinates <= self.limits_x[1]],
            axis=0,
        )
        masked_vertices = vertices[mask_x]
        masked_y_coordinates, masked_x_coordinates = (
            masked_vertices[:, 1],
            masked_vertices[:, 0],
        )

        # linear interpolation of unregular data
        self.x_int = np.linspace(
            masked_x_coordinates.min(),
            masked_x_coordinates.max(),
            num=interpolation_value,
            endpoint=True,
        )
        interpolate = interp1d(masked_x_coordinates, masked_y_coordinates)
        self.y_int = interpolate(self.x_int)

        x_int = self.x_int
        y_int = self.y_int

        # calculate kinking angle between the coordinate points. prepend forces a division by zero warning
        np.seterr(divide="ignore", invalid="ignore")

        kinking_angle_contour = np.degrees(
            np.arctan(
                np.diff(y_int, prepend=y_int[0]) / np.diff(x_int, prepend=x_int[0])
            )
        )
        kinking_angle_raw = np.degrees(
            np.arctan(
                np.diff(y_coordinates, prepend=y_coordinates[0])
                / np.diff(x_coordinates, prepend=x_coordinates[0])
            )
        )

        z_slice_to_xy_plane[self.z_coordinate_slice] = {
            "X_Coordinates_Contour": x_int,
            "Y_Coordinates_Contour": y_int,
            "Kinking_Angle_Contour[°]": kinking_angle_contour,
            "X_Coordinates_Raw": x_coordinates_raw,
            "Y_Coordinates_Raw": y_coordinates_raw,
            "Kinking_Angle_Raw[°]": kinking_angle_raw,
        }

        self.sliced_contour_xy = z_slice_to_xy_plane

        return z_slice_to_xy_plane

    def linear_fit_piecewise(self, number_of_segments: int = 3):
        """Fitting continuous piecewise linear functions to given fractured surface contour .

        Parameters
        ----------
        number_of_segments : int
                desired number of segments (did not test for values other than 3 )

        Returns :
        ----------
        dict_pwlf : dict
                nested dictionary for each sliced x coordinate containing information about the fitted segments.
                This includes:

        'Breakpoints' : ndarray
                z-coordinate of breakpoints between the segments. First and last item define the outermost
                points of the contour
        'Total_Width' : float
                Total width of the evaluated and fitted contour. Equals projected distance between the outermost breakpoint z coordinates
        'Total_Lenght' : float
                Total lenght of the evaluated and fitted contour. Equals sum of all euclidian distances between the breakpints
        'Centerpoint' : float
                z coordinate at centerpoint of contour t/2
        'Central_Slope_Angle' : float
                Slope angle of the segment at centerpoint t/2 in degrees[°]

        'Segment_Index' : dict
                For each segment the dictionary outlines specific properties.
                This includes:

                'Breakpoint_Z' : tuple (float, ..., float), len = number_of_segments+1
                    Minimum and maximum z-coordinates of the linear segment
                'Breakpoint_Y' : tuple (float, ..., float), len = number_of_segments+1
                    Minimum and maximum y-coordinates of the linear segment
                'Segment_Lenght' : float
                    Euclidian distance between the breakpoints of the segment
                'Segment_Width' : float
                    Projected distance between z-coordinates of breakpoints
                'Segment_Width_Percentage' : float
                    Percentage of segment width with regards to the total initial contour width in %
                'Slope_Angle' : float
                    Linear slope angle in degrees [°]
                'R_Squared' : float
                    Correlation coefficient of linear fit R**2
                'Angle_Difference_Central' : float
                    Absolute difference of angles of the respective segment towards the angle at t/2 at the specimen
                    center in degrees [°]
                'Z_Coordinates_Fit' : ndarray
                    Z coordinates of fitted linear segment
                'Y_Coordinates_Fit' : ndarray
                    Y coordinates of fitted linear segment


        """

        dict_pwlf = {}

        self.number_of_segments = number_of_segments

        z_int = self.z_int
        y_int = self.y_int
        x_coordinate = self.x_coordinate

        # piecewise linear fit

        zy_pwlf = pwlf.PiecewiseLinFit(z_int, y_int)
        breaks = zy_pwlf.fit(number_of_segments)

        # find segment index of segment that is centered (z=t/2) and its corresponding central slope angle

        centerpoint = sum(z_int) / float(len(z_int))
        central_segment_index = next(
            j for j, e in enumerate(breaks) if e >= centerpoint
        )
        slope_angle_center = (
            180.0 * np.arctan(zy_pwlf.slopes[central_segment_index - 1]) / np.pi
        )

        # calculate coordinates, lenght and correlation coefficient of each linear segment and save in dict

        segment_dict = {}

        lenseg_all = []

        for j, bp in enumerate(breaks[:-1]):
            dict_breaks = {}

            z_lin = np.linspace(breaks[j], breaks[j + 1], num=30, endpoint=True)
            y_lin = zy_pwlf.intercepts[j] + zy_pwlf.slopes[j] * z_lin

            idx_j = (np.abs(z_lin - breaks[j])).argmin()
            pz_j, py_j = z_lin[idx_j], y_lin[idx_j]

            idx_k = (np.abs(z_lin - breaks[j + 1])).argmin()
            pz_k, py_k = z_lin[idx_k], y_lin[idx_k]

            len_segment = math.dist([pz_j, py_j], [pz_k, py_k])
            lenseg_all.append(len_segment)
            width_segment = math.dist([pz_j], [pz_k])

            percentage_of_t0 = (
                width_segment / (breaks[number_of_segments] - breaks[0])
            ) * 100

            slope_angle = 180.0 * np.arctan(zy_pwlf.slopes[j]) / np.pi

            mask = np.where((z_int < pz_k) & (z_int > pz_j))
            z_seg, y_seg = z_int[mask], y_int[mask]

            # if the breakpoints are too close and the interpolation value too low,
            # linear regression cannot be calculated as segment coordinates array is empty
            try:
                res = linregress(z_seg, y_seg)
                r_squared = (res.rvalue) ** 2
            except ValueError:
                res = 0
                r_squared = 0

            angle_difference = abs(slope_angle) - abs(slope_angle_center)

            dict_breaks[j + 1] = {
                "Breakpoint_Z": (pz_j, pz_k),
                "Breakpoint_Y": (py_j, py_k),
                "Segment_Lenght": len_segment,
                "Segment_Width": width_segment,
                "Segment_Width_Percentage": percentage_of_t0,
                "Slope_Angle": slope_angle,
                "R_Squared": r_squared,
                "Angle_Difference_Central": angle_difference,
                "Z_Coordinates_Fit": z_lin,
                "Y_Coordinates_Fit": y_lin,
            }
            segment_dict.update(dict_breaks)

        dict_pwlf[x_coordinate] = {"Breakpoints": breaks}
        dict_pwlf[x_coordinate].update(
            {
                "Total_Width": breaks[number_of_segments] - breaks[0],
                "Total_Lenght": sum(lenseg_all),
                "Centerpoint": centerpoint,
                "Central_Slope_Angle": slope_angle_center,
            }
        )
        dict_pwlf[x_coordinate].update({"Segment_Index": segment_dict})

        self.pwlf_results = dict_pwlf

        return dict_pwlf

    def shearlip_condition(
        self,
        segment_percentage_min: float,
        segment_percentage_max: float,
        corr_coeff_min: float,
        angle_difference_central_min: float,
    ):
        """Analysis of the occurence of shearlips based on predefined criteria .
        See also:  https://doi.org/10.1111/j.1460-2695.2004.00837.x for further definition.

        Parameters
        ----------
        segment_percentage_min : float
                allowed percentage of shear lip width (t_s1 + t_s2) relative to the total specimen thickness t in percent
        segment_percentage_max : float
                allowed percentage of shear lip width (t_s1 + t_s2) relative to the total specimen thickness t in percent
        corr_coeff_min : float
                lowest limit of accepted correlation coefficient for linear fit
        angle_difference_central_min : lower limit of angle difference of central segment towards the outer
                segments in degrees [°]

        Returns :
        ----------
        shearlip_dict : dict
                nested dictionary for each sliced x coordinate containing information about the occurence of shearlips
                within the segments.
                This includes:

                'Shearlip_Classification' : tuple (str, ..., str), len = number_of_segments
                    Classification for each segment
                'Shearlip_Classification' : tuple (int, ..., int), len = number_of_segments
                    Classification for each segment as categorical value.
                    0 = Undefined,
                    1 = No shearlip,
                    2 = Shearlip

                'Segment_Width' : tuple (float, ..., float), len = number_of_segments
                    Projected distance between z-coordinates of breakpoints
                'Segment_Width_Percentage' : tuple (float, ..., float), len = number_of_segments

                The classification is made in three categories : Undefined, No shearlip, shearlip.
                To facilitate the evaluation, a dictionary contains the cumulative width and percentage of
                cumulative width of the occurence of the categories within the total width.

                Example - applied for each category given above.
                --------

                'Undefined' : dict
                    Dictionary containing information about the total share of the classified width
                    This contains :
                        'Sum_Segment_Width' : float
                            Total sum of all segment widths classified as 'undefined'
                         'Sum_Segment_Percentage' : float
                            Percentage share of segments classified as 'Undefined' in percent

        """

        x_coordinate = self.x_coordinate
        self.segment_percentage_min = segment_percentage_min
        self.segment_percentage_max = segment_percentage_max
        self.corr_coeff_min = corr_coeff_min
        self.angle_difference_central_min = angle_difference_central_min

        pwlf_dict = self.pwlf_results

        shearlip_dict = {}
        shearlip_dict[x_coordinate] = {}

        shearlip_to_category = {"Undefined": 0, "No Shearlip": 1, "Shearlip": 2}

        shearlip_list = []
        shearlip_categorical_list = []
        segment_width = []
        segment_percentage = []

        # setting the condition for case distinction

        for segment in [
            item for item in pwlf_dict[x_coordinate]["Segment_Index"].keys()
        ]:

            if segment == 2:

                shearlip_mode = "No Shearlip"
                shearlip_list.append(shearlip_mode)
                segment_percentage.append(
                    pwlf_dict[x_coordinate]["Segment_Index"][segment][
                        "Segment_Width_Percentage"
                    ]
                )
                segment_width.append(
                    pwlf_dict[x_coordinate]["Segment_Index"][segment]["Segment_Width"]
                )
                shearlip_categorical_list.append(
                    shearlip_to_category.get(str(shearlip_mode))
                )

            else:
                conditions_shearlips = [
                    (
                        pwlf_dict[x_coordinate]["Segment_Index"][segment]["R_Squared"]
                        < corr_coeff_min
                    ),
                    (
                        segment_percentage_min
                        > pwlf_dict[x_coordinate]["Segment_Index"][segment][
                            "Segment_Width_Percentage"
                        ]
                    )
                    | (
                        pwlf_dict[x_coordinate]["Segment_Index"][segment][
                            "Segment_Width_Percentage"
                        ]
                        > segment_percentage_max
                    ),
                    (
                        pwlf_dict[x_coordinate]["Segment_Index"][segment][
                            "Angle_Difference_Central"
                        ]
                        > angle_difference_central_min
                    ),
                ]

                choices_shearlips = ["Undefined", "No Shearlip", "Shearlip"]
                shearlip_mode = np.select(
                    conditions_shearlips, choices_shearlips, default="No Shearlip"
                )
                shearlip_list.append(shearlip_mode)
                shearlip_categorical_list.append(
                    shearlip_to_category.get(str(shearlip_mode))
                )
                segment_percentage.append(
                    pwlf_dict[x_coordinate]["Segment_Index"][segment][
                        "Segment_Width_Percentage"
                    ]
                )
                segment_width.append(
                    pwlf_dict[x_coordinate]["Segment_Index"][segment]["Segment_Width"]
                )

        shearlip_dict[x_coordinate] = {
            "Shearlip_Classification": tuple(shearlip_list),
            "Shearlip_Classification_Category": tuple(shearlip_categorical_list),
            "Segment_Width": tuple(segment_width),
            "Segment_Width_Percentage": tuple(segment_percentage),
        }

        for item in choices_shearlips:
            indices = [
                i
                for i, val in enumerate(
                    shearlip_dict[x_coordinate]["Shearlip_Classification"]
                )
                if val == item
            ]

            total_width_percentage = sum(
                [
                    shearlip_dict[x_coordinate]["Segment_Width_Percentage"][i]
                    for i in indices
                ]
            )
            total_width = sum(
                [shearlip_dict[x_coordinate]["Segment_Width"][i] for i in indices]
            )

            shearlip_dict[x_coordinate].update(
                {
                    item: {
                        "Sum_Segment_Width": total_width,
                        "Sum_Segment_Percentage": total_width_percentage,
                    }
                }
            )

        self.shearlips_results = shearlip_dict

        return shearlip_dict

    def modi_condition(
        self,
        undefined_percentage_max: float,
        shearlip_percentage_max: float,
        reference_angle_max: float,
        reference_angle_min: float,
    ):
        """Analysis of fracture mode based on predefined criteria .

        Parameters
        ----------
        undefined_percentage_max : float
                allowed total percentage of sum of segments being classified as 'undefined' relative
                to the total specimen width in percent
        shearlip_percentage_max : float
                allowed total percentage of sum of segments being classified as 'shearlip' relative
                to the total specimen width in percent
        reference_angle_max : float
                allowed central reference angle mainly for classifying slant and flat mode in degree
        reference_angle_min : float
                allowed central reference angle mainly for classifying slant and flat mode in degree


        Returns :
        ----------
        modi_dict : dict
                dictionary for each sliced x coordinate containing information about the fracture mode
                This includes:

                'Fracture_Mode_Classification' : str
                    Fracture mode
                'Fracture_Mode_Classification_Category' : int
                    Fracture mode as categorical value
                    0 = Flat
                    1 = Slant
                    2 = V-Mode
                    3 = Undefined or Transition

        """

        modi_dict = {}
        x_coordinate = self.x_coordinate

        self.undefined_percentage_max = undefined_percentage_max
        self.shearlip_percentage_max = shearlip_percentage_max
        self.reference_angle_max = reference_angle_max
        self.reference_angle_min = reference_angle_min

        pwlf_dict = self.pwlf_results
        shearlip_dict = self.shearlips_results

        mode_to_category = {
            "Undefined": 3,
            "Flat": 0,
            "V-Mode": 2,
            "Slant": 1,
            "Transition": 3,
        }

        # for v-mode, shearlips have to be at both sides of the specimen - count_of_shearlips = 2

        count_of_shearlips = shearlip_dict[x_coordinate][
            "Shearlip_Classification"
        ].count("Shearlip")

        sum_undefined_percentage = shearlip_dict[x_coordinate]["Undefined"][
            "Sum_Segment_Percentage"
        ]
        sum_shearlip_percentage = shearlip_dict[x_coordinate]["Shearlip"][
            "Sum_Segment_Percentage"
        ]

        # setting the condition for case distinction

        conditions_modi = [
            (
                (sum_undefined_percentage >= undefined_percentage_max)
                & (
                    abs(pwlf_dict[x_coordinate]["Central_Slope_Angle"])
                    > reference_angle_min
                )
            ),
            (
                (sum_shearlip_percentage <= shearlip_percentage_max)
                & (
                    abs(pwlf_dict[x_coordinate]["Central_Slope_Angle"])
                    > reference_angle_min
                )
            ),
            (
                (sum_shearlip_percentage >= shearlip_percentage_max)
                & (
                    abs(pwlf_dict[x_coordinate]["Central_Slope_Angle"])
                    > reference_angle_min
                )
            ),
            (
                (sum_shearlip_percentage <= shearlip_percentage_max)
                & (
                    abs(pwlf_dict[x_coordinate]["Central_Slope_Angle"])
                    < reference_angle_min
                )
            ),
            (
                (count_of_shearlips == 2)
                & (sum_shearlip_percentage >= shearlip_percentage_max)
                & (
                    abs(pwlf_dict[x_coordinate]["Central_Slope_Angle"])
                    < reference_angle_max
                )
            ),
        ]

        choices_modi = ["Undefined", "Slant", "Transition", "Flat", "V-Mode"]
        fracture_mode = np.select(conditions_modi, choices_modi, default="Transition")
        fracture_mode_category = mode_to_category.get(str(fracture_mode))
        modi_dict[x_coordinate] = {
            "Fracture_Mode_Classification": str(fracture_mode),
            "Fracture_Mode_Classification_Category": fracture_mode_category,
        }

        self.modi_results = modi_dict

        print(
            f"Analyzed {self.specimen_name}_{self.side} coordinate x={x_coordinate} mm"
        )

        return modi_dict
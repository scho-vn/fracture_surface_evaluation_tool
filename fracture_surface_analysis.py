import os
from fracture_surface_utils.data_processing import Data_Processing
from fracture_surface_utils.plot import Plotter
from fracture_surface_utils.result_writer import Result_Writer
from fracture_surface_utils.functions import sum_results
from fracture_surface_utils.functions import range_generator

global_path = os.getcwd()

specimen_name = "test_mesh"
side = "right"


##slicing  yz plane

analysis = Data_Processing(
    specimen_name=specimen_name,
    side=side,
    data_path_input=os.path.join(global_path, "01_raw_stl"),
    filename=f"{specimen_name}.stl",
    x_coordinate=None,
    cutoff_y_perc=None,
    cutoff_z_perc=0.15,
    cutoff_y_abs=2,
)

sliced_contour_xy = analysis.slicer_xy_plane(
    z_value=0.85,
    side=side,
    limits_x=[20, 160],
    interpolation_value=3000,
    normal_vector=[0, 0.0053, -1],
)
plotter = Plotter(Result=analysis)

plotter.plot_xy_contour(
    contour_to_plot="Contour", limits_x_plot=[10, 170], limits_angle_plot=[-70, 70]
)
Result_Writer(Result=analysis).write_to_csv_crackpath()


# generate x range for slicing zy plane

x_range = range_generator(side=side, limits=[20, 160], interval=40)




for element in x_range:

    analysis.x_coordinate = element

    sliced_contour = analysis.slicer_yz_plane(interpolation_value=30, reposition=True)
    piecewise_linear_fit = analysis.linear_fit_piecewise(number_of_segments=3)

    shearlip_analysis = analysis.shearlip_condition(
        segment_percentage_min=5,
        segment_percentage_max=40,
        corr_coeff_min=0.75,
        angle_difference_central_min=15,
    )

    modi_analysis = analysis.modi_condition(
        undefined_percentage_max=50,
        shearlip_percentage_max=30,
        reference_angle_max=15,
        reference_angle_min=5,
    )

    # plot results
    plotter = Plotter(Result=analysis)
    plotter.plot_fitted_contour()
    plotter.plot_data_evaluation(legend_off=True)

    # write results
    Result_Writer(Result=analysis).write_to_csv()


# sum it up
summary = sum_results(
    specimen_name=specimen_name,
    side=side,
    result_path=os.path.join(
        global_path, "02_results", f"{specimen_name}", f"{side}", "03_Data_Evaluation"
    ),
    delete_single_files=True,
)

print("done")
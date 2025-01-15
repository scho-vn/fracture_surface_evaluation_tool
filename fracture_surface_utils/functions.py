import os
import csv
import pandas as pd
import numpy as np

def sum_results(specimen_name: str = "not defined", side: str = None, result_path: str = None,
                delete_single_files: bool = True):
    """Loop through single files and sum it up .

    Parameters
    ----------
    specimen_name : str
                self-explaining
        side : str
                side, has to be "left" or "right"

    result_path: str
            Path containing the single .csv files
    delete_single_files: bool, default = True
            Deletes single .csv files after summary - cleanup function :).

    """
    result_path = result_path
    specimen_name = specimen_name
    side = side

    csv_files = [ x for x in list(filter(lambda f: f.endswith('.csv'), os.listdir(result_path))) if "Summary" not in x ]

    df_concat = pd.concat([pd.read_csv(os.path.join(result_path,f)) for f in sorted(csv_files) ], ignore_index=False).set_index('X_Coordinate', drop=False)

    if np.any(delete_single_files):
        for filename in sorted(csv_files):
            file = os.path.join(result_path, filename)
            os.remove(file)
    df_concat.to_csv(os.path.join(result_path, f'{specimen_name}_{side}_Summary.csv'))



def range_generator(side: str = None, limits: list = [10,70], interval: float = 1):

    """Generate the range of x coordinates to be sliced .

    Parameters
    ----------
    side : str
            side, has to be "left" or "right"

    limits: list, [int, int]
            absolute value of maximum and minimum x coordinate
    interval : float
            Interval between slicing planes

    """
    interval = interval
    upper_limit, lower_limit = limits
    if side == 'left':
        limits = [-1 * els for els in limits]

    upper_limit, lower_limit=max(limits), min(limits)

    range = [round(elem, 2) for elem in np.arange(lower_limit,upper_limit,interval)]

    return range




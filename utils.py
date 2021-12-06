import pandas as pd


def save_csv(data_frame, filename):
    """[summary]

    Args:
        data_frame ([type]): [description]
        filename ([type]): [description]
    """

    data_frame.to_csv(filename)

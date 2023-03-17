from math import sin, cos, pi
import numpy as np
from numpy import ndarray

__all__ = ["calculate_polarizer_matrix",
           "calculate_retarder_matrix"]


def calculate_polarizer_matrix(theta):
    """Calculate the polarizer matrix.

    Parameters
    ----------    
    theta: float
        Angle of the transmission axis of the polarizer in degrees

    """
    theta = np.deg2rad(theta)
    POLARIZER_MATRIX = 0.5 * np.asarray(
        [
            [1, cos(2 * theta), sin(2 * theta), 0],
            [
                cos(2 * theta),
                cos(2 * theta) ** 2,
                cos(2 * theta) * sin(2 * theta),
                0,
            ],
            [
                sin(2 * theta),
                cos(2 * theta) * sin(2 * theta),
                sin(2 * theta) ** 2,
                0,
            ],
            [0, 0, 0, 0],
        ]
    )
    return POLARIZER_MATRIX


def calculate_retarder_matrix(phase_difference, theta) -> ndarray:
    """Calculate the retarder matrix for a given phase difference.

    Parameters
    ----------
    phase_difference : float
        Phase difference in degrees.

    theta: float
        Angle of the transmission axis of the retarder in degrees.

    """
    phase_difference = np.deg2rad(phase_difference)
    theta = np.deg2rad(theta)
    retarder_matrix = np.asarray(
        [
            [1, 0, 0, 0],
            [
                0,
                cos(2 * theta) ** 2
                + sin(2 * theta) ** 2 * cos(phase_difference),
                cos(2 * theta)
                * sin(2 * theta)
                * (1 - cos(phase_difference)),
                -sin(2 * theta) * sin(phase_difference),
            ],
            [
                0,
                cos(2 * theta)
                + sin(2 * theta) * (1 - cos(phase_difference)),
                sin(2 * theta) ** 2
                + cos(2 * theta) ** 2 * cos(phase_difference),
                cos(2 * theta) * sin(phase_difference),
            ],
            [
                0,
                sin(2 * theta) * sin(phase_difference),
                -cos(2 * theta) * sin(phase_difference),
                cos(phase_difference),
            ],
        ]
    )

    return retarder_matrix

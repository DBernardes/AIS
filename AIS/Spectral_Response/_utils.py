import cmath
from math import cos, pi, sin

import numpy as np
from numpy import ndarray

__all__ = [
    "calculate_polarizer_matrix",
    "calculate_retarder_matrix",
    "apply_matrix",
    "calculate_depolarizer_matrix",
]


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
            [1, cos(2 * theta), 0, 0],
            [
                cos(2 * theta),
                1,
                0,
                0,
            ],
            [
                0,
                0,
                sin(2 * theta),
                0,
            ],
            [0, 0, 0, sin(2 * theta)],
        ],
        dtype=np.float64,
    )
    return POLARIZER_MATRIX


# def calculate_polarizer_matrix_1(theta):
#     """Calculate the polarizer matrix.

#     Parameters
#     ----------
#     theta: float
#         Angle of the transmission axis of the polarizer in degrees

#     """
#     theta = np.deg2rad(theta)
#     POLARIZER_MATRIX = 0.5 * np.asarray(
#         [
#             [1, cos(2 * theta), sin(2 * theta), 0],
#             [
#                 cos(2 * theta),
#                 cos(2 * theta) ** 2,
#                 cos(2 * theta) * sin(2 * theta),
#                 0,
#             ],
#             [
#                 sin(2 * theta),
#                 cos(2 * theta) * sin(2 * theta),
#                 sin(2 * theta) ** 2,
#                 0,
#             ],
#             [0, 0, 0, 0],
#         ]
#     )
#     return POLARIZER_MATRIX


# def calculate_retarder_matrix_1(phase_difference, theta) -> ndarray:
#     """Calculate the retarder matrix for a given phase difference.

#     Parameters
#     ----------
#     phase_difference : float
#         Phase difference in degrees.

#     theta: float
#         Angle of the transmission axis of the retarder in degrees.

#     """
#     # ! This equation is wrong
#     phase_difference = np.deg2rad(phase_difference)
#     theta = np.deg2rad(theta)
#     G = (1 + cos(phase_difference)) / 2
#     H = (1 - cos(phase_difference)) / 2
#     retarder_matrix = 0.5 * np.asarray(
#         [
#             [1, 0, 0, 0],
#             [
#                 0,
#                 G + H * cos(4 * theta),
#                 H * sin(4 * theta),
#                 -sin(phase_difference) * sin(2 * theta),
#             ],
#             [
#                 0,
#                 H * sin(4 * theta),
#                 G - H * cos(4 * theta),
#                 sin(phase_difference) * cos(2 * theta),
#             ],
#             [
#                 0,
#                 sin(phase_difference) * sin(2 * theta),
#                 -sin(phase_difference) * cos(2 * theta),
#                 cos(phase_difference),
#             ],
#         ],
#         dtype=np.float64,
#     )

#     return retarder_matrix


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
                cos(2 * theta) ** 2 + sin(2 * theta) ** 2 * cos(phase_difference),
                cos(2 * theta) * sin(2 * theta) * (1 - cos(phase_difference)),
                -sin(2 * theta) * sin(phase_difference),
            ],
            [
                0,
                cos(2 * theta) * sin(2 * theta) * (1 - cos(phase_difference)),
                sin(2 * theta) ** 2 + cos(2 * theta) ** 2 * cos(phase_difference),
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


def calculate_depolarizer_matrix(wavelength: float):
    """Calculate the depolarizer matrix

    Parameters
    ----------
        wavelength (float): wavelength in nm

    Returns
    -------
        DEPOLARIZER_MATRIX (ndarray): depolarizer matrix for the provided wavelength
    """

    phase_shift = np.deg2rad(_calculate_phase_shift(wavelength))
    DEPOLARIZER_MATRIX = np.asarray(
        [
            [1, 0, 0, 0],
            [
                0,
                cos(2 * phase_shift),
                sin(2 * phase_shift) * sin(phase_shift),
                sin(2 * phase_shift) * cos(phase_shift),
            ],
            [0, 0, cos(phase_shift), -sin(phase_shift)],
            [
                0,
                -sin(2 * phase_shift),
                cos(2 * phase_shift) * sin(phase_shift),
                cos(2 * phase_shift) * sin(phase_shift),
            ],
        ],
        dtype=np.float64,
    )
    return DEPOLARIZER_MATRIX


def _calculate_phase_shift(wavelength: float):
    a, b = -32878.57142857143, 45656261.9047619
    phase_shift = a * wavelength + b
    return phase_shift


def interpolation_func(x, a, b, c, d):
    return a * x**3 + b * x**2 + c * x + d


def apply_matrix(matrix, sed):
    for idx, _ in enumerate(sed[0]):
        sed[:, idx] = np.transpose(matrix.dot(sed[:, idx]))
    return sed

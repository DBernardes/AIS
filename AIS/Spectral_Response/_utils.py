from math import sin, cos, pi
import numpy as np
from numpy import ndarray

__all__ = ["POLARIZER_MATRIX", "POLARIZER_90_MATRIX",
           "calculate_retarder_matrix"]

_THETA_POL = np.deg2rad(0)
POLARIZER_MATRIX = 0.5 * np.asarray(
    [
        [1, cos(2 * _THETA_POL), sin(2 * _THETA_POL), 0],
        [
            cos(2 * _THETA_POL),
            cos(2 * _THETA_POL) ** 2,
            cos(2 * _THETA_POL) * sin(2 * _THETA_POL),
            0,
        ],
        [
            sin(2 * _THETA_POL),
            cos(2 * _THETA_POL) * sin(2 * _THETA_POL),
            sin(2 * _THETA_POL) ** 2,
            0,
        ],
        [0, 0, 0, 0],
    ]
)

_THETA_POL = np.deg2rad(_THETA_POL + pi / 2)
POLARIZER_90_MATRIX = 0.5 * np.asarray(
    [
        [1, cos(2 * _THETA_POL), sin(2 * _THETA_POL), 0],
        [
            cos(2 * _THETA_POL),
            cos(2 * _THETA_POL) ** 2,
            cos(2 * _THETA_POL) * sin(2 * _THETA_POL),
            0,
        ],
        [
            sin(2 * _THETA_POL),
            cos(2 * _THETA_POL) * sin(2 * _THETA_POL),
            sin(2 * _THETA_POL) ** 2,
            0,
        ],
        [0, 0, 0, 0],
    ]
)


def calculate_retarder_matrix(phase_difference) -> ndarray:
    """Calculate the retarder matrix for a given phase difference.

    Parameters
    ----------
    phase_difference : float
        The phase difference in degrees.
    """
    phase_difference = np.deg2rad(phase_difference)
    retarder_matrix = np.asarray(
        [
            [1, 0, 0, 0],
            [
                0,
                cos(2 * _THETA_POL) ** 2
                + sin(2 * _THETA_POL) ** 2 * cos(phase_difference),
                cos(2 * _THETA_POL)
                * sin(2 * _THETA_POL)
                * (1 - cos(phase_difference)),
                -sin(2 * _THETA_POL) * sin(phase_difference),
            ],
            [
                0,
                cos(2 * _THETA_POL)
                + sin(2 * _THETA_POL) * (1 - cos(phase_difference)),
                sin(2 * _THETA_POL) ** 2
                + cos(2 * _THETA_POL) ** 2 * cos(phase_difference),
                cos(2 * _THETA_POL) * sin(phase_difference),
            ],
            [
                0,
                sin(2 * _THETA_POL) * sin(phase_difference),
                -cos(2 * _THETA_POL) * sin(phase_difference),
                cos(phase_difference),
            ],
        ]
    )

    return retarder_matrix


def main() -> None:
    print(*['', POLARIZER_90_MATRIX, POLARIZER_MATRIX,
          calculate_retarder_matrix(45), ''], sep="\n\n")


if __name__ == "__main__":
    main()

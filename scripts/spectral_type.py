import numpy as np


spectral_type_definitions = {
    "O6": [-5.10, -5.00],
    "O7": [-5.00, -4.75],
    "O8": [-4.75, -4.00],
    "B0": [-4.00, -3.00],
    "B1": [-3.00, -2.10],
    "B2": [-2.10, -1.45],
    "B3": [-1.45, -0.90],
    "B5": [-0.90, -0.30],
    "B6": [-0.30, 0.10],
    "B7": [0.10, 0.45],
    "B8": [0.45, 0.70],
    "B9": [0.70, 0.95],
    "A0": [0.95, 1.20],
    "A1": [1.20, 1.40],
    "A2": [1.40, 1.85],
    "A5": [1.85, 2.45],
    "A8": [2.45, 2.85]
}


def AbsMag(self, m_V, distance=3.0):
    # distance is in kpc: calculated mean/median distance to cluster
    M_V = m_V - 5 * np.log10(100 * distance)

    return M_V


def GetSpectralType(self, m_V, distance=3.0):
    M_V = self.AbsMag(m_V, distance)

    for elem in list(spectral_type_definitions):
        if spectral_type_definitions[elem][0] <= M_V < spectral_type_definitions[elem][1]:
            spectral_type = elem
            break
    else:
        # print("Error: Spectral type not found for a magnitude of M_V = " + str(M_V))
        spectral_type = "--"

    return spectral_type

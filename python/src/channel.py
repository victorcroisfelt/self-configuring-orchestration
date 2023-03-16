import numpy as np
from scipy.constants import speed_of_light


def scenario(wavelength, scenario_size, M, N):
    """

    :param scenario_size:
    :param M:
    :param N:
    :return:
    """

    # Coordinates of the RIS
    pos_ris = 0

    # Coordinates of RIS elements
    pos_ris_els = pos_ris + (wavelength / 2) * np.arange(N)
    pos_ris_els = pos_ris_els - ((pos_ris_els[-1] - pos_ris_els[0]) / 2)

    # Convert to vector form
    pos_ris = np.array([pos_ris.real, pos_ris.imag])
    pos_ris_els = np.array([pos_ris_els.real, pos_ris_els.imag]).T

    # Coordinates of the BS
    pos_bs = - 1 * scenario_size + 1j * scenario_size / 2

    # Coordinates of BS elements
    pos_bs_els = pos_bs + 1j * (wavelength / 2) * np.arange(M)
    pos_bs_els = pos_bs_els - ((pos_bs_els[-1] - pos_bs_els[0]) / 2)

    # Convert to vector form
    pos_bs = np.array([pos_bs.real, pos_bs.imag])
    pos_bs_els = np.array([pos_bs_els.real, pos_bs_els.imag]).T

    # Array steering vectors
    bs_ris_steering = array_steering_vector(wavelength, pos_ris, pos_bs, pos_bs_els)
    ris_bs_steering = array_steering_vector(wavelength, pos_bs, pos_ris, pos_ris_els)

    # BS-HRIS pathloss
    bs_ris_pathloss = pathloss(2, pos_bs, pos_ris)

    # BS-HRIS channel
    bs_ris_channels = np.sqrt(bs_ris_pathloss) * np.outer(bs_ris_steering, ris_bs_steering.conj())

    return pos_bs, pos_bs_els, pos_ris, pos_ris_els, bs_ris_channels, ris_bs_steering


def array_steering_vector(wavelength, pos1, pos2, pos_els):
    """

    :param wavelength:
    :param pos1:
    :param pos2:
    :param pos_els:
    :return:
    """

    # Compute wave vector
    wave_vector = (2 * np.pi / wavelength) * (pos1 - pos2) / np.linalg.norm(pos2 - pos1)

    # Compute inner product
    if len(pos1.shape) == 1:
        angles = wave_vector * (pos_els - pos2)
        angles = angles.sum(axis=-1)
    else:
        angles = wave_vector[:, None, :] * (pos_els - pos2)
        angles = angles.sum(axis=-1)

    # Compute steering vector
    steering = np.exp(1j * angles)

    return steering


def pathloss(pathloss_exp, pos1, pos2):
    """

    :param pathloss_exp:
    :param pos1:
    :param pos2:
    :return:
    """

    # Constant term
    # constant_term_db = -35.3
    # constant_term = 10**(constant_term_db/10)

    # Compute distances
    distances = np.linalg.norm(pos1 - pos2, axis=-1)

    # Compute pathloss
    #pathloss = constant_term * (1 / distances) ** pathloss_exp
    pathloss = (1 / distances) ** pathloss_exp

    return pathloss


def generate_channel_realizations(wavelength, pos_bs, pos_bs_els, pos_ris, pos_ris_els, pos_ues, sigma2_dr, sigma2_rr, n_channels):
    """

    :param wavelength:
    :param pathloss_exp:
    :param pos_bs:
    :param pos_bs_els:
    :param pos_ris:
    :param pos_ris_els:
    :param pos_ues:
    :return:
    """

    # BS-UE channels
    bs_ue_steering = array_steering_vector(wavelength, pos_ues, pos_bs, pos_bs_els)
    bs_ue_pathloss = pathloss(3.76, pos_bs, pos_ues)

    los_bs_ue_channels = np.sqrt(bs_ue_pathloss)[:, None] * bs_ue_steering

    # RIS-UE channels
    ris_ue_steering = array_steering_vector(wavelength, pos_ues, pos_ris, pos_ris_els)
    ris_ue_pathloss = pathloss(2, pos_ris, pos_ues)

    los_ris_ue_channels = np.sqrt(ris_ue_pathloss)[:, None] * ris_ue_steering

    # Generate realizations
    bs_ue_channels = los_bs_ue_channels[:, :, None] + np.sqrt(sigma2_dr/2) * (np.random.randn(n_channels) + 1j * np.random.randn(n_channels))
    ris_ue_channels = los_ris_ue_channels[:, :, None] + np.sqrt(sigma2_rr/2) * (np.random.randn(n_channels) + 1j * np.random.randn(n_channels))

    return bs_ue_channels, los_bs_ue_channels, ris_ue_channels, los_ris_ue_channels


def drop_ues(K, pos_ris, scenario_size, guard_distance=0.0):

    # Position RIS to complex
    pos_ris = pos_ris[0] + 1j * pos_ris[1]

    # Prepare to save UE positions
    pos_ues = np.zeros(K, dtype=np.complex_)

    # Drop UEs
    for kk in range(K):

        # Drop UE
        while True:
            pos_ue = (scenario_size * np.random.rand() + 1j * scenario_size * np.random.rand()) + 100
            test = np.abs(pos_ue - pos_ris) > guard_distance

            if test:
                pos_ues[kk] = pos_ue
                break

    # Conversion to vector form
    pos_ues = np.array([pos_ues.real, pos_ues.imag]).T

    return pos_ues


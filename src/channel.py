import numpy as np
from scipy.constants import speed_of_light


def scenario(wavelength, M, N):
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

    # Guard distance from the RIS
    guard_distance_ris = 2 / wavelength * (2 * pos_ris_els[-1, 0])**2

    # Coordinates of BS elements
    pos_bs_els = (wavelength / 2) * np.arange(M)
    pos_bs_els = pos_bs_els - ((pos_bs_els[-1] - pos_bs_els[0]) / 2)

    # Guard distance from the BS
    guard_distance_bs = 2 / wavelength * (2 * pos_bs_els[-1].imag) ** 2

    # Coordinates of the BS
    pos_bs = 1000 * (- np.cos(np.pi/2) + 1j * np.sin(np.pi/2))

    # Convert to vector form
    pos_bs = np.array([pos_bs.real, pos_bs.imag])
    pos_bs_els = np.array([pos_bs_els.real, pos_bs_els.imag]).T
    pos_bs_els += pos_bs

    # Array steering vectors
    bs_ris_steering = array_steering_vector(wavelength, pos_ris, pos_bs, pos_bs_els)
    ris_bs_steering = array_steering_vector(wavelength, pos_bs, pos_ris, pos_ris_els)
     
    # BS-HRIS pathloss
    bs_ris_pathloss = pathloss(2, pos_bs, pos_ris)

    # BS-HRIS channel
    bs_ris_channels = np.sqrt(bs_ris_pathloss) * np.outer(bs_ris_steering, ris_bs_steering.conj())

    return pos_bs, pos_bs_els, pos_ris, pos_ris_els, bs_ris_channels, ris_bs_steering, guard_distance_ris


def drop_ues(K, pos_ris, dmax, dmin=None, guard_distance_ris=0.0):

    # Position RIS to complex
    pos_ris = pos_ris[0] + 1j * pos_ris[1]

    if dmin is None:
        dmin = guard_distance_ris

    # Generate distances
    distances = np.sqrt(np.random.rand(K) * (dmax**2 - dmin**2) + dmin**2)

    # Generate angles
    angles = np.pi/2 * np.random.rand(K)

    # Conversion to vector form
    pos_ues = distances[:, None] * np.array([np.cos(angles), np.sin(angles)]).T

    return pos_ues


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

    # BS-UE LoS channels
    bs_ue_steering = array_steering_vector(wavelength, pos_ues, pos_bs, pos_bs_els)
    bs_ue_pathloss = pathloss(3.76, pos_bs, pos_ues)

    los_bs_ue_channels = np.sqrt(bs_ue_pathloss)[:, None] * bs_ue_steering

    # RIS-UE LoS channels
    ris_ue_steering = array_steering_vector(wavelength, pos_ues, pos_ris, pos_ris_els)
    ris_ue_pathloss = pathloss(2, pos_ris, pos_ues)

    los_ris_ue_channels = np.sqrt(ris_ue_pathloss)[:, None] * ris_ue_steering
    
    # Generate realizations
    random_ = np.sqrt(sigma2_dr/2) *  (np.random.randn(n_channels) + 1j * np.random.randn(n_channels))
    bs_ue_channels = los_bs_ue_channels[:, :, None] + random_[None, None, :]

    random_ = np.sqrt(sigma2_rr/2) * (np.random.randn(n_channels) + 1j * np.random.randn(n_channels))
    ris_ue_channels = los_ris_ue_channels[:, :, None] + random_[None, None, :]

    # Correct
    bs_ue_channels = bs_ue_channels.transpose((1, 2, 0))
    ris_ue_channels = ris_ue_channels.transpose((1, 2, 0))

    return bs_ue_channels, ris_ue_channels

def generate_los_channel_realizations(wavelength, pos_bs, pos_bs_els, pos_ris, pos_ris_els, pos_ues):
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

    # BS-UE LoS channels
    bs_ue_steering = array_steering_vector(wavelength, pos_ues, pos_bs, pos_bs_els)
    bs_ue_pathloss = pathloss(3.76, pos_bs, pos_ues)

    los_bs_ue_channels = np.sqrt(bs_ue_pathloss)[:, None] * bs_ue_steering

    # RIS-UE LoS channels
    ris_ue_steering = array_steering_vector(wavelength, pos_ues, pos_ris, pos_ris_els)
    ris_ue_pathloss = pathloss(2, pos_ris, pos_ues)

    los_ris_ue_channels = np.sqrt(ris_ue_pathloss)[:, None] * ris_ue_steering
    
    # Correct
    los_bs_ue_channels = los_bs_ue_channels.transpose((1, 0))
    los_ris_ue_channels = los_ris_ue_channels.transpose((1, 0))

    return los_bs_ue_channels, los_ris_ue_channels

def array_steering_vector(wavelength, pos1, pos2, pos_els):
    """

    :param wavelength:
    :param pos1:
    :param pos2:
    :param pos_els:
    :return:
    """

    # Compute inner product
    if len(pos1.shape) == 1:
        wave_vector = (2 * np.pi / wavelength) * (pos1 - pos2) / np.linalg.norm(pos2 - pos1)
        angles = wave_vector * (pos_els - pos2)
        angles = angles.sum(axis=-1)
    else:
        wave_vector = (2 * np.pi / wavelength) * (pos1 - pos2) / np.linalg.norm(pos2 - pos1, axis=1)[:, None]
        angles = wave_vector[:, None, :] * (pos_els - pos2)
        angles = angles.sum(axis=-1 )
    
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


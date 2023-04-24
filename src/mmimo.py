import numpy as np

from src.channel import array_steering_vector
from scipy.stats.distributions import chi2


def bs_rx_chest(P_ue, n_pilots, sigma2_n, n_pilot_subblocks, n_pilot_subblocks_probe,
                eq_channels_probe, eq_channels):
    """
    Return estimated equivalent channel.

    :param P_ue:
    :param n_pilots:
    :param sigma2_n:
    :param n_pilot_subblocks:
    :param n_pilot_subblocks_probe:
    :param eq_channels_probe:
    :param eq_channels:
    :return:
    """

    # Extract useful constants
    K, M, n_channels = eq_channels.shape

    # Generate estimation noise
    noise = np.random.randn(K, M, n_channels) + 1j * np.random.randn(K, M, n_channels)
    noise *= np.sqrt(sigma2_n / 2 / n_pilot_subblocks / P_ue / n_pilots)

    # Compute average equivalent channel
    avg_eq_channels = eq_channels_probe + (n_pilot_subblocks - n_pilot_subblocks_probe) * eq_channels
    avg_eq_channels *= (1 / n_pilot_subblocks)

    # Get equivalent channel estimates
    hat_eq_channels = avg_eq_channels + noise

    return hat_eq_channels


def bs_rx_chest_no_probe(P_ue, n_pilots, sigma2_n, n_pilot_subblocks, n_pilot_subblocks_probe, eq_channels):
    """
    Return estimated equivalent channel.

    :param P_ue:
    :param n_pilots:
    :param sigma2_n:
    :param n_pilot_subblocks:
    :param n_pilot_subblocks_probe:
    :param eq_channels_probe:
    :param eq_channels:
    :return:
    """

    # Extract useful constants
    K, M, n_channels = eq_channels.shape

    # Generate estimation noise
    noise = np.random.randn(K, M, n_channels) + 1j * np.random.randn(K, M, n_channels)
    noise *= np.sqrt(sigma2_n / 2 / (n_pilot_subblocks-n_pilot_subblocks_probe) / P_ue / n_pilots)

    # Get equivalent channel estimates
    hat_eq_channels = eq_channels + noise

    return hat_eq_channels


def bs_comm(P_ue, sigma2_n, channels, hat_channels, method='MR'):
    """

    :param P_ue:
    :param sigma2_n:
    :param channels:
    :param hat_channels:
    :return:
    """

    # Extract useful constants
    K, M, n_channels = channels.shape

    if method == 'MR':
        comb_vec = hat_channels
    elif method == 'ZF':
        comb_vec = hat_channels.transpose(2, 1, 0)
        comb_vec = np.matmul(comb_vec.conj().transpose(0, 2, 1), comb_vec)
        comb_vec = np.matmul(hat_channels.transpose(2, 1, 0), np.linalg.inv(comb_vec))
        comb_vec = comb_vec.transpose(2, 1, 0)

    # Compute numerator of the SINR
    num = P_ue * np.abs((comb_vec.conj() * channels).sum(axis=1))**2

    # Compute denominator 1
    den1 = P_ue * np.abs((comb_vec[:, None, :, :].conj() * channels[None, :, :, :]).sum(axis=2))**2
    den1 = den1.sum(axis=1)
    den1 -= num

    # Generate some noise
    noise = np.sqrt(sigma2_n / 2) * (np.random.randn(K, M, n_channels) + 1j * np.random.randn(K, M, n_channels))

    # Compute denominator 2
    den2 = np.abs((comb_vec.conj() * noise).sum(axis=1))**2
    #breakpoint()
    # Compute SINR
    sinr = num / (den1 + den2)

    # Compute SE
    se = np.log2(1 + sinr)

    return se, num, den1, den2







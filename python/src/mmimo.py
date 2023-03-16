import numpy as np

from src.channel import array_steering_vector
from scipy.stats.distributions import chi2


def bs_rx_chest(eta, P_ue, n_pilots, sigma2_n, n_pilot_subblocks, n_pilot_subblocks_probe,
                eq_channels_probe, eq_channels):

    """

    :param eta:
    :param P_ue:
    :param n_pilots:
    :param sigma2_n:
    :param n_pilot_subblocks:
    :param n_pilot_subblocks_probe:
    :param bs_ue_channels:
    :param bs_ris_channels:
    :param ris_ue_channels:
    :param reflection_configs:
    :param pow_probe_configs:
    :return:
    """

    # Extract useful constants
    K, M, n_channels = eq_channels.shape

    # Generate noise after decorrelation
    noise = (np.random.randn(n_pilot_subblocks, K, M, n_channels) + 1j * np.random.randn(n_pilot_subblocks, K, M, n_channels))
    noise *= np.sqrt((n_pilots * sigma2_n) / 2)

    # Compute average equivalent channel
    avg_eq_channels = (1/n_pilot_subblocks) * (eq_channels_probe + (n_pilot_subblocks - n_pilot_subblocks_probe) * eq_channels)

    # Prepare to save received signal
    rx_signal = np.sqrt(P_ue) * n_pilots * avg_eq_channels + noise

    # Take the average over subblocks
    average_subblocks = np.mean(rx_signal, axis=0)

    # Get equivalent channel estimates
    hat_eq_channels = average_subblocks / np.sqrt(P_ue) / n_pilots

    return hat_eq_channels


def bs_comm(P_ue, sigma2_n, channels, hat_channels):
    """

    :param P_ue:
    :param sigma2_n:
    :param channels:
    :param hat_channels:
    :return:
    """

    # Extract useful constants
    K, M, n_channels = channels.shape

    # Compute numerator of the SINR
    num = P_ue * np.abs((hat_channels.conj() * channels).sum(axis=1))**2

    # Compute denominator 1
    den1 = P_ue * np.abs((hat_channels[:, None, :, :].conj() * channels[None, :, :, :]).sum(axis=2))**2
    den1 = den1.sum(axis=1)
    den1 -= num

    # Generate some noise
    noise = np.sqrt(sigma2_n / 2) * (np.random.randn(K, M, n_channels) + 1j * np.random.randn(K, M, n_channels))

    # Compute denominator 2
    den2 = np.abs((hat_channels.conj() * noise).sum(axis=1))**2

    # Compute SINR
    sinr = num / (den1 + den2)

    # Compute SE
    se = np.log2(1 + sinr)

    return se, num, den1, den2







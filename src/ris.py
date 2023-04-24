import numpy as np

from src.channel import array_steering_vector
from scipy.stats.distributions import chi2


def pow_ris_config_codebook(wavelength, n_probe_pilot_subblocks, pos_ris, pos_ris_els):
    """

    :param wavelength:
    :param n_probe_pilot_subblocks:
    :param pos_ris:
    :param pos_ris_els:
    :return:
    """

    # Slice vector space into C sectors from 0 to pi/2
    if n_probe_pilot_subblocks == 1:
        sliced_angular_space = np.array([np.pi / 2])
    else:
        sliced_angular_space = np.linspace(0, np.pi, n_probe_pilot_subblocks)

    # Compute position vector with respect to the sliced angular space
    pos_slices = np.array([np.cos(sliced_angular_space), np.sin(sliced_angular_space)]).T

    # Compute probe_configurations
    power_probe_configs = array_steering_vector(wavelength, pos_slices, pos_ris, pos_ris_els)

    return power_probe_configs


def ris_rx_chest(eta, P_ue, n_pilots, sigma2_n, n_probe_pilot_subblocks, ris_ue_channels, power_probe_configs=None):
    """

    :param eta:
    :param P_ue:
    :param n_pilots:
    :param sigma2_n:
    :param n_probe_pilot_subblocks:
    :param ris_ue_channels:
    :param power_probe_configs:
    :return:
    """

    # Extract useful constants
    K, N, n_channels = ris_ue_channels.shape

    if power_probe_configs is not None:

        # Generate noise
        noise = np.random.randn(n_probe_pilot_subblocks, K, n_channels) + 1j * np.random.randn(n_probe_pilot_subblocks, K, n_channels)
        noise *= np.sqrt(N * sigma2_n / 2)

        # Compute equivalent channel
        equivalent_channel = (power_probe_configs.conj()[:, None, :, None] * ris_ue_channels[None, :, :, :]).sum(axis=2)

        # Compute received signal
        ris_rx_chest = np.sqrt((1 - eta)) * np.sqrt(P_ue) * np.sqrt(n_pilots) * equivalent_channel + noise

    else:

        # Generate noise
        noise = (np.random.randn(n_probe_pilot_subblocks, K, N, n_channels) + 1j * np.random.randn(n_probe_pilot_subblocks, K, N, n_channels))
        noise *= np.sqrt(n_pilots * sigma2_n / 2)

        # Compute received signal
        ris_rx_chest = np.sqrt((1 - eta)) * np.sqrt(P_ue) * n_pilots * ris_ue_channels[None, :, :, :] + noise

    return ris_rx_chest


def gen_ris_probe(ris_ue_channels):

    # Compute weights
    weights = np.linalg.norm(ris_ue_channels, axis=1)
    weights = weights / weights.sum(axis=0)

    # Compute reflection configuration
    reflection_configs = weights[:, None, :] * (np.exp(1j * np.angle(ris_ue_channels)).conj())
    reflection_configs = reflection_configs.sum(axis=0)

    return reflection_configs, weights


def pow_ris_probe(N, sigma2_n, proba_false_alarm, power_ris_rx_chest, power_probe_configs):
    """

    :param N:
    :param sigma2_n:
    :param proba_false_alarm:
    :param power_ris_rx_chest:
    :param power_probe_configs:
    :return:
    """

    # Compute measured power
    measured_power = np.abs(power_ris_rx_chest)**2

    # Compute threshold
    threshold = - (2 * N * sigma2_n) * np.log(proba_false_alarm)

    # Perform test.py
    test = measured_power > threshold

    # Apply OR over subblock dimension
    or_test = np.sum(test, axis=0)

    # Detected UEs
    detected_ues = or_test > 0

    # Average probability of detection
    pd = detected_ues.mean()

    # Store detected UEs
    measured_power[:, ~detected_ues] = 0

    # Get the maximum
    max_measured_power = measured_power.max(axis=0)

    # Compute weights
    weights = np.sqrt(max_measured_power) / np.sum(np.sqrt(max_measured_power), axis=0)

    # Extract relevant constants
    K, n_channels = detected_ues.shape

    # Prepare to store results
    hat_aoa = np.zeros((K, N, n_channels), dtype=np.complex_)
    reflection_configs = np.zeros((N, n_channels), dtype=np.complex_)

    hat_aoa[:] = np.nan + 1j*np.nan

    # Go through all channel realizations
    for rr in range(n_channels):

        # Go through all detected UEs
        for dd in range(K):

            # Check if UE was detected
            if detected_ues[dd, rr]:

                # Find the best configuration
                hat_aoa[dd, :, rr] = power_probe_configs[np.argmax(measured_power[:, dd, rr]), :]

                # Compute reflection configuration
                reflection_configs[:, rr] += weights[dd, rr] * hat_aoa[dd, :, rr].conj()

    return reflection_configs, weights, hat_aoa, pd


def sig_ris_probe(n_pilots, sigma2_n, proba_false_alarm, signal_ris_rx_chest):
    """

    :param sigma2_n:
    :param proba_false_alarm:
    :param signal_ris_rx_chest:
    :return:
    """

    # Extract useful constants
    n_probe_pilot_subblocks, K, N, n_channels = signal_ris_rx_chest.shape

    # Average signal over probe subblocks
    avg_subblocks = np.mean(signal_ris_rx_chest, axis=0)

    # Average signal over antennas
    avg_antennas = np.mean(avg_subblocks, axis=1)

    # Get threshold value
    threshold = chi2.ppf(1 - proba_false_alarm, df=4)

    # Compute mu (test.py variable)
    mu = np.abs(avg_antennas) ** 2
    mu /= ((n_pilots * sigma2_n) / (2 * N * n_probe_pilot_subblocks))

    # Perform test.py
    detected_ues = mu > threshold

    # Average probability of detection
    pd = np.mean(detected_ues)

    # Update avg antennas
    avg_antennas[~detected_ues] = 0

    # Compute weights
    weights = np.linalg.norm(avg_subblocks, axis=1)
    weights = weights / weights.sum(axis=0)

    # Prepare to store the reflection configuration
    hat_aoa = np.zeros((K, N, n_channels), dtype=np.complex_)
    reflection_configs = np.zeros((N, n_channels), dtype=np.complex_)

    hat_aoa[:] = np.nan

    # Go through all channel realizations
    for rr in range(n_channels):

        # Go through all detected UEs
        for dd in range(K):

            # Check if UE was detected
            if detected_ues[dd, rr]:

                # Find the best configuration
                hat_aoa[dd, :, rr] = np.exp(1j * np.angle(avg_subblocks[dd, :, rr]))

                # Compute reflection configuration
                reflection_configs[:, rr] += weights[dd, rr] * hat_aoa[dd, :, rr].conj()

    return reflection_configs, weights, hat_aoa, pd


# def pow_ris_probe(N, sigma2_n, proba_false_alarm, power_ris_rx_chest, power_probe_configs):
#     """
#
#     :param N:
#     :param sigma2_n:
#     :param proba_false_alarm:
#     :param power_ris_rx_chest:
#     :param power_probe_configs:
#     :return:
#     """
#
#     # Compute measured power
#     measured_power = np.abs(power_ris_rx_chest)**2
#
#     # Compute threshold
#     threshold = - (2 * N * sigma2_n) * np.log(proba_false_alarm)
#
#     # Perform test.py
#     test = measured_power > threshold
#
#     # Apply OR over subblock dimension
#     or_test = np.sum(test, axis=0)
#
#     # Detected UEs
#     detected_ues = or_test > 0
#
#     # Extract relevant constants
#     K, n_channels = detected_ues.shape
#
#     # Prepare to store the reflection configuration
#     reflection_configs = np.ones((N, n_channels), dtype=np.complex_)
#
#     # Go through all channel realizations
#     for cc in range(n_channels):
#
#         # Go through all pilots
#         for pp in range(K):
#
#             if detected_ues[pp, cc]:
#
#                 # Find the best configuration
#                 hat_aoa = power_probe_configs[np.argmax(measured_power[:, pp, cc]), :]
#
#                 # Compute reflection configuration
#                 reflection_configs[:, cc] *= hat_aoa.conj()
#
#     return reflection_configs.conj()

# def pow_ris_probe3(N, sigma2_n, proba_false_alarm, power_ris_rx_chest, power_probe_configs):
#     """
#
#     :param N:
#     :param sigma2_n:
#     :param proba_false_alarm:
#     :param power_ris_rx_chest:
#     :param power_probe_configs:
#     :return:
#     """
#
#     # Compute measured power
#     measured_power = np.abs(power_ris_rx_chest)**2
#
#     # Compute threshold
#     threshold = - (2 * N * sigma2_n) * np.log(proba_false_alarm)
#
#     # Perform test.py
#     test = measured_power > threshold
#
#     # Apply OR over subblock dimension
#     or_test = np.sum(test, axis=0)
#
#     # Detected UEs
#     detected_ues = or_test > 0
#
#     # Extract relevant constants
#     K, n_channels = detected_ues.shape
#
#     # Prepare to store the reflection configuration
#     reflection_configs = np.zeros((N, n_channels), dtype=np.complex_)
#
#     # Number of elements per UE
#     N_per_ue = N // K
#
#     # Go through all channel realizations
#     for cc in range(n_channels):
#
#         # Go through all pilots
#         for pp in range(K):
#
#             if detected_ues[pp, cc]:
#
#                 # Find the best configuration
#                 hat_aoa = power_probe_configs[np.argmax(measured_power[:, pp, cc]), :]
#
#                 # Compute reflection configuration
#                 reflection_configs[pp * N_per_ue:(1 + pp) * N_per_ue, cc] = hat_aoa[pp * N_per_ue:(1 + pp) * N_per_ue].conj()
#
#     return reflection_configs


#

# def sig_ris_probe(n_pilots, sigma2_n, proba_false_alarm, signal_ris_rx_chest):
#     """
#
#     :param sigma2_n:
#     :param proba_false_alarm:
#     :param signal_ris_rx_chest:
#     :return:
#     """
#
#     # Extract useful constants
#     n_probe_pilot_subblocks, K, N, n_channels = signal_ris_rx_chest.shape
#
#     # Average signal over probe subblocks
#     avg_subblocks = np.mean(signal_ris_rx_chest, axis=0)
#
#     # Average signal over antennas
#     avg_antennas = np.mean(avg_subblocks, axis=1)
#
#     # Get threshold value
#     threshold = chi2.ppf(1 - proba_false_alarm, df=4)
#
#     # Compute mu (test.py variable)
#     mu = np.abs(avg_antennas) ** 2
#     mu /= ((n_pilots * sigma2_n) / (2 * N * n_probe_pilot_subblocks))
#
#     # Perform test.py
#     detected_ues = mu > threshold
#
#     # Prepare to store the reflection configuration
#     reflection_configs = np.ones((N, n_channels), dtype=np.complex_)
#
#     # Go through all channel realizations
#     for cc in range(n_channels):
#
#         # Go through all pilots
#         for pp in range(K):
#
#             if detected_ues[pp, cc]:
#
#                 # Find the best configuration
#                 hat_aoa = np.exp(1j * np.angle(avg_subblocks[pp, :, cc]))
#
#                 # Compute reflection configuration
#                 reflection_configs[:, cc] *= hat_aoa.conj()
#
#     return reflection_configs

# def sig_ris_probe3(n_pilots, sigma2_n, proba_false_alarm, signal_ris_rx_chest):
#     """
#
#     :param sigma2_n:
#     :param proba_false_alarm:
#     :param signal_ris_rx_chest:
#     :return:
#     """
#
#     # Extract useful constants
#     n_probe_pilot_subblocks, K, N, n_channels = signal_ris_rx_chest.shape
#
#     # Average signal over probe subblocks
#     avg_subblocks = np.mean(signal_ris_rx_chest, axis=0)
#
#     # Average signal over antennas
#     avg_antennas = np.mean(avg_subblocks, axis=1)
#
#     # Get threshold value
#     threshold = chi2.ppf(1 - proba_false_alarm, df=4)
#
#     # Compute mu (test.py variable)
#     mu = np.abs(avg_antennas) ** 2
#     mu /= ((n_pilots * sigma2_n) / (2 * N * n_probe_pilot_subblocks))
#
#     # Perform test.py
#     detected_ues = mu > threshold
#
#     # Prepare to store the reflection configuration
#     reflection_configs = np.zeros((N, n_channels), dtype=np.complex_)
#
#     # Number of elements per UE
#     N_per_ue = N // K
#
#     # Go through all channel realizations
#     for cc in range(n_channels):
#
#         # Go through all pilots
#         for pp in range(K):
#
#             if detected_ues[pp, cc]:
#
#                 # Find the best configuration
#                 hat_aoa = np.exp(1j * np.angle(avg_subblocks[pp, :, cc]))
#
#                 # Compute reflection configuration
#                 reflection_configs[pp * N_per_ue:(1 + pp) *  N_per_ue, cc] = hat_aoa[pp * N_per_ue:(1 + pp) * N_per_ue].conj()
#
#     return reflection_configs
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
        sliced_angular_space = np.linspace(0, np.pi/2, n_probe_pilot_subblocks)

    # Compute position vector with respect to the sliced angular space
    pos_slices = np.array([np.cos(sliced_angular_space), np.sin(sliced_angular_space)]).T

    # Compute probe_configurations
    pd_probe_configs = array_steering_vector(wavelength, pos_slices, pos_ris, pos_ris_els)

    return pd_probe_configs


def ris_rx_chest(eta, P_ue, n_pilots, sigma2_n, n_probe_pilot_subblocks, ris_ue_channels, pd_probe_configs=None):
    """

    :param eta:
    :param P_ue:
    :param n_pilots:
    :param sigma2_n:
    :param n_probe_pilot_subblocks:
    :param ris_ue_channels:
    :param pd_probe_configs:
    :return:
    """

    # Extract useful constants
    N, n_channels, K = ris_ue_channels.shape

    if pd_probe_configs is not None:

        # Generate noise
        noise = np.random.randn(n_probe_pilot_subblocks, n_channels, K) + 1j * np.random.randn(n_probe_pilot_subblocks, n_channels, K)
        noise *= np.sqrt(N * sigma2_n / 2)

        # Compute equivalent channel
        equivalent_channel = (pd_probe_configs.conj()[:, :, None, None] * ris_ue_channels[None, :, :, :]).sum(axis=1)

        # Compute received signal
        ris_rx_chest = noise
        ris_rx_chest[:, :, :] += np.sqrt((1 - eta)) * np.sqrt(P_ue) * np.sqrt(n_pilots) * equivalent_channel[:, :, :]

    else:

        # Generate noise
        noise = (np.random.randn(n_probe_pilot_subblocks, N, n_channels, K) + 1j * np.random.randn(n_probe_pilot_subblocks, N, n_channels, K))
        noise *= np.sqrt(n_pilots * sigma2_n / 2)

        # Compute received signal
        ris_rx_chest = noise
        ris_rx_chest[:, :, :, :] += np.sqrt((1 - eta)) * np.sqrt(P_ue) * n_pilots * ris_ue_channels[None, :, :, :]

    return ris_rx_chest


def ris_rx_chest_with_choice(eta, P_ue, n_pilots, sigma2_n, n_probe_pilot_subblocks, ris_ue_channels, mask, pd_probe_configs=None):
    """

    :param eta:
    :param P_ue:
    :param n_pilots:
    :param sigma2_n:
    :param n_probe_pilot_subblocks:
    :param ris_ue_channels:
    :param pd_probe_configs:
    :return:
    """

    # Extract useful constants
    N, n_channels, K = ris_ue_channels.shape

    if pd_probe_configs is not None:

        # Generate noise
        noise = np.random.randn(n_probe_pilot_subblocks, n_channels, K) + 1j * np.random.randn(n_probe_pilot_subblocks, n_channels, K)
        noise *= np.sqrt(N * sigma2_n / 2)

        # Compute equivalent channel
        equivalent_channel = (pd_probe_configs.conj()[:, :, None, None] * ris_ue_channels[None, :, :, :]).sum(axis=1)

        # Compute received signal
        ris_rx_chest = noise
        ris_rx_chest[:, :, mask] += np.sqrt((1 - eta)) * np.sqrt(P_ue) * np.sqrt(n_pilots) * equivalent_channel[:, :, mask]

    else:

        # Generate noise
        noise = (np.random.randn(n_probe_pilot_subblocks, N, n_channels, K) + 1j * np.random.randn(n_probe_pilot_subblocks, N, n_channels, K))
        noise *= np.sqrt(n_pilots * sigma2_n / 2)

        # Compute received signal
        ris_rx_chest = noise
        ris_rx_chest[:, :, :, mask] += np.sqrt((1 - eta)) * np.sqrt(P_ue) * n_pilots * ris_ue_channels[None, :, :, mask]

    return ris_rx_chest


def gen_ris_probe(ris_ue_channels):

    # Compute weights
    weights = np.linalg.norm(ris_ue_channels, axis=0)
    weights = weights / weights.sum(axis=-1)[:, None]

    # Compute reflection configuration
    reflection_configs = weights[None, :, :] * (np.exp(1j * np.angle(ris_ue_channels)).conj())
    reflection_configs = reflection_configs.sum(axis=-1)

    return reflection_configs, weights


def pow_ris_probe(N, sigma2_n, proba_false_alarm, pd_ris_rx_chest, pd_probe_configs):
    """

    :param N:
    :param sigma2_n:
    :param proba_false_alarm:
    :param pd_ris_rx_chest:
    :param pd_probe_configs:
    :return:
    """

    # Compute measured power
    measured_power = np.abs(pd_ris_rx_chest)**2

    # Get best probing configurations according to highest received power
    best_indices = np.argmax(measured_power, axis=0)
    measured_power = measured_power.max(axis=0)

    # Compute threshold
    threshold = - (2 * N * sigma2_n) * np.log(proba_false_alarm)

    # Perform test
    detected_ues = measured_power > threshold

    # Update measured power 
    measured_power[~detected_ues] = 0.0

    # Compute weights
    weights = np.sqrt(measured_power) 
    weights = weights / np.sum(np.sqrt(measured_power), axis=-1)[:, None]

    # Extract relevant constants
    n_channels, K = detected_ues.shape

    # Prepare to store results
    hat_aoa = np.zeros((N, n_channels, K), dtype=np.complex_)
    reflection_configs = np.zeros((N, n_channels), dtype=np.complex_)
    
    hat_aoa[:] = np.nan + 1j * np.nan

    # Go through all channel realizations
    for rr in range(n_channels):

        # Go through all UEs
        for dd in range(K):

            # Check if UE was detected
            if detected_ues[rr, dd]:

                # Find the best configuration
                hat_aoa[:, rr, dd] = pd_probe_configs[best_indices[rr, dd]]

                # Compute reflection configuration
                reflection_configs[:, rr] += weights[rr, dd] * hat_aoa[:, rr, dd].conj()

    return reflection_configs, weights, hat_aoa


def sig_ris_probe(n_pilots, sigma2_n, proba_false_alarm, dsp_ris_rx_chest):
    """

    :param sigma2_n:
    :param proba_false_alarm:
    :param dsp_ris_rx_chest:
    :return:
    """

    # Extract useful constants
    n_probe_pilot_subblocks, N, n_channels, K = dsp_ris_rx_chest.shape

    # Average signal over probe subblocks
    avg_subblocks = np.mean(dsp_ris_rx_chest, axis=0)

    # Compute test 
    ynorms = np.linalg.norm(avg_subblocks, axis=0)
    test = 2 * N * n_probe_pilot_subblocks / n_pilots / sigma2_n * ynorms**2

    # Get threshold value
    threshold = chi2.ppf(1 - proba_false_alarm, df=2*N)

    # Perform test
    detected_ues = test > threshold

    # Update avg antennas
    avg_subblocks[:, ~detected_ues] = 0.0

    # Compute weights
    weights = ynorms
    weights = weights / weights.sum(axis=-1)[:, None]

    # Prepare to store the reflection configuration
    hat_aoa = np.zeros((N, n_channels, K), dtype=np.complex_)
    reflection_configs = np.zeros((N, n_channels), dtype=np.complex_)

    hat_aoa[:] = np.nan

    # Go through all channel realizations
    for rr in range(n_channels):

        # Go through all detected UEs
        for dd in range(K):

            # Check if UE was detected
            if detected_ues[rr, dd]:

                # Find the best configuration
                hat_aoa[:, rr, dd] = np.exp(1j * np.angle(avg_subblocks[:, rr, dd]))

                # Compute reflection configuration
                reflection_configs[:, rr] += weights[rr, dd] * hat_aoa[:, rr, dd].conj()

    return reflection_configs, weights, hat_aoa
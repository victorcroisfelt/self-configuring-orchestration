import numpy as np
from scipy.constants import speed_of_light

from tqdm import tqdm

from src.channel import generate_channel_realizations, scenario, drop_ues
from src.ris import pow_ris_config_codebook, ris_rx_chest, gen_ris_probe, pow_ris_probe, sig_ris_probe

# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    # Set random seed
    np.random.seed(42)

    ##################################################
    # BS Parameters
    ##################################################

    # Number of BS antennas
    M = 64

    ##################################################
    # HRIS Parameters
    ##################################################

    # Number of RIS elements
    N = 32

    ##################################################
    # UE Parameters
    ##################################################

    # Number of UEs
    K = 4

    # Transmit power at the UE = 10 dBm
    P_ue = 10 ** ((0 - 30) / 10)

    ##################################################
    # System Parameters
    ##################################################

    # Number of pilots
    n_pilots = K

    # Size of the coherence block
    size_coherence_block = 128

    # Number of pilot subblocks
    n_pilot_subblocks = 16

    ##################################################
    # Scenario Parameters
    ##################################################

    # Physical parameters
    freq = 28 * 10 ** 9
    wavelength = speed_of_light / freq

    # NLoS variances
    sigma2_dr = 0.1 * 6.1848 * 1e-12
    sigma2_rr = 0.1 * 5.9603 * 1e-4

    # Noise power
    sigma2_n_bs = 10 ** ((-94 - 30) / 10)
    sigma2_n_ris = 10 ** ((-91 - 30) / 10)

    # Generate scenario
    pos_bs, pos_bs_els, pos_ris, pos_ris_els, bs_ris_channels, ris_bs_steering, guard_distance_ris = scenario(wavelength, M, N)

    ##################################################
    # Simulation Parameters
    ##################################################

    # Define number of setups
    n_setups = 128

    # Define number of channel realizations
    n_channels = 64

    # Define number of noise realizations
    n_noise = 64

    # Number of pilot subbblocks for probe
    n_pilot_subblocks_probe_range = np.arange(1, n_pilot_subblocks + 1)
    n_xaxis = len(n_pilot_subblocks_probe_range)

    # Range of HRIS reflection parameter
    eta = 0.99

    # Define probability of false alarm
    proba_false_alarm = 0.001

    ##################################################
    # Simulation
    ##################################################

    # Prepare to save results
    pow_nmse1 = np.zeros((n_xaxis, n_setups, n_noise))
    sig_nmse1 = np.zeros((n_xaxis, n_setups, n_noise))

    pow_nmse1[:] = np.nan
    sig_nmse1[:] = np.nan

    pow_nmse2 = np.zeros((n_xaxis, n_setups, n_noise))
    sig_nmse2 = np.zeros((n_xaxis, n_setups, n_noise))

    pow_nmse2[:] = np.nan
    sig_nmse2[:] = np.nan

    pow_pd = np.zeros((n_xaxis, n_setups, n_noise))
    sig_pd = np.zeros((n_xaxis, n_setups, n_noise))

    pow_pd[:] = np.nan
    sig_pd[:] = np.nan

    # Go through all setups
    for ss in tqdm(range(n_setups)):

        # Drop the UEs over the area of interest
        pos_ues = drop_ues(K, pos_ris, dmax=100, guard_distance_ris=guard_distance_ris)

        # Generate UE channels
        bs_ue_channels, los_bs_ue_channels, ris_ue_channels, los_ris_ue_channels = generate_channel_realizations(
            wavelength, pos_bs, pos_bs_els, pos_ris, pos_ris_els, pos_ues, sigma2_dr, sigma2_rr, n_channels)

        # True reflection configuration
        gen_reflection_configs, gen_weights = gen_ris_probe(ris_ue_channels)

        # Go through all possible values for C
        for cc, n_pilot_subblocks_probe in enumerate(n_pilot_subblocks_probe_range):

            # Generate power-RIS configuration codebook
            pow_probe_configs = pow_ris_config_codebook(wavelength, n_pilot_subblocks_probe, pos_ris, pos_ris_els)

            # Go through noise realizations
            for nn in range(n_noise):

                # Compute received pilot signals
                pow_ris_rx_chest = ris_rx_chest(eta, P_ue, n_pilots, sigma2_n_ris, n_pilot_subblocks_probe, ris_ue_channels, pow_probe_configs)
                sig_ris_rx_chest = ris_rx_chest(eta, P_ue, n_pilots, sigma2_n_ris, n_pilot_subblocks_probe, ris_ue_channels)

                # HRIS probe
                pow_reflection_configs, pow_weights, pow_hat_aoa, pow_pd_ = pow_ris_probe(N, sigma2_n_ris, proba_false_alarm, pow_ris_rx_chest, pow_probe_configs)
                sig_reflection_configs, sig_weights, sig_hat_aoa, sig_pd_ = sig_ris_probe(n_pilots, sigma2_n_ris, proba_false_alarm, sig_ris_rx_chest)

                # Evaluate reflection configs
                pow_nmse_1 = np.linalg.norm(pow_reflection_configs - gen_reflection_configs,
                                           axis=0) ** 2 / np.linalg.norm(gen_reflection_configs, axis=0) ** 2

                sig_nmse_1 = np.linalg.norm(sig_reflection_configs - gen_reflection_configs,
                                           axis=0) ** 2 / np.linalg.norm(gen_reflection_configs, axis=0) ** 2

                # Compute true AoA information
                true_aoa = np.exp(1j * np.angle(ris_ue_channels))

                # Evaluate reflection configs
                pow_nmse_2 = np.linalg.norm(pow_hat_aoa - true_aoa,
                                           axis=0) ** 2 / np.linalg.norm(true_aoa, axis=0) ** 2

                sig_nmse_2 = np.linalg.norm(sig_hat_aoa - true_aoa,
                                           axis=0) ** 2 / np.linalg.norm(true_aoa, axis=0) ** 2

                # Store
                pow_nmse1[cc, ss, nn] = np.nanmean(pow_nmse_1)
                sig_nmse1[cc, ss, nn] = np.nanmean(sig_nmse_1)

                pow_nmse2[cc, ss, nn] = np.nanmean(pow_nmse_2)
                sig_nmse2[cc, ss, nn] = np.nanmean(sig_nmse_2)

                pow_pd[cc, ss, nn] = pow_pd_
                sig_pd[cc, ss, nn] = sig_pd_

    np.savez('data/figure6a.npz',
             n_pilot_subblocks_probe_range=n_pilot_subblocks_probe_range,
             pow_nmse1=pow_nmse1,
             pow_nmse2=pow_nmse2,
             pow_pd=pow_pd,
             sig_nmse1=sig_nmse1,
             sig_nmse2=sig_nmse2,
             sig_pd=sig_pd,
             )
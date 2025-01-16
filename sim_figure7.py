import numpy as np
from scipy.constants import speed_of_light

from tqdm import tqdm

from src.channel import generate_channel_realizations, scenario, drop_ues
from src.ris import pow_ris_config_codebook, ris_rx_chest_with_choice, gen_ris_probe, pow_ris_probe, sig_ris_probe

from scipy.stats.distributions import chi2

import matplotlib.pyplot as plt

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

    # Transmit power at the UE = 0 dBm
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
    sigma2_dr = 0.1 * 9.07 * 1e-9
    sigma2_rr = 0.1 * 1.12 * 1e-6

    # Noise power
    sigma2_n_bs = 10 ** ((-94 - 30) / 10)
    sigma2_n_ris = 10 ** ((-91 - 30) / 10)

    # Generate scenario
    pos_bs, pos_bs_els, pos_ris, pos_ris_els, bs_ris_channels, ris_bs_steering, guard_distance_ris = scenario(wavelength, M, N)

    ##################################################
    # Simulation Parameters
    ##################################################

    # Define number of setups
    n_setups = 100

    # Define number of channel realizations
    n_channels = 100

    # Define number of noise realizations
    n_noise = 100

    # Number of pilot subbblocks for probe
    n_pilot_subblocks_probe_range = np.arange(1, n_pilot_subblocks + 1)
    n_probe = len(n_pilot_subblocks_probe_range)

    # Range of HRIS reflection parameter
    eta_range = np.array([0.9, 0.99, 0.999, 0.9999])
    n_etas = len(eta_range)

    # Define probability of false alarm
    proba_false_alarm_range = np.array([0.01, 0.001])
    n_probas = len(proba_false_alarm_range)

    ##################################################
    # Simulation
    ##################################################

    # Prepare to save simulation results
    pow_test = np.zeros((n_etas, n_probe, n_noise, n_channels, n_setups))
    sig_test = np.zeros((n_etas, n_probe, n_noise, n_channels, n_setups))

    # Drop the UEs over the area of interest
    pos_ues = drop_ues(n_setups, pos_ris, dmax=1000, guard_distance_ris=900)

    # Generate UE channels
    bs_ue_channels, ris_ue_channels = generate_channel_realizations(wavelength, pos_bs, pos_bs_els, pos_ris, pos_ris_els, pos_ues, sigma2_dr, sigma2_rr, n_channels)
    
    # Generate UE's choice
    ue_choice = np.random.rand(n_setups)

    # Get mask
    mask = ue_choice > 0.5

    # Convert into True and False
    ue_choice[mask] = True
    ue_choice[~mask] = False

    #
    # Generate signals
    #

    # Go through all possible values for C
    for cc, n_pilot_subblocks_probe in tqdm(enumerate(n_pilot_subblocks_probe_range), total=n_pilot_subblocks):

        # Generate power-RIS configuration codebook
        pow_probe_configs = pow_ris_config_codebook(wavelength, n_pilot_subblocks_probe, pos_ris, pos_ris_els)

        # Go through all values of eta
        for ee in range(n_etas):

            # Get current value of eta
            eta = eta_range[ee]

            # Go through noise realizations
            for nn in range(n_noise):

                # Compute received pilot signals
                pow_ris_rx_chest = ris_rx_chest_with_choice(eta, P_ue, n_pilots, sigma2_n_ris, n_pilot_subblocks_probe, ris_ue_channels, mask, pow_probe_configs)
                sig_ris_rx_chest = ris_rx_chest_with_choice(eta, P_ue, n_pilots, sigma2_n_ris, n_pilot_subblocks_probe, ris_ue_channels, mask)
                
                # Power-based HRIS probe
                pow_test_curr = np.abs(pow_ris_rx_chest) ** 2
                pow_test_curr = np.max(pow_test_curr, axis=0)

                pow_test[ee, cc, nn, :, :] = pow_test_curr

                # Signal-based HRIS probe
                avg_subblocks = np.mean(sig_ris_rx_chest, axis=0)
                avg_antennas = np.mean(avg_subblocks, axis=0)
                sig_test_curr = np.abs(avg_antennas) ** 2
                sig_test_curr /= ((n_pilots * sigma2_n_ris) / (2 * N * n_pilot_subblocks_probe))

                sig_test[ee, cc, nn, :, :] = sig_test_curr

    print('time to test')

    #
    # Perform test
    #

    pow_threshold = np.zeros((n_probas))
    sig_threshold = np.zeros((n_probas))

    # Prepare to save test results
    pow_detection = np.zeros((n_probas, n_etas, n_probe, n_noise, n_channels))
    sig_detection = np.zeros((n_probas, n_etas, n_probe, n_noise, n_channels))

    # Prepare to save results
    pow_nmse = np.zeros((n_probas, n_etas, n_probe, n_noise, n_channels))
    sig_nmse = np.zeros((n_probas, n_etas, n_probe, n_noise, n_channels))
    
    pow_nmse[:] = np.nan
    sig_nmse[:] = np.nan

    # Genie reflection configuration
    gen_reflection_configs, gen_weights = gen_ris_probe(ris_ue_channels)

    # Go through all possible values for C
    for cc, n_pilot_subblocks_probe in tqdm(enumerate(n_pilot_subblocks_probe_range), total=n_pilot_subblocks):

        # Generate power-RIS configuration codebook
        pow_probe_configs = pow_ris_config_codebook(wavelength, n_pilot_subblocks_probe, pos_ris, pos_ris_els)

        # Go through all false alarm probabilities
        for pp in range(n_probas):

            # Get current probability of false alarm
            proba_false_alarm = proba_false_alarm_range[pp]

            # Get threshold values
            pow_threshold[pp] = - (2 * N * sigma2_n_ris) * np.log(proba_false_alarm)
            sig_threshold[pp] = chi2.ppf(1 - proba_false_alarm, df=4)

            # Go through all values of eta
            for ee in range(n_etas):

                # Go through the channels
                for ch in range(n_channels):

                    # Go through noise realizations
                    for nn in range(n_noise):

                        # Perform test
                        pow_mask = pow_test[ee, cc, nn, ch, :] > pow_threshold[pp]
                        sig_mask = sig_test[ee, cc, nn, ch, :] > sig_threshold[pp]
    
                        pow_detection[pp, ee, cc, nn, ch] = np.mean(pow_mask == mask)
                        sig_detection[pp, ee, cc, nn, ch] = np.mean(sig_mask == mask)

                        #####
                        pow_ris_rx_chest = ris_rx_chest_with_choice(eta, P_ue, n_pilots, sigma2_n_ris, n_pilot_subblocks_probe, ris_ue_channels, pow_mask, pow_probe_configs)
                        sig_ris_rx_chest = ris_rx_chest_with_choice(eta, P_ue, n_pilots, sigma2_n_ris, n_pilot_subblocks_probe, ris_ue_channels, sig_mask)

                        #####
                        pow_reflection_configs, pow_weights, pow_hat_aoa, pow_pd_ = pow_ris_probe(N, sigma2_n_ris, proba_false_alarm, pow_ris_rx_chest, pow_probe_configs)
                        sig_reflection_configs, sig_weights, sig_hat_aoa, sig_pd_ = sig_ris_probe(n_pilots, sigma2_n_ris, proba_false_alarm, sig_ris_rx_chest)

                        # Evaluate reflection configs
                        pow_nmse = np.linalg.norm(pow_reflection_configs - gen_reflection_configs,
                                                axis=0) ** 2 / np.linalg.norm(gen_reflection_configs, axis=0) ** 2

                        sig_nmse = np.linalg.norm(sig_reflection_configs - gen_reflection_configs,
                                                axis=0) ** 2 / np.linalg.norm(gen_reflection_configs, axis=0) ** 2

    np.savez('data/new_figure6a_xx.npz',
             pow_test=pow_test,
             sig_test=sig_test,
             pow_threshold=pow_threshold,
             sig_threshold=sig_threshold,
             pow_detection=pow_detection,
             sig_detection=sig_detection,
             pow_nmse=pow_nmse,
             sig_nmse=sig_nmse
             )

#     # Store values
#     pow_pd[pp, ee, cc, ss, nn] = pow_pd_
#     sig_pd[pp, ee, cc, ss, nn] = sig_pd_
# # Store
# pow_nmse1[cc, ss, nn] = np.nanmean(pow_nmse_1)
# sig_nmse1[cc, ss, nn] = np.nanmean(sig_nmse_1)

# pow_nmse2[cc, ss, nn] = np.nanmean(pow_nmse_2)
# sig_nmse2[cc, ss, nn] = np.nanmean(sig_nmse_2)
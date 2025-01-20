import numpy as np
from scipy.constants import speed_of_light

from tqdm import tqdm

from src.channel import generate_channel_realizations, scenario, drop_ues
from src.ris import pow_ris_config_codebook, ris_rx_chest, gen_ris_probe, pow_ris_probe, sig_ris_probe
from src.mmimo import bs_rx_chest, bs_comm, bs_rx_chest_no_probe

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

    # Coherence interval length
    coherence_interval_length = 128

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

    # Define probability of false alarm
    proba_false_alarm = 0.01

    # HRIS reflection parameter
    eta = 0.9999

    # Number of pilots
    n_pilots = K

    # Number of pilot subblocks
    n_pilot_subblocks = int(64 // K)

    # Number of probe pilot subbblocks
    n_pilot_subblocks_probe_range = np.arange(1, n_pilot_subblocks + 1)

    # Calculate pre-log term
    pre_log_term = (coherence_interval_length - n_pilot_subblocks * n_pilots) / coherence_interval_length

    ##################################################
    # Simulation
    ##################################################

    # Prepare to save results
    gen_avg_nmse = np.zeros((n_pilot_subblocks, n_setups, n_noise, n_channels, K))
    pow_avg_nmse = np.zeros((n_pilot_subblocks, n_setups, n_noise, n_channels, K))
    sig_avg_nmse = np.zeros((n_pilot_subblocks, n_setups, n_noise, n_channels, K))

    gen_avg_se = np.zeros((2, n_pilot_subblocks, n_setups, n_noise, n_channels, K))
    pow_avg_se = np.zeros((2, n_pilot_subblocks, n_setups, n_noise, n_channels, K))
    sig_avg_se = np.zeros((2, n_pilot_subblocks, n_setups, n_noise, n_channels, K))

    gen_avg_se_pos = np.zeros((2, n_pilot_subblocks, n_setups, n_noise, n_channels, K))
    pow_avg_se_pos = np.zeros((2, n_pilot_subblocks, n_setups, n_noise, n_channels, K))
    sig_avg_se_pos = np.zeros((2, n_pilot_subblocks, n_setups, n_noise, n_channels, K))

    gen_avg_nmse[:] = np.nan
    pow_avg_nmse[:] = np.nan
    sig_avg_nmse[:] = np.nan

    gen_avg_se[:] = np.nan
    pow_avg_se[:] = np.nan
    sig_avg_se[:] = np.nan

    gen_avg_se_pos[:] = np.nan
    pow_avg_se_pos[:] = np.nan
    sig_avg_se_pos[:] = np.nan

    # Go through all setups
    for ss in tqdm(range(n_setups), desc='setups'):

        # Drop the UEs over the area of interest
        pos_ues = drop_ues(K, pos_ris, dmax=1000, guard_distance_ris=900)

        # Generate UE channels
        bs_ue_channels, ris_ue_channels = generate_channel_realizations(
            wavelength, pos_bs, pos_bs_els, pos_ris, pos_ris_els, pos_ues, sigma2_dr, sigma2_rr, n_channels)
        
        # Genie reflection configuration
        gen_reflection_configs, gen_weights = gen_ris_probe(ris_ue_channels)

        # Go through all points in the x-dimension
        for cc, n_pilot_subblocks_probe in enumerate(n_pilot_subblocks_probe_range):
            
            if n_pilot_subblocks_probe == n_pilot_subblocks:
                break

            # Generate power-RIS configuration codebook
            pow_probe_configs = pow_ris_config_codebook(wavelength, n_pilot_subblocks_probe, pos_ris, pos_ris_els)

            # Go through noise realizations
            for nn in range(n_noise):

                # Compute received pilot signals
                pow_ris_rx_chest = ris_rx_chest(eta, P_ue, n_pilots, sigma2_n_ris, n_pilot_subblocks_probe, ris_ue_channels, pow_probe_configs)
                sig_ris_rx_chest = ris_rx_chest(eta, P_ue, n_pilots, sigma2_n_ris, n_pilot_subblocks_probe, ris_ue_channels)

                # HRIS probe
                pow_reflection_configs, pow_weights, pow_hat_aoa = pow_ris_probe(N, sigma2_n_ris, proba_false_alarm, pow_ris_rx_chest, pow_probe_configs)
                sig_reflection_configs, sig_weights, sig_hat_aoa = sig_ris_probe(n_pilots, sigma2_n_ris, proba_false_alarm, sig_ris_rx_chest)

                # Complete reflection configurations by inserting BS-RIS knowledge
                gen_reflection_configs *= ris_bs_steering[:, None]
                pow_reflection_configs *= ris_bs_steering[:, None]
                sig_reflection_configs *= ris_bs_steering[:, None]

                # Compute reflected channels during probe
                pow_refl_channels_probe = pow_probe_configs[:, :, None, None] * ris_ue_channels[None, :, :, :]
                pow_refl_channels_probe = bs_ris_channels[None, :, :, None, None] * pow_refl_channels_probe[:, None, :, :, :]
                pow_refl_channels_probe = pow_refl_channels_probe.sum(axis=2)
                pow_refl_channels_probe = pow_refl_channels_probe.sum(axis=0)

                ones = np.ones_like(pow_probe_configs)
                sig_refl_channels_probe = ones[:, :, None, None] * ris_ue_channels[None, :, :, :]
                sig_refl_channels_probe = bs_ris_channels[None, :, :, None, None] * sig_refl_channels_probe[:, None, :, :, :]
                sig_refl_channels_probe = sig_refl_channels_probe.sum(axis=2)
                sig_refl_channels_probe = sig_refl_channels_probe.sum(axis=0)

                # Compute equivalent channels during probe
                pow_eq_channels_probe = bs_ue_channels + np.sqrt(eta) * pow_refl_channels_probe
                sig_eq_channels_probe = bs_ue_channels + np.sqrt(eta) * sig_refl_channels_probe

                # Compute reflected channels during communication
                gen_refl_channels = gen_reflection_configs[:, :, None] * ris_ue_channels
                gen_refl_channels = bs_ris_channels[:, :, None, None] * gen_refl_channels[None, :, :, :]
                gen_refl_channels = gen_refl_channels.sum(axis=1)

                pow_refl_channels = pow_reflection_configs[:, :, None] * ris_ue_channels
                pow_refl_channels = bs_ris_channels[:, :, None, None] * pow_refl_channels[None, :, :, :]
                pow_refl_channels = pow_refl_channels.sum(axis=1)

                sig_refl_channels = sig_reflection_configs[:, :, None] * ris_ue_channels
                sig_refl_channels = bs_ris_channels[:, :, None, None] * sig_refl_channels[None, :, :, :]
                sig_refl_channels = sig_refl_channels.sum(axis=1)

                # Compute equivalent channels during communication
                gen_eq_channels = bs_ue_channels + np.sqrt(eta) * gen_refl_channels
                pow_eq_channels = bs_ue_channels + np.sqrt(eta) * pow_refl_channels
                sig_eq_channels = bs_ue_channels + np.sqrt(eta) * sig_refl_channels

                # Get channel estimates
                gen_hat_eq_channels = bs_rx_chest_no_probe(P_ue, n_pilots, sigma2_n_bs, n_pilot_subblocks,
                                                  n_pilot_subblocks_probe, gen_eq_channels)

                pow_hat_eq_channels = bs_rx_chest(P_ue, n_pilots, sigma2_n_bs, n_pilot_subblocks,
                                                  n_pilot_subblocks_probe, pow_eq_channels_probe, pow_eq_channels)

                sig_hat_eq_channels = bs_rx_chest(P_ue, n_pilots, sigma2_n_bs, n_pilot_subblocks,
                                                  n_pilot_subblocks_probe, sig_eq_channels_probe, sig_eq_channels)

                # Compute normalized mean squared error
                diff = gen_hat_eq_channels - gen_eq_channels
                gen_avg_nmse[cc, ss, nn, :, :] = (
                            np.linalg.norm(diff, axis=0) ** 2 / np.linalg.norm(gen_eq_channels, axis=0))

                diff = pow_hat_eq_channels - pow_eq_channels
                pow_avg_nmse[cc, ss, nn, :, :] = (
                            np.linalg.norm(diff, axis=0) ** 2 / np.linalg.norm(pow_eq_channels, axis=0))

                diff = sig_hat_eq_channels - sig_eq_channels
                sig_avg_nmse[cc, ss, nn, :, :] = (
                            np.linalg.norm(diff, axis=0) ** 2 / np.linalg.norm(sig_eq_channels, axis=0))

                ##################################################
                # Communication Phase
                ##################################################
                gen_se, _, _, _ = bs_comm(P_ue, sigma2_n_bs, gen_eq_channels, gen_hat_eq_channels)
                pow_se, _, _, _ = bs_comm(P_ue, sigma2_n_bs, pow_eq_channels, pow_hat_eq_channels)
                sig_se, _, _, _ = bs_comm(P_ue, sigma2_n_bs, sig_eq_channels, sig_hat_eq_channels)

                # Store results for MR
                gen_avg_se[0, cc, ss, nn, :, :] = gen_se
                pow_avg_se[0, cc, ss, nn, :, :] = pow_se
                sig_avg_se[0, cc, ss, nn, :, :] = sig_se

                gen_avg_se_pos[0, cc, ss, nn, :, :] = pre_log_term * gen_se
                pow_avg_se_pos[0, cc, ss, nn, :, :] = pre_log_term * pow_se
                sig_avg_se_pos[0, cc, ss, nn, :, :] = pre_log_term * sig_se

                gen_se, _, _, _ = bs_comm(P_ue, sigma2_n_bs, gen_eq_channels, gen_hat_eq_channels, method='ZF')
                pow_se, _, _, _ = bs_comm(P_ue, sigma2_n_bs, pow_eq_channels, pow_hat_eq_channels, method='ZF')
                sig_se, _, _, _ = bs_comm(P_ue, sigma2_n_bs, sig_eq_channels, sig_hat_eq_channels, method='ZF')

                # Store results for ZF
                gen_avg_se[1, cc, ss, nn, :, :] = gen_se
                pow_avg_se[1, cc, ss, nn, :, :] = pow_se
                sig_avg_se[1, cc, ss, nn, :, :] = sig_se

                gen_avg_se_pos[1, cc, ss, nn, :, :] = pre_log_term * gen_se
                pow_avg_se_pos[1, cc, ss, nn, :, :] = pre_log_term * pow_se
                sig_avg_se_pos[1, cc, ss, nn, :, :] = pre_log_term * sig_se

    np.savez('data/figure9_gen-ris_K' + str(K) + '_N' + str(N) + '.npz',
            n_pilot_subblocks=n_pilot_subblocks, 
            n_pilot_subblocks_probe_range=n_pilot_subblocks_probe_range,
            gen_avg_nmse=gen_avg_nmse,
            gen_avg_se=gen_avg_se,
            gen_avg_se_pos=gen_avg_se_pos
            )

    np.savez('data/figure9_pow-ris_K' + str(K) + '_N' + str(N) + '.npz',
            n_pilot_subblocks=n_pilot_subblocks, 
            n_pilot_subblocks_probe_range=n_pilot_subblocks_probe_range,
            pow_avg_nmse=pow_avg_nmse,
            pow_avg_se=pow_avg_se,
            pow_avg_se_pos=pow_avg_se_pos
            )

    np.savez('data/figure9_sig-ris_K' + str(K) + '_N' + str(N) + '.npz',
            n_pilot_subblocks=n_pilot_subblocks, 
            n_pilot_subblocks_probe_range=n_pilot_subblocks_probe_range,
            sig_avg_nmse=sig_avg_nmse,
            sig_avg_se=sig_avg_se,
            sig_avg_se_pos=sig_avg_se_pos
            )
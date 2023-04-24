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
    tau_c = 128

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

    # Maximum distance
    distance_max = 100

    ##################################################
    # Simulation Parameters
    ##################################################

    # Define number of setups
    n_setups = 128

    # Define number of channel realizations
    n_channels = 64

    # Define number of noise realizations
    n_noise = 64

    # Define probability of false alarm
    proba_false_alarm = 0.1

    # HRIS reflection parameter
    eta = 0.9

    # Krange
    K_range = np.array([1, 2, 4, 8, 16])
    n_yaxis = len(K_range)

    #
    n_xaxis = 64

    ##################################################
    # Simulation
    ##################################################

    # Prepare to save results
    gen_avg_nmse = np.zeros((n_yaxis, n_xaxis, n_setups, n_channels, n_noise))
    pow_avg_nmse = np.zeros((n_yaxis, n_xaxis, n_setups, n_channels, n_noise))
    sig_avg_nmse = np.zeros((n_yaxis, n_xaxis, n_setups, n_channels, n_noise))

    pow_avg_pd = np.zeros((n_yaxis, n_xaxis, n_setups, n_channels, n_noise))
    sig_avg_pd = np.zeros((n_yaxis, n_xaxis, n_setups, n_channels, n_noise))

    gen_avg_se = np.zeros((2, n_yaxis, n_xaxis, n_setups, n_channels, n_noise))
    pow_avg_se = np.zeros((2, n_yaxis, n_xaxis, n_setups, n_channels, n_noise))
    sig_avg_se = np.zeros((2, n_yaxis, n_xaxis, n_setups, n_channels, n_noise))

    gen_avg_num = np.zeros((2, n_yaxis, n_xaxis, n_setups, n_channels, n_noise))
    gen_avg_den1 = np.zeros((2, n_yaxis, n_xaxis, n_setups, n_channels, n_noise))
    gen_avg_den2 = np.zeros((2, n_yaxis, n_xaxis, n_setups, n_channels, n_noise))

    pow_avg_num = np.zeros((2, n_yaxis, n_xaxis, n_setups, n_channels, n_noise))
    pow_avg_den1 = np.zeros((2, n_yaxis, n_xaxis, n_setups, n_channels, n_noise))
    pow_avg_den2 = np.zeros((2, n_yaxis, n_xaxis, n_setups, n_channels, n_noise))

    sig_avg_num = np.zeros((2, n_yaxis, n_xaxis, n_setups, n_channels, n_noise))
    sig_avg_den1 = np.zeros((2, n_yaxis, n_xaxis, n_setups, n_channels, n_noise))
    sig_avg_den2 = np.zeros((2, n_yaxis, n_xaxis, n_setups, n_channels, n_noise))

    gen_avg_sir = np.zeros((2, n_yaxis, n_xaxis, n_setups, n_channels, n_noise))
    pow_avg_sir = np.zeros((2, n_yaxis, n_xaxis, n_setups, n_channels, n_noise))
    sig_avg_sir = np.zeros((2, n_yaxis, n_xaxis, n_setups, n_channels, n_noise))

    gen_avg_nmse[:] = np.nan
    pow_avg_nmse[:] = np.nan
    sig_avg_nmse[:] = np.nan

    pow_avg_pd[:] = np.nan
    sig_avg_pd[:] = np.nan

    gen_avg_se[:] = np.nan
    pow_avg_se[:] = np.nan
    sig_avg_se[:] = np.nan

    gen_avg_num[:] = np.nan
    gen_avg_den1[:] = np.nan
    gen_avg_den2[:] = np.nan

    pow_avg_num[:] = np.nan
    pow_avg_den1[:] = np.nan
    pow_avg_den2[:] = np.nan

    sig_avg_num[:] = np.nan
    sig_avg_den1[:] = np.nan
    sig_avg_den2[:] = np.nan

    gen_avg_sir[:] = np.nan
    pow_avg_sir[:] = np.nan
    sig_avg_sir[:] = np.nan

    # Go through all setups
    for ss in tqdm(range(n_setups)):

        # Go through all number of UEs
        for kk, K in enumerate(K_range):

            K = int(K)

            # Number of pilots
            n_pilots = K

            # Number of pilot subblocks
            n_pilot_subblocks = int(64 // K)

            # Number of probe pilot subbblocks
            n_pilot_subblocks_probe_range = np.arange(1, n_pilot_subblocks + 1)

            # Calculate pre-log term
            pre_log_term = (tau_c - n_pilot_subblocks * n_pilots) / tau_c

            # Drop the UEs over the area of interest
            pos_ues = drop_ues(K, pos_ris, dmax=distance_max, guard_distance_ris=guard_distance_ris)

            # Generate UE channels
            bs_ue_channels, los_bs_ue_channels, ris_ue_channels, los_ris_ue_channels = generate_channel_realizations(
                wavelength, pos_bs, pos_bs_els, pos_ris, pos_ris_els, pos_ues, sigma2_dr, sigma2_rr, n_channels)

            # Genie reflection configuration
            gen_reflection_configs, gen_weights = gen_ris_probe(ris_ue_channels)

            # Go through all points in the x-dimension
            for cc, n_pilot_subblocks_probe in enumerate(n_pilot_subblocks_probe_range):

                # Generate power-RIS configuration codebook
                pow_probe_configs = pow_ris_config_codebook(wavelength, n_pilot_subblocks_probe, pos_ris, pos_ris_els)

                # Go through noise realizations
                for nn in range(n_noise):

                    # Compute received pilot signals
                    pow_ris_rx_chest = ris_rx_chest(eta, P_ue, n_pilots, sigma2_n_ris, n_pilot_subblocks_probe, ris_ue_channels, pow_probe_configs)
                    sig_ris_rx_chest = ris_rx_chest(eta, P_ue, n_pilots, sigma2_n_ris, n_pilot_subblocks_probe, ris_ue_channels)

                    # HRIS probe
                    pow_reflection_configs, pow_weights, pow_hat_aoa, pow_avg_pd[kk, cc, ss, :, nn] = pow_ris_probe(N, sigma2_n_ris, proba_false_alarm, pow_ris_rx_chest, pow_probe_configs)
                    sig_reflection_configs, sig_weights, sig_hat_aoa, pow_avg_pd[kk, cc, ss, :, nn] = sig_ris_probe(n_pilots, sigma2_n_ris, proba_false_alarm, sig_ris_rx_chest)

                    # Complete reflection configurations by inserting BS-RIS knowledge
                    gen_reflection_configs *= ris_bs_steering[:, None]
                    pow_reflection_configs *= ris_bs_steering[:, None]
                    sig_reflection_configs *= ris_bs_steering[:, None]

                    # Compute reflected channels during probe
                    pow_refl_channels_probe = pow_probe_configs[:, None, :, None] * ris_ue_channels[None, :, :, :]
                    pow_refl_channels_probe = bs_ris_channels[None, None, :, :, None] * pow_refl_channels_probe[:, :, None,
                                                                                        :, :]
                    pow_refl_channels_probe = pow_refl_channels_probe.sum(axis=3)
                    pow_refl_channels_probe = pow_refl_channels_probe.sum(axis=0)

                    ones = np.ones_like(pow_probe_configs)
                    sig_refl_channels_probe = ones[:, None, :, None] * ris_ue_channels[None, :, :, :]
                    sig_refl_channels_probe = bs_ris_channels[None, None, :, :, None] * sig_refl_channels_probe[:, :, None,
                                                                                        :, :]
                    sig_refl_channels_probe = sig_refl_channels_probe.sum(axis=3)
                    sig_refl_channels_probe = sig_refl_channels_probe.sum(axis=0)

                    # Compute equivalent channels during probe
                    pow_eq_channels_probe = bs_ue_channels + np.sqrt(eta) * pow_refl_channels_probe
                    sig_eq_channels_probe = bs_ue_channels + np.sqrt(eta) * sig_refl_channels_probe

                    # Compute reflected channels during communication
                    gen_refl_channels = ((bs_ris_channels[:, :, None] * gen_reflection_configs[None, :, :])[None, :, :, :] *
                                         ris_ue_channels[:, None, :, :]).sum(axis=2)

                    pow_refl_channels = ((bs_ris_channels[:, :, None] * pow_reflection_configs[None, :, :])[None, :, :, :] *
                                         ris_ue_channels[:, None, :, :]).sum(axis=2)

                    sig_refl_channels = ((bs_ris_channels[:, :, None] * sig_reflection_configs[None, :, :])[None, :, :, :] *
                                         ris_ue_channels[:, None, :, :]).sum(axis=2)

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
                    gen_avg_nmse[kk, cc, ss, :, nn] = (
                                np.linalg.norm(diff, axis=1) ** 2 / np.linalg.norm(gen_eq_channels, axis=1)).mean(axis=0)

                    diff = pow_hat_eq_channels - pow_eq_channels
                    pow_avg_nmse[kk, cc, ss, :, nn] = (
                                np.linalg.norm(diff, axis=1) ** 2 / np.linalg.norm(pow_eq_channels, axis=1)).mean(axis=0)

                    diff = sig_hat_eq_channels - sig_eq_channels
                    sig_avg_nmse[kk, cc, ss, :, nn] = (
                                np.linalg.norm(diff, axis=1) ** 2 / np.linalg.norm(sig_eq_channels, axis=1)).mean(axis=0)

                    ##################################################
                    # Communication Phase
                    ##################################################
                    gen_se, gen_num, gen_den1, gen_den2 = bs_comm(P_ue, sigma2_n_bs, gen_eq_channels, gen_hat_eq_channels)
                    pow_se, pow_num, pow_den1, pow_den2 = bs_comm(P_ue, sigma2_n_bs, pow_eq_channels, pow_hat_eq_channels)
                    sig_se, sig_num, sig_den1, sig_den2 = bs_comm(P_ue, sigma2_n_bs, sig_eq_channels, sig_hat_eq_channels)

                    # Store results
                    gen_avg_se[0, kk, cc, ss, :, nn] = pre_log_term * gen_se.sum(axis=0)
                    pow_avg_se[0, kk, cc, ss, :, nn] = pre_log_term * pow_se.sum(axis=0)
                    sig_avg_se[0, kk, cc, ss, :, nn] = pre_log_term * sig_se.sum(axis=0)

                    gen_avg_num[0, kk, cc, ss, :, nn] = gen_num.mean(axis=0)
                    pow_avg_num[0, kk, cc, ss, :, nn] = pow_num.mean(axis=0)
                    sig_avg_num[0, kk, cc, ss, :, nn] = sig_num.mean(axis=0)

                    gen_avg_den1[0, kk, cc, ss, :, nn] = gen_den1.mean(axis=0)
                    pow_avg_den1[0, kk, cc, ss, :, nn] = pow_den1.mean(axis=0)
                    sig_avg_den1[0, kk, cc, ss, :, nn] = sig_den1.mean(axis=0)

                    gen_avg_den2[0, kk, cc, ss, :, nn] = gen_den2.mean(axis=0)
                    pow_avg_den2[0, kk, cc, ss, :, nn] = pow_den2.mean(axis=0)
                    sig_avg_den2[0, kk, cc, ss, :, nn] = sig_den2.mean(axis=0)

                    gen_avg_sir[0, kk, cc, ss, :, nn] = (gen_num / gen_den1).mean(axis=0)
                    pow_avg_sir[0, kk, cc, ss, :, nn] = (pow_num / pow_den1).mean(axis=0)
                    sig_avg_sir[0, kk, cc, ss, :, nn] = (sig_num / sig_den1).mean(axis=0)

                    gen_se, gen_num, gen_den1, gen_den2 = bs_comm(P_ue, sigma2_n_bs, gen_eq_channels, gen_hat_eq_channels, method='ZF')
                    pow_se, pow_num, pow_den1, pow_den2 = bs_comm(P_ue, sigma2_n_bs, pow_eq_channels, pow_hat_eq_channels, method='ZF')
                    sig_se, sig_num, sig_den1, sig_den2 = bs_comm(P_ue, sigma2_n_bs, sig_eq_channels, sig_hat_eq_channels, method='ZF')

                    # Store results
                    gen_avg_se[1, kk, cc, ss, :, nn] = pre_log_term * gen_se.sum(axis=0)
                    pow_avg_se[1, kk, cc, ss, :, nn] = pre_log_term * pow_se.sum(axis=0)
                    sig_avg_se[1, kk, cc, ss, :, nn] = pre_log_term * sig_se.sum(axis=0)

                    gen_avg_num[1, kk, cc, ss, :, nn] = gen_num.mean(axis=0)
                    pow_avg_num[1, kk, cc, ss, :, nn] = pow_num.mean(axis=0)
                    sig_avg_num[1, kk, cc, ss, :, nn] = sig_num.mean(axis=0)

                    gen_avg_den1[1, kk, cc, ss, :, nn] = gen_den1.mean(axis=0)
                    pow_avg_den1[1, kk, cc, ss, :, nn] = pow_den1.mean(axis=0)
                    sig_avg_den1[1, kk, cc, ss, :, nn] = sig_den1.mean(axis=0)

                    gen_avg_den2[1, kk, cc, ss, :, nn] = gen_den2.mean(axis=0)
                    pow_avg_den2[1, kk, cc, ss, :, nn] = pow_den2.mean(axis=0)
                    sig_avg_den2[1, kk, cc, ss, :, nn] = sig_den2.mean(axis=0)

                    gen_avg_sir[1, kk, cc, ss, :, nn] = (gen_num / gen_den1).mean(axis=0)
                    pow_avg_sir[1, kk, cc, ss, :, nn] = (pow_num / pow_den1).mean(axis=0)
                    sig_avg_sir[1, kk, cc, ss, :, nn] = (sig_num / sig_den1).mean(axis=0)

    np.savez('data/figure7_gen-ris' + str(K) + '_N' + str(N) + '.npz',
             K_range=K_range,
             gen_avg_nmse=gen_avg_nmse,
             gen_avg_se=gen_avg_se,
             gen_avg_num=gen_avg_num,
             gen_avg_den1=gen_avg_den1,
             gen_avg_den2=gen_avg_den2,
             gen_avg_sir=gen_avg_sir
             )

    np.savez('data/figure7_pow-ris' + str(K) + '_N' + str(N) + '.npz',
             K_range=K_range,
             pow_avg_nmse=pow_avg_nmse,
             pow_avg_pd=pow_avg_pd,
             pow_avg_se=pow_avg_se,
             pow_avg_num=pow_avg_num,
             pow_avg_den1=pow_avg_den1,
             pow_avg_den2=pow_avg_den2,
             pow_avg_sir=pow_avg_sir
             )

    np.savez('data/figure7_sig-ris_K' + str(K) + '_N' + str(N) + '.npz',
             K_range=K_range,
             sig_avg_nmse=sig_avg_nmse,
             sig_avg_pd=sig_avg_pd,
             sig_avg_se=sig_avg_se,
             sig_avg_num=sig_avg_num,
             sig_avg_den1=sig_avg_den1,
             sig_avg_den2=sig_avg_den2,
             sig_avg_sir=sig_avg_sir
             )

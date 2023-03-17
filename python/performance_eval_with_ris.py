import numpy as np
from scipy.constants import speed_of_light

import matplotlib.pyplot as plt

from tqdm import tqdm

from src.channel import generate_channel_realizations, scenario, drop_ues
from src.ris import pow_ris_config_codebook, ris_rx_chest, pow_ris_probe, sig_ris_probe
from src.mmimo import bs_rx_chest, bs_comm

from matplotlib import rc

rc('font', **{'family': 'sans serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press ⌘F8 to toggle the breakpoint.


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

    # HRIS reflection parameter
    eta = 0.8

    ##################################################
    # UE Parameters
    ##################################################

    # Number of UEs
    K = 4

    # Transmit power at the UE = 10 dBm
    P_ue = 10 ** ((10 - 30) / 10)

    ##################################################
    # System Parameters
    ##################################################

    # Number of pilots
    n_pilots = 8

    # Number of pilot subblocks
    n_pilot_subblocks = 16

    # Channel estimation length
    tau_chest = n_pilot_subblocks * n_pilots

    # Coherence interval length
    tau_c = 2 * tau_chest

    # Communication length
    tau_comm = tau_chest

    ##################################################
    # Scenario Parameters
    ##################################################

    # Physical parameters
    freq = 3 * 10 ** 9
    wavelength = speed_of_light / freq

    # NLoS variances
    sigma2_dr = 1e-10
    sigma2_rr = 1e-6

    # Noise power
    sigma2_n = 10 ** ((-94 - 30) / 10)

    # Parameterized variable to define area of interest
    scenario_size = 100

    # Guard distance
    guard_distance = 30

    # Generate scenario
    pos_bs, pos_bs_els, pos_ris, pos_ris_els, bs_ris_channels, ris_bs_steering = scenario(wavelength, scenario_size, M, N)

    ##################################################
    # Simulation Parameters
    ##################################################

    # Define number of setups
    n_setups = 10

    # Define number of channel realizations
    n_channels = 100

    # Define number of noise realizations
    n_noise = 100

    # Number of probe pilot subbblocks
    n_pilot_subblocks_probe_range = np.arange(1, n_pilot_subblocks + 1)
    #n_pilot_subblocks_probe_range = np.array([8])
    n_xaxis = len(n_pilot_subblocks_probe_range)

    # Define probability of false alarm
    proba_false_alarm = 0.1

    ##################################################
    # Simulation
    ##################################################

    # Prepare to save results
    pow_avg_nmse = np.zeros((n_xaxis, n_setups, n_channels, n_noise))
    sig_avg_nmse = np.zeros((n_xaxis, n_setups, n_channels, n_noise))

    pow_avg_se = np.zeros((n_xaxis, n_setups, n_channels, n_noise))
    sig_avg_se = np.zeros((n_xaxis, n_setups, n_channels, n_noise))

    pow_avg_num = np.zeros((n_xaxis, n_setups, n_channels, n_noise))
    pow_avg_den1 = np.zeros((n_xaxis, n_setups, n_channels, n_noise))
    pow_avg_den2 = np.zeros((n_xaxis, n_setups, n_channels, n_noise))

    sig_avg_num = np.zeros((n_xaxis, n_setups, n_channels, n_noise))
    sig_avg_den1 = np.zeros((n_xaxis, n_setups, n_channels, n_noise))
    sig_avg_den2 = np.zeros((n_xaxis, n_setups, n_channels, n_noise))

    # Go through all setups
    for ss in tqdm(range(n_setups)):

        # Drop the UEs over the area of interest
        pos_ues = drop_ues(K, pos_ris, scenario_size, guard_distance)

        # fig, ax = plt.subplots()
        #
        # ax.plot(pos_bs[0], pos_bs[1], 's', label='BS')
        # ax.plot(pos_ris[0], pos_ris[1], 'x', label='HRIS')
        # ax.plot(pos_ues[:, 0], pos_ues[:, 1], '.', label='UE')
        #
        # ax.legend()
        #
        # ax.set_xlabel('$x$ [m]')
        # ax.set_ylabel('$y$ [m]')
        #
        # plt.tight_layout()
        #
        # plt.show()
        #breakpoint()

        # Generate UE channels
        bs_ue_channels, los_bs_ue_channels, ris_ue_channels, los_ris_ue_channels = generate_channel_realizations(
            wavelength, pos_bs, pos_bs_els, pos_ris, pos_ris_els, pos_ues, sigma2_dr, sigma2_rr, n_channels)

        # Go through all points in the x-dimension
        for cc, n_pilot_subblocks_probe in enumerate(n_pilot_subblocks_probe_range):

            # Generate power-RIS configuration codebook
            pow_probe_configs = pow_ris_config_codebook(wavelength, n_pilot_subblocks_probe, pos_ris, pos_ris_els)

            # Go through noise realizations
            for nn in range(n_noise):

                ##################################################
                # Chest Phase
                ##################################################

                # Compute received signals at the power-RIS
                pow_ris_rx_chest = ris_rx_chest(eta, P_ue, n_pilots, sigma2_n, n_pilot_subblocks_probe, ris_ue_channels,
                                                pow_probe_configs)
                sig_ris_rx_chest = ris_rx_chest(eta, P_ue, n_pilots, sigma2_n, n_pilot_subblocks_probe, ris_ue_channels)

                # HRIS probe
                pow_reflection_configs = pow_ris_probe(N, sigma2_n, proba_false_alarm, pow_ris_rx_chest,
                                                       pow_probe_configs)
                sig_reflection_configs = sig_ris_probe(n_pilots, sigma2_n, proba_false_alarm, sig_ris_rx_chest)

                #####

                # # Print angles
                # print('los = ', np.angle(los_ris_ue_channels.squeeze()))
                # print('realization = ', np.angle(ris_ue_channels[:, :, 0].squeeze()))
                #
                # print('pow = ', np.angle(pow_reflection_configs[:, 0].squeeze()))
                # print('sig = ', np.angle(sig_reflection_configs[:, 0].squeeze()))

                #breakpoint()

                # Complete reflection configurations by inserting BS-RIS knowledge
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
                pow_refl_channels = ((bs_ris_channels[:, :, None] * pow_reflection_configs[None, :, :])[None, :, :, :] *
                                     ris_ue_channels[:, None, :, :]).sum(axis=2)

                sig_refl_channels = ((bs_ris_channels[:, :, None] * sig_reflection_configs[None, :, :])[None, :, :, :] *
                                     ris_ue_channels[:, None, :, :]).sum(axis=2)

                # Compute equivalent channels during communication
                pow_eq_channels = bs_ue_channels + np.sqrt(eta) * pow_refl_channels
                sig_eq_channels = bs_ue_channels + np.sqrt(eta) * sig_refl_channels

                #####

                # # Check average power ratio
                # pow_avg_power_ratio = np.mean(np.abs(bs_ue_channels) ** 2) / np.mean(np.abs(pow_refl_channels) ** 2)
                # sig_avg_power_ratio = np.mean(np.abs(bs_ue_channels) ** 2) / np.mean(np.abs(sig_refl_channels) ** 2)
                #
                # print((np.abs(los_bs_ue_channels)**2).mean())
                # print((np.abs(los_ris_ue_channels)**2).mean())
                #
                # print((np.abs(bs_ue_channels) ** 2).mean())
                # print((np.abs(ris_ue_channels) ** 2).mean())
                #
                # print(pow_avg_power_ratio)
                # print(sig_avg_power_ratio)
                #
                # breakpoint()

                #####

                # Get channel estimates
                pow_hat_eq_channels = bs_rx_chest(eta, P_ue, n_pilots, sigma2_n, n_pilot_subblocks,
                                                  n_pilot_subblocks_probe, pow_eq_channels_probe, pow_eq_channels)

                sig_hat_eq_channels = bs_rx_chest(eta, P_ue, n_pilots, sigma2_n, n_pilot_subblocks,
                                                  n_pilot_subblocks_probe, sig_eq_channels_probe, sig_eq_channels)

                # Compute normalized mean squared error
                diff = pow_hat_eq_channels - pow_eq_channels
                pow_avg_nmse[cc, ss, :, nn] = (
                            np.linalg.norm(diff, axis=1) ** 2 / np.linalg.norm(pow_eq_channels, axis=1)).mean(axis=0)

                diff = sig_hat_eq_channels - sig_eq_channels
                sig_avg_nmse[cc, ss, :, nn] = (
                            np.linalg.norm(diff, axis=1) ** 2 / np.linalg.norm(sig_eq_channels, axis=1)).mean(axis=0)

                ##################################################
                # Communication Phase
                ##################################################
                pow_se, pow_num, pow_den1, pow_den2 = bs_comm(P_ue, sigma2_n, pow_eq_channels, pow_hat_eq_channels)
                sig_se, sig_num, sig_den1, sig_den2 = bs_comm(P_ue, sigma2_n, sig_eq_channels, sig_hat_eq_channels)

                # Store results
                pow_avg_se[cc, ss, :, nn] = pow_se.mean(axis=0)
                sig_avg_se[cc, ss, :, nn] = sig_se.mean(axis=0)

                pow_avg_num[cc, ss, :, nn] = pow_num.mean(axis=0)
                sig_avg_num[cc, ss, :, nn] = sig_num.mean(axis=0)

                pow_avg_den1[cc, ss, :, nn] = pow_den1.mean(axis=0)
                sig_avg_den1[cc, ss, :, nn] = sig_den1.mean(axis=0)

                pow_avg_den2[cc, ss, :, nn] = pow_den2.mean(axis=0)
                sig_avg_den2[cc, ss, :, nn] = sig_den2.mean(axis=0)

    # Compute relative probe time
    relative_probe_time = n_pilot_subblocks_probe_range / n_pilot_subblocks

    np.savez('data/pow-ris_K' + str(K) + '.npz',
             relative_probe_time=relative_probe_time,
             pow_avg_nmse=pow_avg_nmse,
             pow_avg_se=pow_avg_se,
             pow_avg_num=pow_avg_num,
             pow_avg_den1=pow_avg_den1,
             pow_avg_den2=pow_avg_den2
             )

    np.savez('data/sig-ris_K' + str(K) + '.npz',
             relative_probe_time=relative_probe_time,
             sig_avg_nmse=sig_avg_nmse,
             sig_avg_se=sig_avg_se,
             sig_avg_num=sig_avg_num,
             sig_avg_den1=sig_avg_den1,
             sig_avg_den2=sig_avg_den2
             )

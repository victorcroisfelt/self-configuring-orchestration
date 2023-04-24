import numpy as np
from scipy.constants import speed_of_light

from tqdm import tqdm

from src.channel import array_steering_vector, pathloss, scenario, drop_ues
from src.mmimo import bs_comm


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
    # UE Parameters
    ##################################################

    # Transmit power at the UE = 10 dBm
    P_ue = 10 ** ((0 - 30) / 10)

    ##################################################
    # System Parameters
    ##################################################

    # Coherence interval length
    tau_c = 128

    ##################################################
    # Geometry
    ##################################################

    # Physical parameters
    freq = 28 * 10 ** 9
    wavelength = speed_of_light / freq

    # NLoS variances
    sigma2_dr = 0.1 * 6.1848 * 1e-12
    sigma2_rr = 0.1 * 5.9603 * 1e-4

    # Noise power
    sigma2_n = 10 ** ((-94 - 30) / 10)

    # Generate scenario
    pos_bs, pos_bs_els, pos_ris, _, _, _, guard_distance_ris = scenario(wavelength, M, 32)

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

    # Krange
    K_range = np.array([1, 2, 4, 8, 16])
    n_yaxis = len(K_range)

    ##################################################
    # Simulation
    ##################################################

    # Prepare to save results
    avg_nmse = np.zeros((n_yaxis, n_setups, n_channels, n_noise))

    avg_se = np.zeros((2, n_yaxis, n_setups, n_channels, n_noise))

    avg_num = np.zeros((2, n_yaxis, n_setups, n_channels, n_noise))
    avg_den1 = np.zeros((2, n_yaxis, n_setups, n_channels, n_noise))
    avg_den2 = np.zeros((2, n_yaxis, n_setups, n_channels, n_noise))

    avg_sir = np.zeros((2, n_yaxis, n_setups, n_channels, n_noise))

    avg_nmse[:] = np.nan

    avg_se[:] = np.nan

    avg_num[:] = np.nan
    avg_den1[:] = np.nan
    avg_den2[:] = np.nan

    avg_sir[:] = np.nan

    # Go through all setups
    for ss in tqdm(range(n_setups)):

        # Go through all number of UEs
        for kk, K in enumerate(K_range):

            K = int(K)

            # Number of pilots
            n_pilots = K

            # Number of pilot subblocks
            n_pilot_subblocks = int(64 // K)

            # Calculate pre-log term
            pre_log_term = (tau_c - n_pilot_subblocks * n_pilots) / tau_c

            # Drop the UEs over the area of interest
            pos_ues = drop_ues(K, pos_ris, dmax=distance_max, guard_distance_ris=guard_distance_ris)

            # Compute LoS components of UE related channels
            bs_ue_steering = array_steering_vector(wavelength, pos_ues, pos_bs, pos_bs_els)
            bs_ue_pathloss = pathloss(3.76, pos_bs, pos_ues)

            los_bs_ue_channels = np.sqrt(bs_ue_pathloss)[:, None] * bs_ue_steering

            # Generate fading channels
            bs_ue_channels = np.sqrt(sigma2_dr / 2) * (np.random.randn(n_channels) + 1j * np.random.randn(n_channels))
            bs_ue_channels = los_bs_ue_channels[:, :, None] + bs_ue_channels[None, None, :]

            # Go through noise realizations
            for nn in range(n_noise):

                # Generate estimation noise
                noise = np.random.randn(K, M, n_channels) + 1j * np.random.randn(K, M, n_channels)
                noise *= np.sqrt(sigma2_n / 2 / n_pilot_subblocks / P_ue / n_pilots)

                # Get equivalent channel estimates
                hat_bs_ue_channels = bs_ue_channels + noise

                # Compute normalized mean squared error
                diff = hat_bs_ue_channels - bs_ue_channels
                avg_nmse[kk, ss, :, nn] = (np.linalg.norm(diff, axis=1) ** 2 / np.linalg.norm(bs_ue_channels, axis=1)).mean(axis=0)

                # Communication phase with MR
                se, num, den1, den2 = bs_comm(P_ue, sigma2_n, bs_ue_channels, hat_bs_ue_channels, method='MR')

                # Store results
                avg_se[0, kk, ss, :, nn] = pre_log_term * se.sum(axis=0)

                avg_num[0, kk, ss, :, nn] = num.mean(axis=0)
                avg_den1[0, kk, ss, :, nn] = den1.mean(axis=0)
                avg_den2[0, kk, ss, :, nn] = den2.mean(axis=0)

                avg_sir[0, kk, ss, :, nn] = (num/den1).mean(axis=0)

                # Communication phase with ZF
                se, num, den1, den2 = bs_comm(P_ue, sigma2_n, bs_ue_channels, hat_bs_ue_channels, method='ZF')

                # Store results
                avg_se[1, ss, :, nn] = pre_log_term * se.sum(axis=0)

                avg_num[1, kk, ss, :, nn] = num.mean(axis=0)
                avg_den1[1, kk, ss, :, nn] = den1.mean(axis=0)
                avg_den2[1, kk, ss, :, nn] = den2.mean(axis=0)

                avg_sir[1, kk, ss, :, nn] = (num / den1).mean(axis=0)

    np.savez('data/figure7_mmimo_K' + str(K) + '.npz',
             K_range=K_range,
             avg_nmse=avg_nmse,
             avg_se=avg_se,
             avg_num=avg_num,
             avg_den1=avg_den1,
             avg_den2=avg_den2,
             avg_sir=avg_sir
             )

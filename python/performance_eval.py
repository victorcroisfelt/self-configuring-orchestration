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

    # Transmit power at the BS = 20 dBm
    P_bs = 10 ** ((20 - 30) / 10)

    ##################################################
    # UE Parameters
    ##################################################

    # Number of UEs
    K = 2

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
    # Geometry
    ##################################################

    # Physical parameters
    freq = 3 * 10 ** 9
    wavelength = speed_of_light / freq

    # NLoS variances
    sigma2_dr = 1e-10

    # Noise power
    sigma2_n = 10 ** ((-94 - 30) / 10)

    # Parameterized variable to define area of interest
    scenario_size = 100

    # Guard distance
    guard_distance = 30

    # Generate scenario
    pos_bs, pos_bs_els, pos_ris, _, _, _ = scenario(wavelength, scenario_size, M, 32)

    ##################################################
    # Simulation Parameters
    ##################################################

    # Define number of setups
    n_setups = 10

    # Define number of channel realizations
    n_channels = 100

    # Define number of noise realizations
    n_noise = 100

    ##################################################
    # Simulation
    ##################################################

    # Prepare to save results
    avg_nmse = np.zeros((n_setups, n_channels, n_noise))

    avg_se = np.zeros((n_setups, n_channels, n_noise))

    avg_num = np.zeros((n_setups, n_channels, n_noise))
    avg_den1 = np.zeros((n_setups, n_channels, n_noise))
    avg_den2 = np.zeros((n_setups, n_channels, n_noise))

    # Go through all setups
    for ss in tqdm(range(n_setups)):

        # Drop the UEs over the area of interest
        pos_ues = drop_ues(K, pos_ris, scenario_size, guard_distance)

        # Compute LoS components of UE related channels
        bs_ue_steering = array_steering_vector(wavelength, pos_ues, pos_bs, pos_bs_els)
        bs_ue_pathloss = pathloss(3.76, pos_bs, pos_ues)

        los_bs_ue_channels = np.sqrt(bs_ue_pathloss)[:, None] * bs_ue_steering

        # Generate fading channels
        bs_ue_channels = los_bs_ue_channels[:, :, None] + \
                         np.sqrt(sigma2_dr / 2) * (np.random.randn(n_channels) + 1j * np.random.randn(n_channels))

        # Go through noise realizations
        for nn in range(n_noise):

            ##################################################
            # Chest Phase
            ##################################################

            # Generate noise after decorrelation
            noise = (np.random.randn(n_pilot_subblocks, K, M, n_channels) + 1j * np.random.randn(n_pilot_subblocks, K, M, n_channels))
            noise *= np.sqrt((n_pilots * sigma2_n) / 2)

            # Compute RX signal over subblocks
            rx_signal = np.sqrt(P_ue) * n_pilots * bs_ue_channels + noise

            # Take the average over subblocks
            average_subblocks = np.mean(rx_signal, axis=0)

            # Get equivalent channel estimates
            hat_bs_ue_channels = average_subblocks / np.sqrt(P_ue) / n_pilots

            # Compute normalized mean squared error
            diff = hat_bs_ue_channels - bs_ue_channels
            avg_nmse[ss, :, nn] = (np.linalg.norm(diff, axis=1) ** 2 / np.linalg.norm(bs_ue_channels, axis=1)).mean(
                axis=0)

            # Communication phase
            se, num, den1, den2 = bs_comm(P_ue, sigma2_n, bs_ue_channels, hat_bs_ue_channels)

            # Store results
            avg_se[ss, :, nn] = se.mean(axis=0)
            avg_num[ss, :, nn] = num.mean(axis=0)
            avg_den1[ss, :, nn] = den1.mean(axis=0)
            avg_den2[ss, :, nn] = den2.mean(axis=0)

    np.savez('data/mmimo_K' + str(K) + '.npz',
             avg_nmse=avg_nmse,
             avg_se=avg_se,
             avg_num=avg_num,
             avg_den1=avg_den1,
             avg_den2=avg_den2
             )

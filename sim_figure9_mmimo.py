import numpy as np
from scipy.constants import speed_of_light

from tqdm import tqdm

from src.channel import scenario, drop_ues, generate_channel_realizations
from src.mmimo import bs_comm

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

    # Transmit power at the UE = 0 dBm
    P_ue = 10 ** ((0 - 30) / 10)

    ##################################################
    # System Parameters
    ##################################################

    # Coherence interval length
    coherence_interval_length = 128

    ##################################################
    # Geometry
    ##################################################

    # Physical parameters
    freq = 28 * 10 ** 9
    wavelength = speed_of_light / freq

    # NLoS variances
    sigma2_dr = 0.1 * 9.07 * 1e-9
    sigma2_rr = 0.1 * 1.12 * 1e-6

    # Noise power
    sigma2_n = 10 ** ((-94 - 30) / 10)

    # Generate scenario
    pos_bs, pos_bs_els, pos_ris, pos_ris_els, bs_ris_channels, ris_bs_steering, guard_distance_ris = scenario(wavelength, M, N=32)

    ##################################################
    # Simulation Parameters
    ##################################################

    # Define number of setups
    n_setups = 100

    # Define number of channel realizations
    n_channels = 100

    # Define number of noise realizations
    n_noise = 100

    # Define the number of UEs
    K = 4

    # Number of pilots
    n_pilots = int(K)

    # Number of pilot subblocks
    n_pilot_subblocks = int(64 // K)

    # Calculate pre-log term
    pre_log_term = (coherence_interval_length - n_pilot_subblocks * n_pilots) / coherence_interval_length

    ##################################################
    # Simulation
    ##################################################

    # Prepare to save results
    avg_nmse = np.zeros((n_setups, n_noise, n_channels, K))
    avg_se = np.zeros((2, n_setups, n_noise, n_channels, K))
    avg_se_pos = np.zeros((2, n_setups, n_noise, n_channels, K))

    avg_nmse[:] = np.nan
    avg_se[:] = np.nan
    avg_se_pos[:] = np.nan

    # Go through all setups
    for ss in tqdm(range(n_setups)):

        # Drop the UEs over the area of interest
        pos_ues = drop_ues(K, pos_ris, dmax=1000, guard_distance_ris=900)

        # Generate UE channels
        bs_ue_channels, _ = generate_channel_realizations(wavelength, pos_bs, pos_bs_els, pos_ris, pos_ris_els, pos_ues, sigma2_dr, sigma2_rr, n_channels)

        # Go through noise realizations
        for nn in range(n_noise):

            # Generate estimation noise
            noise = np.random.randn(M, n_channels, K) + 1j * np.random.randn(M, n_channels, K)
            noise *= np.sqrt(sigma2_n / 2 / n_pilot_subblocks / P_ue / n_pilots)

            # Get equivalent channel estimates
            hat_bs_ue_channels = bs_ue_channels + noise

            # Compute normalized mean squared error
            diff = hat_bs_ue_channels - bs_ue_channels
            avg_nmse[ss, nn, :, :] = (np.linalg.norm(diff, axis=0) ** 2 / np.linalg.norm(bs_ue_channels, axis=0))

            # Communication phase with MR
            se, _, _, _ = bs_comm(P_ue, sigma2_n, bs_ue_channels, hat_bs_ue_channels, method='MR')

            # Store results
            avg_se[0, ss, :, nn] = se
            avg_se_pos[0, ss, :, nn] = pre_log_term * se

            # Communication phase with ZF
            se, _, _, _ = bs_comm(P_ue, sigma2_n, bs_ue_channels, hat_bs_ue_channels, method='ZF')

            # Store results
            avg_se[1, ss, :, nn] = se
            avg_se_pos[1, ss, :, nn] = pre_log_term * se

    np.savez('data/figure9_mmimo_K' + str(K) + '.npz',
             avg_nmse=avg_nmse,
             avg_se=avg_se,
             avg_se_pos=avg_se_pos
             )

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
    # Scenario Parameters
    ##################################################

    # Physical parameters
    freq = 28 * 10 ** 9
    wavelength = speed_of_light / freq

    # NLoS variances
    sigma2_dr = 0.1 * 9.08 * 1e-7
    sigma2_rr = 0.1 * 1.11 * 1e-6

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

    # Number of pilots
    n_pilots = K

    ##################################################
    # Simulation
    ##################################################

    # Prepare to save results
    avg_se = np.zeros((2, n_setups, n_noise, n_channels, K))
    avg_se[:] = np.nan

    # Go through all setups
    for ss in tqdm(range(n_setups), desc='setups'):

        # Drop the UEs over the area of interest
        pos_ues = drop_ues(K, pos_ris, dmax=1000, guard_distance_ris=900)

        # Generate UE channels
        bs_ue_channels, ris_ue_channels = generate_channel_realizations(
            wavelength, pos_bs, pos_bs_els, pos_ris, pos_ris_els, pos_ues, sigma2_dr, sigma2_rr, n_channels)
        
        # Genie reflection configuration
        gen_reflection_configs, gen_weights = gen_ris_probe(ris_ue_channels)
        gen_reflection_configs *= ris_bs_steering[:, None]

        # Go through noise realizations
        for nn in range(n_noise):

            # Compute reflected channels during communication
            gen_refl_channels = gen_reflection_configs[:, :, None] * ris_ue_channels
            gen_refl_channels = bs_ris_channels[:, :, None, None] * gen_refl_channels[None, :, :, :]
            gen_refl_channels = gen_refl_channels.sum(axis=1)

            # Compute equivalent channels during communication
            gen_eq_channels = bs_ue_channels + gen_refl_channels

            ##################################################
            # Communication Phase
            ##################################################
            gen_se, _, _, _ = bs_comm(P_ue, sigma2_n_bs, gen_eq_channels, gen_eq_channels)

            # Store results
            avg_se[0, ss, nn, :, :] = gen_se

    np.savez('data/figure9_controlled-ris_K' + str(K) + '_N' + str(N) + '.npz',
            avg_se=avg_se
            )

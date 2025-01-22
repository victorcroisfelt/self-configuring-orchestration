import numpy as np

from tqdm import tqdm

from scipy.constants import speed_of_light
from scipy.stats.distributions import chi2

from src.channel import generate_channel_realizations, scenario, drop_ues
from src.ris import pow_ris_config_codebook, ris_rx_chest_with_choice

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

    # Transmit power; 0 dBm
    P_ue = 10 ** ((0 - 30) / 10)

    ##################################################
    # System Parameters
    ##################################################

    # Number of pilots
    n_pilots = K

    # Number of pilot subblocks
    n_pilot_subblocks = 16

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
    sigma2_n_hris = 10 ** ((-91 - 30) / 10)

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
    pd_proba_detection = np.zeros((n_probe, n_etas, n_probas))
    dsp_proba_detection = np.zeros((n_probe, n_etas, n_probas))

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
                pow_ris_rx_chest = ris_rx_chest_with_choice(eta, P_ue, n_pilots, sigma2_n_hris, n_pilot_subblocks_probe, ris_ue_channels, mask, pow_probe_configs)
                sig_ris_rx_chest = ris_rx_chest_with_choice(eta, P_ue, n_pilots, sigma2_n_hris, n_pilot_subblocks_probe, ris_ue_channels, mask)

                #
                # PD-based HRIS probe
                #

                # Compute measured power
                pow_test = np.abs(pow_ris_rx_chest)**2
                pow_test = pow_test.max(axis=0)

                #
                # DSP-based HRIS probe
                #

                # Compute test
                sig_test = np.linalg.norm(sig_ris_rx_chest, axis=0)**2
                sig_test *= 2 * n_pilot_subblocks_probe / n_pilots / sigma2_n_hris

                # Go through all false alarm probabilities
                for pp in range(n_probas):

                    # Get current value of eta
                    proba_false_alarm = proba_false_alarm_range[pp]

                    #
                    # PD-based HRIS probe
                    #

                    # Compute threshold
                    pow_threshold = - (2 * N * sigma2_n_hris) * np.log(proba_false_alarm)

                    # Perform test
                    pow_detected_ues = pow_test > pow_threshold

                    # Estimate probability of detection
                    pd_proba_detection_curr = (ue_choice == pow_detected_ues)

                    # Store results
                    pd_proba_detection[cc, ee, pp] += np.mean(pd_proba_detection_curr)

                    #
                    # DSP-based HRIS probe
                    #

                    # Get threshold value
                    sig_threshold = chi2.ppf((1 - proba_false_alarm), df=2*N)

                    # Perform test
                    sig_detected_ues = sig_test > sig_threshold
           
                    # Estimate probability of detection
                    dsp_proba_detection_curr = (ue_choice == sig_detected_ues)
       
                    # Store results
                    dsp_proba_detection[cc, ee, pp] += np.mean(dsp_proba_detection_curr)

    np.savez('data/figure7.npz',
        n_pilot_subblocks=n_pilot_subblocks,
        n_pilot_subblocks_probe_range=n_pilot_subblocks_probe_range,
        eta_range=eta_range,
        proba_false_alarm_range=proba_false_alarm_range,
        pd_proba_detection=pd_proba_detection,
        dsp_proba_detection=dsp_proba_detection
    )
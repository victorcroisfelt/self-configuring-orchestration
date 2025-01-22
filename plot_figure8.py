import numpy as np
from scipy.constants import speed_of_light

import matplotlib.pyplot as plt

from src.channel import scenario, generate_los_channel_realizations, array_steering_vector
from src.ris import pow_ris_config_codebook, ris_rx_chest, gen_ris_probe, pow_ris_probe, sig_ris_probe

from matplotlib import rc
rc('font', **{'family': 'sans serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

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
    eta = 0.999

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

    # Number of pilot subblocks
    n_pilot_subblocks = 16

    # Number of probe pilot subblocks
    n_probe_pilot_subblocks = 8

    ##################################################
    # Scenario Parameters
    ##################################################

    # Physical parameters
    freq = 28 * 10 ** 9
    wavelength = speed_of_light / freq

    # Noise power
    sigma2_n_bs = 10 ** ((-94 - 30) / 10)
    sigma2_n_ris = 10 ** ((-91 - 30) / 10)

    # Generate scenario
    pos_bs, pos_bs_els, pos_ris, pos_ris_els, bs_ris_channels, ris_bs_steering, guard_distance_ris = scenario(wavelength, M, N)

    ##################################################
    # Simulation Parameters
    ##################################################

    # Define probability of false alarm
    proba_false_alarm = 0.01

    ##################################################
    # Simulation
    ##################################################

    # Drop the UEs over the area of interest
    angles_ues_deg = np.array([30, 75])
    angles_ues = np.deg2rad(angles_ues_deg)
    distance_ues = np.array([10, 20])

    # Get UE positions
    pos_ues = distance_ues[:, None] * np.array([np.cos(angles_ues), np.sin(angles_ues)]).T

    # Generate UE channels
    los_bs_ue_channels, los_ris_ue_channels = generate_los_channel_realizations(wavelength, pos_bs, pos_bs_els, pos_ris, pos_ris_els, pos_ues)

    # Generate PD-enabled HRIS configuration codebook
    pow_probe_configs = pow_ris_config_codebook(wavelength, 8, pos_ris, pos_ris_els)

    ##################################################
    # Chest Phase
    ##################################################
    
    # Compute received signals during probe
    pow_ris_rx_chest = ris_rx_chest(eta, P_ue, n_pilots, sigma2_n_ris, 8, los_ris_ue_channels[:, None, :], pow_probe_configs)
    sig_ris_rx_chest = ris_rx_chest(eta, P_ue, n_pilots, sigma2_n_ris, 8, los_ris_ue_channels[:, None, :])

    # HRIS probe mode
    gen_reflection_configs, gen_weights = gen_ris_probe(los_ris_ue_channels[:, None, :])
    pow_reflection_configs, pow_weights, pow_hat_aoa = pow_ris_probe(N, sigma2_n_ris, proba_false_alarm, pow_ris_rx_chest, pow_probe_configs)
    sig_reflection_configs, sig_weights, sig_hat_aoa = sig_ris_probe(n_probe_pilot_subblocks, n_pilots, sigma2_n_ris, proba_false_alarm, sig_ris_rx_chest)

    # Average configurations
    gen_reflection_configs = np.squeeze(gen_reflection_configs)
    pow_reflection_configs = np.squeeze(pow_reflection_configs)
    sig_reflection_configs = np.squeeze(sig_reflection_configs)

    ####
    # Directing
    ####

    # Get angular position of the UEs
    ang_ues = np.angle(pos_ues[:, 0] + 1j * pos_ues[:, 1])

    # Create AoA positions
    angles = np.linspace(0, np.pi/2, 101)

    # Get respective positions
    pos = np.zeros((angles.size, 2))

    pos[:, 0] = np.cos(angles)
    pos[:, 1] = np.sin(angles)

    # Compute steering vector with respect to positions
    ris_angle_steering = array_steering_vector(wavelength, pos, pos_ris, pos_ris_els)
    ris_angle_steering = ris_angle_steering.T
    
    gen_inner = gen_reflection_configs[:, np.newaxis] * ris_angle_steering
    gen_inner = gen_inner.sum(axis=0)
    
    pow_inner = pow_reflection_configs[:, np.newaxis] * ris_angle_steering
    pow_inner = pow_inner.sum(axis=0)

    sig_inner = sig_reflection_configs[:, np.newaxis] * ris_angle_steering
    sig_inner = sig_inner.sum(axis=0)

    # Plot for PD-enabled HRIS
    fig, axes = plt.subplots()
    markers = ['1', '2', '3', '4']

    pow_abs = (np.abs(pow_inner)**2)
    angles_deg = np.rad2deg(angles)

    # Save the array to a text file
    pow_res = np.vstack((angles_deg, pow_abs)).T
    np.savetxt('txts/fig8_pow_res.txt', pow_res, fmt='%.4f', comments='')

    axes.plot(angles_deg, pow_abs, linewidth=1.5, color='black')

    axes.set_xticks([0, 30, 60, 90])
    axes.grid(color='#E9E9E9', linestyle=':', linewidth=1.0, alpha=0.5)

    axes.set_xlabel('UE angular position')
    axes.set_ylabel('channel gain')

    plt.tight_layout()

    # Plot for DSP-enabled HRIS
    fig, axes = plt.subplots()

    sig_abs = (np.abs(sig_inner)**2)

    sig_res = np.vstack((angles_deg, sig_abs)).T
    np.savetxt('txts/fig8_sig_res.txt', sig_res, fmt='%.4f', comments='')

    axes.plot(np.rad2deg(angles), sig_abs, linewidth=1.5, color='black')
    axes.set_xticks([0, 30, 60, 90])
    axes.grid(color='#E9E9E9', linestyle=':', linewidth=1.0, alpha=0.5)

    axes.set_xlabel('UE angular position')
    axes.set_ylabel('channel gain ')

    plt.tight_layout()

    print('--------------------')
    print('NMSE')
    print('--------------------')

    nmse_pow = np.linalg.norm(pow_reflection_configs - gen_reflection_configs)**2 / np.linalg.norm(gen_reflection_configs)**2 
    nmse_sig = np.linalg.norm(sig_reflection_configs - gen_reflection_configs)**2 / np.linalg.norm(gen_reflection_configs)**2 

    print('PD-enabled NMSE =', nmse_pow)
    print('DSP-enabled NMSE =', nmse_sig)

    plt.show()
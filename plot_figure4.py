import numpy as np
from scipy.constants import speed_of_light

import matplotlib.pyplot as plt

from src.channel import scenario, generate_los_channel_realizations
from src.ris import pow_ris_config_codebook, ris_rx_chest, gen_ris_probe, pow_ris_probe, sig_ris_probe

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
    eta = 0.9999

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
    n_probe_pilot_subblocks = 16

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

    # Drop a single UE
    angles_ues_deg = np.array([45, ])
    angles_ues = np.deg2rad(angles_ues_deg)
    distance_ues = np.array([100, ])

    # Get UE position
    pos_ues = distance_ues[:, None] * np.array([np.cos(angles_ues), np.sin(angles_ues)]).T

    # Generate UE channel realizations
    los_bs_ue_channels, los_ris_ue_channels = generate_los_channel_realizations(wavelength, pos_bs, pos_bs_els, pos_ris, pos_ris_els, pos_ues)

    # Generate PD-enabled HRIS configuration codebook
    pow_probe_configs = pow_ris_config_codebook(wavelength, n_probe_pilot_subblocks, pos_ris, pos_ris_els)

    ##################################################
    # Chest Phase
    ##################################################

    # Compute received signals during probe
    pow_ris_rx_chest = ris_rx_chest(eta, P_ue, n_pilots, sigma2_n_ris, n_probe_pilot_subblocks, los_ris_ue_channels[:, None, :], pow_probe_configs)
    sig_ris_rx_chest = ris_rx_chest(eta, P_ue, n_pilots, sigma2_n_ris, n_probe_pilot_subblocks, los_ris_ue_channels[:, None, :])

    # HRIS probe
    gen_reflection_configs, gen_weights = gen_ris_probe(los_ris_ue_channels[:, :, None])
    pow_reflection_configs, pow_weights, pow_hat_aoa = pow_ris_probe(N, sigma2_n_ris, proba_false_alarm, pow_ris_rx_chest, pow_probe_configs)
    sig_reflection_configs, sig_weights, sig_hat_aoa = sig_ris_probe(n_pilots, sigma2_n_ris, proba_false_alarm, sig_ris_rx_chest)

    # Complete reflection configurations by inserting BS-RIS knowledge
    gen_reflection_configs *= ris_bs_steering[:, None]
    pow_reflection_configs *= ris_bs_steering[:, None]
    sig_reflection_configs *= ris_bs_steering[:, None]

    # Compute reflected channels during probe
    pow_refl_channels_probe = pow_probe_configs[:, :, None] * los_ris_ue_channels[None, :, :]
    pow_refl_channels_probe = bs_ris_channels[None, :, :, None] * pow_refl_channels_probe[:, None, :, :]
    
    pow_refl_channels_probe = pow_refl_channels_probe.sum(axis=2)
    
    ones = np.ones_like(pow_probe_configs)
    sig_refl_channels_probe = ones[:, :, None] * los_ris_ue_channels[None, :, :]
    sig_refl_channels_probe = bs_ris_channels[None, :, :, None] * sig_refl_channels_probe[:, None, :, :]
    
    sig_refl_channels_probe = sig_refl_channels_probe.sum(axis=2)
    
    # Compute equivalent channels during probe
    pow_eq_channels_probe = los_bs_ue_channels[None, :, :] + np.sqrt(eta) * pow_refl_channels_probe
    sig_eq_channels_probe = los_bs_ue_channels[None, :, :] + np.sqrt(eta) * sig_refl_channels_probe
  
    # Compute reflected channels during communication
    gen_refl_channels = gen_reflection_configs[:, :] * los_ris_ue_channels 
    gen_refl_channels = bs_ris_channels[:, :, None] * gen_refl_channels[None, :, :]
    gen_refl_channels = gen_refl_channels.sum(axis=1)

    pow_refl_channels = pow_reflection_configs[:, :] * los_ris_ue_channels
    pow_refl_channels = bs_ris_channels[:, :, None] * pow_refl_channels[None, :, :]
    pow_refl_channels = pow_refl_channels.sum(axis=1)
    
    sig_refl_channels = sig_reflection_configs[:, :] * los_ris_ue_channels
    sig_refl_channels = bs_ris_channels[:, :, None] * sig_refl_channels[None, :, :]
    sig_refl_channels = sig_refl_channels.sum(axis=1)

    # Compute equivalent channels during communication
    gen_eq_channels = los_bs_ue_channels + np.sqrt(eta) * gen_refl_channels
    pow_eq_channels = los_bs_ue_channels + np.sqrt(eta) * pow_refl_channels
    sig_eq_channels = los_bs_ue_channels + np.sqrt(eta) * sig_refl_channels

    # Plot figure
    fig, axes = plt.subplots()
    
    # Squeeze out user's dimension
    probe = np.squeeze(pow_eq_channels_probe)
    refle = np.squeeze(pow_eq_channels)

    Y_a = 10 * np.log10(np.linalg.norm(probe, axis=1)**2)
    Y_b = 10 * np.log10(np.linalg.norm(refle)**2)

    Y2_a = 10 * np.log10(np.linalg.norm(los_bs_ue_channels[0])**2)

    X_ = np.arange(1, 32+1)
    Y_ = np.concatenate((Y_a, Y_b * np.ones(16)))
 
    # Save the array to a text file
    figure4 = np.vstack((X_, Y_)).T
    np.savetxt('txts/fig4_hris-assisted.txt', figure4, fmt='%.4f', comments='')

    Y2_ = Y2_a * np.ones(32)
    figure4 = np.vstack((X_, Y2_)).T
    np.savetxt('txts/fig4_standalone.txt', figure4, fmt='%.4f', comments='')

    axes.plot(X_[0:16], Y_[0:16], linewidth=1.5, color='red', label='Probing')
    axes.plot(X_[16:], Y_[16:], linewidth=1.5, color='blue', label='Reflecting')

    axes.plot(X_, Y2_, linewidth=1.5, color='black', label='Standalone')

    axes.set_ylabel('Equivalent BS-UE channel gain [dB]')
    axes.set_xlabel('Subblocks, $s$')

    plt.legend()

    plt.tight_layout()

    plt.show()
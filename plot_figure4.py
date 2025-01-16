import numpy as np
from scipy.constants import speed_of_light

import matplotlib.pyplot as plt

from src.channel import scenario, drop_ues, generate_channel_realizations2, array_steering_vector
from src.ris import pow_ris_config_codebook, ris_rx_chest, gen_ris_probe, pow_ris_probe, sig_ris_probe

#from matplotlib import rc

#import tikzplotlib

# rc('font', **{'family': 'sans serif', 'serif': ['Computer Modern']})
# rc('text', usetex=True)

# Press the green button in the gutter to run the script.
if __name__ == '__main__':

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
    # sigma2_n = 10 ** ((-94 - 30) / 10)

    # Maximum distances of the UEs
    maximum_distance = 100

    # NLoS variances
    sigma2_dr = 0.1 * 9.07 * 1e-9
    sigma2_rr = 0.1 * 1.12 * 1e-6

    # Noise power
    sigma2_n_bs = 10 ** ((-94 - 30) / 10)
    sigma2_n_ris = 10 ** ((-91 - 30) / 10)

    # Generate scenario
    pos_bs, pos_bs_els, pos_ris, pos_ris_els, bs_ris_channels, ris_bs_steering, guard_distance_ris = scenario(wavelength, 64, N)

    ##################################################
    # Simulation Parameters
    ##################################################

    # Define probability of false alarm
    proba_false_alarm = 0.01

    ##################################################
    # Simulation
    ##################################################

    # Drop the UEs over the area of interest
    angles_ues_deg = np.array([30, 80])
    angles_ues = np.deg2rad(angles_ues_deg)
    distance_ues = np.array([100, 100])

    # Get UE positions
    pos_ues = distance_ues[:, None] * np.array([np.cos(angles_ues), np.sin(angles_ues)]).T

    # Generate UE channels
    bs_ue_channels, los_bs_ue_channels, ris_ue_channels, los_ris_ue_channels = generate_channel_realizations2(
        wavelength, pos_bs, pos_bs_els, pos_ris, pos_ris_els, pos_ues, sigma2_dr, sigma2_rr, 1)

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
    pow_reflection_configs, pow_weights, pow_hat_aoa, pow_pd = pow_ris_probe(N, sigma2_n_ris, proba_false_alarm, pow_ris_rx_chest, pow_probe_configs)
    sig_reflection_configs, sig_weights, sig_hat_aoa, sig_pd = sig_ris_probe(n_pilots, sigma2_n_ris, proba_false_alarm, sig_ris_rx_chest)

    # Complete reflection configurations by inserting BS-RIS knowledge
    gen_reflection_configs *= ris_bs_steering[:, None]
    pow_reflection_configs *= ris_bs_steering[:, None]
    sig_reflection_configs *= ris_bs_steering[:, None]

    # Compute reflected channels during probe
    pow_refl_channels_probe = pow_probe_configs[:, :, None] * los_ris_ue_channels[None, :, :]
    pow_refl_channels_probe = bs_ris_channels[None, :, :, None] * pow_refl_channels_probe[:, None, :, :]
    
    pow_refl_channels_probe = pow_refl_channels_probe.sum(axis=2)
    #pow_refl_channels_probe = pow_refl_channels_probe.sum(axis=0)
    
    ones = np.ones_like(pow_probe_configs)
    sig_refl_channels_probe = ones[:, :, None] * los_ris_ue_channels[None, :, :]
    sig_refl_channels_probe = bs_ris_channels[None, :, :, None] * sig_refl_channels_probe[:, None, :, :]
    
    sig_refl_channels_probe = sig_refl_channels_probe.sum(axis=2)
    #sig_refl_channels_probe = sig_refl_channels_probe.sum(axis=0)
    
    # Compute equivalent channels during probe
    pow_eq_channels_probe = los_bs_ue_channels[None, :, :] + np.sqrt(eta) * pow_refl_channels_probe
    sig_eq_channels_probe = los_bs_ue_channels[None, :, :] + np.sqrt(eta) * sig_refl_channels_probe

    # Compute reflected channels during communication
    gen_refl_channels = gen_reflection_configs[:, :] * los_ris_ue_channels 
    gen_refl_channels = bs_ris_channels[:, :, None] * gen_refl_channels[None, :, :]
    gen_refl_channels = gen_refl_channels.sum(axis=1)

    pow_refl_channels = pow_reflection_configs[:, :] * los_ris_ue_channels
    #breakpoint()
    pow_refl_channels = bs_ris_channels[:, :, None] * pow_refl_channels[None, :, :]
    #breakpoint()
    pow_refl_channels = pow_refl_channels.sum(axis=1)
    
    sig_refl_channels = sig_reflection_configs[:, :] * los_ris_ue_channels
    sig_refl_channels = bs_ris_channels[:, :, None] * sig_refl_channels[None, :, :]
    sig_refl_channels = sig_refl_channels.sum(axis=1)

    # Compute equivalent channels during communication
    gen_eq_channels = los_bs_ue_channels + np.sqrt(eta) * gen_refl_channels
    pow_eq_channels = los_bs_ue_channels + np.sqrt(eta) * pow_refl_channels
    sig_eq_channels = los_bs_ue_channels + np.sqrt(eta) * sig_refl_channels

    #breakpoint()

    ####
    # Directing
    ####

    # # Get angular position of the UEs
    # ang_ues = np.angle(pos_ues[:, 0] + 1j * pos_ues[:, 1])

    # # Create AoA positions
    # angles = np.linspace(0, np.pi/2, 101)

    # # Get respective positions
    # pos = np.zeros((angles.size, 2))

    # pos[:, 0] = np.cos(angles)
    # pos[:, 1] = np.sin(angles)

    # # Compute steering vector with respect to positions
    # ris_angle_steering = array_steering_vector(wavelength, pos, pos_ris, pos_ris_els)
    # ris_angle_steering = ris_angle_steering.T
    
    # gen_inner = gen_reflection_configs[:, np.newaxis] * ris_angle_steering
    # gen_inner = gen_inner.sum(axis=0)
    
    # pow_inner = pow_reflection_configs[:, np.newaxis] * ris_angle_steering
    # pow_inner = pow_inner.sum(axis=0)

    # sig_inner = sig_reflection_configs[:, np.newaxis] * ris_angle_steering
    # sig_inner = sig_inner.sum(axis=0)
    #breakpoint()
    # Plot figure
    fig, axes = plt.subplots()
    # markers = ['1', '2', '3', '4']
    
    probe = pow_eq_channels_probe[:, 0]
    refle = pow_eq_channels[:, 0]

    Y_a = 10 * np.log10(np.linalg.norm(probe, axis=1)**2)
    Y_b = 10 * np.log10(np.linalg.norm(refle)**2)
    Y2_a = 10 * np.log10(np.linalg.norm(los_bs_ue_channels[0])**2)

    breakpoint()
    X_ = np.arange(1, 32+1)
    Y_ = np.concatenate((Y_a, Y_b * np.ones(16)))
    breakpoint()    
    # Save the array to a text file
    figure4 = np.vstack((X_, Y_)).T
    np.savetxt('fig4_a.txt', figure4, fmt='%.4f', comments='')

    Y2_ = Y2_a * np.ones(32)
    figure4 = np.vstack((X_, Y2_)).T
    np.savetxt('fig4_b.txt', figure4, fmt='%.4f', comments='')


    axes.plot(X_, Y_, linewidth=1.5, color='black')
    axes.plot(X_, Y2_, linewidth=1.5, color='black')
    plt.show()

    # axes.set_xticks([0, 30, 60, 90, 120, 150, 180])

    #axes.plot(np.rad2deg(ang_ues[0]), np.min(np.abs(pow_inner)) ** 2 * np.ones_like(ang_ues[0]), linewidth=0.0, marker=markers[0], markersize=7, label=r'UE 1, $d_{1}=10$ m, $\theta_1=30^{\circ}$')
    #axes.plot(np.rad2deg(ang_ues[1]), np.min(np.abs(pow_inner)) ** 2 * np.ones_like(ang_ues[1]), linewidth=0.0, marker=markers[1], markersize=7, label=r'UE 2, $d_{2}=20$ m, $\theta_1=80^{\circ}$')
    #axes.plot(np.rad2deg(ang_ues[2]), np.min(np.abs(pow_inner)) ** 2 * np.ones_like(ang_ues[2]), linewidth=0.0, marker=markers[2], markersize=7, label=r'UE 3, $d_{3}=20$ m, $\theta_1=110^{\circ}$')
    #axes.plot(np.rad2deg(ang_ues[3]), np.min(np.abs(pow_inner)) ** 2 * np.ones_like(ang_ues[3]), linewidth=0.0, marker=markers[3], markersize=7, label=r'UE 4, $d_{4}=10$ m, $\theta_1=150^{\circ}$')

    # axes.grid(color='#E9E9E9', linestyle=':', linewidth=1.0, alpha=0.5)

    # axes.legend(framealpha=0.5)

    # axes.set_xlabel('UE angular position')
    # axes.set_ylabel('channel gain')

    # plt.tight_layout()

    # #tikzplotlib.save("tikz/figure6b.tex")

    # fig, axes = plt.subplots()

    # #sig_abs = 10 * np.log10(np.abs(sig_inner)**2)
    # sig_abs = (np.abs(sig_inner)**2)

    # sig_res = np.vstack((angles_deg, sig_abs)).T
    # np.savetxt('sig_res.txt', sig_res, fmt='%.4f', comments='')

    # axes.plot(np.rad2deg(angles), sig_abs, linewidth=1.5, color='black')

    # #axes.plot(np.rad2deg(ang_ues[0]), np.min(np.abs(sig_inner)) ** 2 * np.ones_like(ang_ues[0]), marker=markers[0], markersize=7)
    # #axes.plot(np.rad2deg(ang_ues[1]), np.min(np.abs(sig_inner)) ** 2 * np.ones_like(ang_ues[1]), marker=markers[1], markersize=7)
    # #axes.plot(np.rad2deg(ang_ues[2]), np.min(np.abs(sig_inner)) ** 2 * np.ones_like(ang_ues[2]), marker=markers[2], markersize=7)
    # #axes.plot(np.rad2deg(ang_ues[3]), np.min(np.abs(sig_inner)) ** 2 * np.ones_like(ang_ues[3]), marker=markers[3], markersize=7)

    # axes.set_xticks([0, 30, 60, 90, 120, 150, 180])

    # axes.grid(color='#E9E9E9', linestyle=':', linewidth=1.0, alpha=0.5)

    # axes.set_xlabel('UE angular position')
    # axes.set_ylabel('channel gain ')

    # plt.tight_layout()

    # # tikzplotlib.save("tikz/figure6c.tex")

    # #plt.show()

    # nmse_pow = np.linalg.norm(pow_reflection_configs - gen_reflection_configs)**2 / np.linalg.norm(gen_reflection_configs)**2 
    # nmse_sig = np.linalg.norm(sig_reflection_configs - gen_reflection_configs)**2 / np.linalg.norm(gen_reflection_configs)**2 

    # print(nmse_pow)
    # print(nmse_sig)
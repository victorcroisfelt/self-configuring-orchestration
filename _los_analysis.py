import numpy as np
from scipy.constants import speed_of_light

from tqdm import tqdm

from src.channel import scenario, drop_ues, generate_los_channel_realizations

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
    # Geometry
    ##################################################

    # Physical parameters
    freq = 6 * 10 ** 9
    wavelength = speed_of_light / freq

    # Generate scenario
    pos_bs, pos_bs_els, pos_ris, pos_ris_els, bs_ris_channels, ris_bs_steering, guard_distance_ris = scenario(wavelength, M, N)

    ##################################################
    # Simulation Parameters
    ##################################################

    # Define number of setups
    n_setups = 10000

    ##################################################
    # Simulation
    ##################################################

    # Prepare to save results
    los_bs_ue_channels = np.zeros((n_setups, M), dtype=np.complex128)
    los_ris_ue_channels = np.zeros((n_setups, N), dtype=np.complex128)

    # Go through all setups
    for ss in range(n_setups):

        # Drop the UEs over the area of interest
        pos_ues = drop_ues(1, pos_ris, dmax=1000, guard_distance_ris=900)

        # Generate LoS channels
        bs_ue, ris_ue = generate_los_channel_realizations(wavelength, pos_bs, pos_bs_els, pos_ris, pos_ris_els, pos_ues)

        # Store LoS channels
        los_bs_ue_channels[ss, :] = np.squeeze(bs_ue)
        los_ris_ue_channels[ss, :] = np.squeeze(ris_ue)

    print('--------------------')
    print('BS-UE power LoS stats')
    print('--------------------')
    power_bs_ue = np.abs(los_bs_ue_channels)**2

    print('Mean = {:.2e}'.format(power_bs_ue.mean()))
    print('Max. = {:.2e}'.format(power_bs_ue.max()))
    print('Min. = {:.2e}'.format(power_bs_ue.min()))

    print('--------------------')
    print('RIS-UE power LoS stats')
    print('--------------------')
    power_ris_ue = np.abs(los_ris_ue_channels)**2

    print('Mean = {:.2e}'.format(power_ris_ue.mean()))
    print('Max. = {:.2e}'.format(power_ris_ue.max()))
    print('Min. = {:.2e}'.format(power_ris_ue.min()))

    print('\n')

    print('Values are given per antennas / element.')

breakpoint()
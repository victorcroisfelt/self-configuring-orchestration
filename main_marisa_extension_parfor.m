% close all
% clear all
% 
% addpath('functions')
% 
% plot_flag = false();

%% Physical Parameters
c = 3*10^8;   %physconst('LightSpeed'); % speed of light
freq = 28*10^9; % [Hz] signal frequency
lambda = c/freq; % lambda [m]

%% Channel Parameters

% Pathloss exponent
beta = 2;

% NLoS variance 
nlos_var = 1;

% Noise power
sigma2n = 10^((-94 - 30)/10);

%% BS Parameters

% Number of BS antennas
M = 64;

% Transmit power at the BS = 20 dBm (-26)
P_bs = 10^((20 - 30)/10);

% RIS/BSantenna elements interdistance normalized wrt lambda
el_dist = 0.5;

%% HRIS Parameters

% Number of RIS elements on the x-axis
Nx = 32;

% Number of RIS elements on the y-axis
Ny = 1;

% Number of RIS elements
N = Nx * Ny;

% HRIS absorption parameter
eta = 0.5;

%% UE Parameters

% Number of UEs
K = 16;

% Transmit power at the UE = 10 dBm
P_ue = 10^((10 - 30)/10);

%% System Parameters

% Channel estimation relative length
L = 64;

% Channel estimation length
tau_est = L * K;

% Coherence interval length
tau_c = 2 * tau_est;

% Communication length
tau_com = tau_est;

%% Simulation Parameters

% Range of probing relative lenght
C_vec = 2.^(1:log2(L)); % iterate

% Range of probing length
tau_pro = C_vec * K;

% Range of reflection length
tau_ref = tau_c - tau_pro;

% Number of setups (simulation runs)
N_setup = 50;
N_channel_realizations = 50;

% Parameterized variable to define area of interest
scenario_size = 100;

%% Simulation Setup

% Limits of the area of the interest
x_lim = [0, scenario_size];
y_lim = [-scenario_size, scenario_size];

% Coordinates of RIS
pos_ris = 0 + 1j * y_lim(1);

% Coordinates of RIS elements
pos_ris_els = pos_ris + lambda./2. * (0:N-1);         
pos_ris_els = pos_ris_els - ((pos_ris_els(end) - pos_ris_els(1))/2);

% BS coordinates
pos_bs = -scenario_size + 1j * 0;

% Coordinates of BS elements
pos_bs_els = pos_bs + lambda./2. * (0:M-1);         
pos_bs_els = pos_bs_els - ((pos_bs_els(end) - pos_bs_els(1))/2);

% Pilots
pilots = sqrt(K) * diag(ones(K, 1));

% False_alarm values
false_alarm_prob_vec = [0.001, 0.01, 0.1];

%% Results
hat_detected_ue_pow = zeros(numel(C_vec), numel(false_alarm_prob_vec), N_channel_realizations, N_setup);
hat_detected_ue_sig = zeros(numel(C_vec), numel(false_alarm_prob_vec), N_channel_realizations, N_setup);

th_detected_ue_pow = zeros(numel(C_vec), numel(false_alarm_prob_vec), N_channel_realizations, N_setup);
th_detected_ue_sig = zeros(numel(C_vec), numel(false_alarm_prob_vec), N_channel_realizations, N_setup);

distance_pow = zeros(numel(C_vec), numel(false_alarm_prob_vec), N_channel_realizations, N_setup);
distance_sig = zeros(numel(C_vec), numel(false_alarm_prob_vec), N_channel_realizations, N_setup);

MSE_cha_est_pow = zeros(numel(C_vec), numel(false_alarm_prob_vec), N_channel_realizations, N_setup);
MSE_cha_est_sig = zeros(numel(C_vec), numel(false_alarm_prob_vec), N_channel_realizations, N_setup);

SE_est_sig = zeros(K, numel(C_vec), numel(false_alarm_prob_vec), N_channel_realizations, N_setup);
SINR_est_sig = zeros(K, numel(C_vec), numel(false_alarm_prob_vec), N_channel_realizations, N_setup);

SE_est_pow = zeros(K, numel(C_vec), numel(false_alarm_prob_vec), N_channel_realizations, N_setup);
SINR_est_pow = zeros(K, numel(C_vec), numel(false_alarm_prob_vec), N_channel_realizations, N_setup);

SE_bound_sig = zeros(K, numel(C_vec), numel(false_alarm_prob_vec), N_channel_realizations, N_setup);
SINR_bound_sig = zeros(K, numel(C_vec), numel(false_alarm_prob_vec), N_channel_realizations, N_setup);

SE_bound_pow = zeros(K, numel(C_vec), numel(false_alarm_prob_vec), N_channel_realizations, N_setup);
SINR_bound_pow = zeros(K, numel(C_vec), numel(false_alarm_prob_vec), N_channel_realizations, N_setup);

%% Simulation

rng(49)

for ind_setup = 1:N_setup
    
    message = ['run', num2str(ind_setup), '/',num2str(N_setup)];
    disp(message)

    tic

    % Drop UEs
    pos_ues = (x_lim(2)-x_lim(1)) * rand(K, 1) + x_lim(1); 
    pos_ues = pos_ues + 1j * ((y_lim(2)-y_lim(1)) * rand(K, 1) + y_lim(1));
   
    % Compute LOS components 
    los_components(lambda, M, N, K, pos_bs, pos_bs_els, pos_ris, pos_ris_els, pos_ues)



    % Prepare to save parfor results
    par_hat_detected_ue_pow = zeros(numel(C_vec), numel(false_alarm_prob_vec), N_channel_realizations);
    par_hat_detected_ue_sig = zeros(numel(C_vec), numel(false_alarm_prob_vec), N_channel_realizations);

    par_th_detected_ue_pow = zeros(numel(C_vec), numel(false_alarm_prob_vec), N_channel_realizations);
    par_th_detected_ue_sig = zeros(numel(C_vec), numel(false_alarm_prob_vec), N_channel_realizations);

    par_distance_pow = zeros(numel(C_vec), numel(false_alarm_prob_vec), N_channel_realizations);
    par_distance_sig = zeros(numel(C_vec), numel(false_alarm_prob_vec), N_channel_realizations);

    par_MSE_cha_est_pow = zeros(numel(C_vec), numel(false_alarm_prob_vec), N_channel_realizations);
    par_MSE_cha_est_sig = zeros(numel(C_vec), numel(false_alarm_prob_vec), N_channel_realizations);

    par_SE_est_sig = zeros(K, numel(C_vec), numel(false_alarm_prob_vec), N_channel_realizations);
    par_SINR_est_sig = zeros(K, numel(C_vec), numel(false_alarm_prob_vec), N_channel_realizations);

    par_SE_est_pow = zeros(K, numel(C_vec), numel(false_alarm_prob_vec), N_channel_realizations);    
    par_SINR_est_pow = zeros(K, numel(C_vec), numel(false_alarm_prob_vec), N_channel_realizations);

    par_SE_bound_sig = zeros(K, numel(C_vec), numel(false_alarm_prob_vec), N_channel_realizations);    
    par_SINR_bound_sig = zeros(K, numel(C_vec), numel(false_alarm_prob_vec), N_channel_realizations);

    par_SE_bound_pow = zeros(K, numel(C_vec), numel(false_alarm_prob_vec), N_channel_realizations);
    par_SINR_bound_pow = zeros(K, numel(C_vec), numel(false_alarm_prob_vec), N_channel_realizations);


    % Number of simulation points to go over
    length_Cvec = length(C_vec);
    length_fa = length(false_alarm_prob_vec);

    %parfor
    for ind_ch =1:N_channel_realizations
        
        % Compute blockage probability of paths
        [block_u_0] = blockage_path(d_u_0, lambda_b, radius_B, height_B, height_U, hheight_A); % RIS block
        [block_Bu] = blockage_path(d_Bu, lambda_b, radius_B, height_B, height_U, hheight_A); % BS block
        
        block_u_0(:) = 0;

        % Channels
        [a_R_los, ~, G_los, h_los, h_D_los, a_BSU_los, a_BSR_los] = compute_channels(d_0_G, d_u_0, d_Bu, [beta,beta_nlos], N, M, phiB_0_a, phiB_0_d, phiU_0, phiUB,el_dist, block_u_0, block_Bu);

        for ind_d = 1:length_Cvec

            D = C_vec(ind_d);
            C = D;
            
            %message = ['run', num2str(ind_setup), '/',num2str(N_setup),' Chan',num2str(ind_ch),'/',num2str(N_channel_realizations),' C',num2str(ind_d),'/',num2str(length(C_vec))];
            %disp(message)
            
            % Generate configuration codebook according to Eq. (12)
            [~, theta_in, phi_vec] = generate_codebook(C, N);
            
            % Probing codebook definition
            Theta_prob_sig = repmat(eye(N), [1, 1, C]);
            Theta_prob_pow = zeros(N, N, C);
            for t = 1:C
                Theta_prob_pow(:, :, t) = diag(theta_in(:, t));
            end
            
            % Received signal at RIS
            [Y_r_pow] = received_signal_RIS(Theta_prob_pow, h_los, eta, pilots, P_ue, sigma2n);
            [Y_r_sig] = received_signal_RIS(Theta_prob_sig, h_los, eta, pilots, P_ue, sigma2n);

            for ind_prob = 1:length_fa
   
                false_alarm_prob = false_alarm_prob_vec(ind_prob);

                [Theta_opt_sig, hat_det_rate_sig, th_det_rate_sig, dist_sig] = MARISA_EXTENSION(Y_r_sig, theta_in, false_alarm_prob, phiB_0_a, sigma2n, 'signal');
                [Theta_opt_pow, hat_det_rate_pow, th_det_rate_pow, dist_pow] = MARISA_EXTENSION(Y_r_pow, theta_in, false_alarm_prob, phiB_0_a, sigma2n, 'power');

                % Equivalent BS-UE channel
                Theta_over_blocks_pow = cat(3, Theta_prob_pow, repmat(Theta_opt_pow,[1,1,L-C])); % to change after hris optimization
                                                            
                [SINR_pow, SE_pow, SINR_b_pow, SE_b_pow] = channel_estimation_MMIMO(Theta_prob_pow, Theta_opt_pow, M, C, L, K, tau_est, tau_c, sigma2n, P_ue, G_los, h_los, h_D_los, eta);
                [SINR_sig, SE_sig, SINR_b_sig, SE_b_sig] = channel_estimation_MMIMO(Theta_prob_sig, Theta_opt_sig, M, C, L, K, tau_est, tau_c, sigma2n, P_ue, G_los, h_los, h_D_los, eta);

                % Channel estimation mse
                [MSE_sig] = channel_estimation_MSE(M, C, L, K, sigma2n, P_ue, G_los, Theta_prob_sig, Theta_opt_sig, h_los, h_D_los, eta);
                [MSE_pow] = channel_estimation_MSE(M, C, L, K, sigma2n, P_ue, G_los, Theta_prob_pow, Theta_opt_pow, h_los, h_D_los, eta);

                % Save simulation results
                par_hat_detected_ue_sig(ind_d, ind_prob, ind_ch) = hat_det_rate_sig;
                par_hat_detected_ue_pow(ind_d, ind_prob, ind_ch) = hat_det_rate_pow;
                
                par_th_detected_ue_sig(ind_d, ind_prob, ind_ch) = th_det_rate_sig;
                par_th_detected_ue_pow(ind_d, ind_prob, ind_ch) = th_det_rate_pow;

                par_distance_sig(ind_d, ind_prob, ind_ch) = dist_sig;
                par_distance_pow(ind_d, ind_prob, ind_ch) = dist_pow;

                par_MSE_cha_est_sig(ind_d, ind_prob, ind_ch) = MSE_sig;
                par_MSE_cha_est_pow(ind_d, ind_prob, ind_ch) = MSE_pow;
                
                par_SINR_est_sig(:, ind_d, ind_prob, ind_ch) = SINR_sig;
                par_SINR_est_pow(:, ind_d, ind_prob, ind_ch) = SINR_pow;
                
                par_SE_est_sig(:, ind_d, ind_prob, ind_ch) = SE_sig;
                par_SE_est_pow(:, ind_d, ind_prob, ind_ch) = SE_pow;

                par_SINR_bound_sig(:, ind_d,ind_prob, ind_ch) = SINR_b_sig;
                par_SINR_bound_pow(:, ind_d,ind_prob, ind_ch) = SINR_b_pow;
                
                par_SE_bound_sig(:, ind_d,ind_prob, ind_ch) = SE_b_sig;
                par_SE_bound_pow(:, ind_d,ind_prob, ind_ch) = SE_b_pow;
                
            end
        end
    end

    % Save simulation results
    hat_detected_ue_sig(:, :, :, ind_setup) = par_hat_detected_ue_sig;
    hat_detected_ue_pow(:, :, :, ind_setup) = par_hat_detected_ue_pow;
    
    th_detected_ue_sig(:, :, :, ind_setup) = par_th_detected_ue_sig;
    th_detected_ue_pow(:, :, :, ind_setup) = par_th_detected_ue_pow;

    distance_sig(:, :, :, ind_setup) = par_distance_sig;
    distance_pow(:, :, :, ind_setup) = par_distance_pow;    
    
    MSE_cha_est_sig(:, :, :, ind_setup) = par_MSE_cha_est_sig;
    MSE_cha_est_pow(:, :, :, ind_setup) = par_MSE_cha_est_pow;
    
    SINR_est_sig(:, :, :, :, ind_setup) = par_SINR_est_sig;
    SINR_est_pow(:, :, :, :, ind_setup) = par_SINR_est_pow;
    
    SE_est_sig(:, :, :, :, ind_setup) = par_SE_est_sig;
    SE_est_pow(:, :, :, :, ind_setup) = par_SE_est_pow;

    SINR_bound_sig(:, :, :, :, ind_setup) = par_SINR_bound_sig;
    SINR_bound_pow(:, :, :, :, ind_setup) = par_SINR_bound_pow;
    
    SE_bound_sig(:, :, :, :, ind_setup) = par_SE_bound_sig;
    SE_bound_pow(:, :, :, :, ind_setup) = par_SE_bound_pow;

    toc

end


%% Save results
save(['RESULTS','_M', num2str(M),'_N', num2str(N),'_K', num2str(K),'_L', num2str(L),'_',num2str(scenario_size),'x',num2str(scenario_size)])



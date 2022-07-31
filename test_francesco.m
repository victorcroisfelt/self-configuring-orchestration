%close all
%clear all

%addpath('functions')

%plot_flag = false();

%% Physical Parameters
c = 3*10^8;   %physconst('LightSpeed'); % speed of light
freq = 28*10^9; % [Hz] signal frequency
lambda = c/freq; % lambda [m]

%% Channel Parameters

% Blockage parameters
lambda_b = 1.0; % density of pedestrians
radius_B = 0.6; % [m] radius of pedestrians
height_B = 1.7; % [m] height of pedestrians
height_U = 1.5; % [m] height of UEs
hheight_A = 6;  % [m] height of NR

% Pathloss exponent
beta = 2; % LoS
beta_nlos = 4; % NLoS

% Noise power (-80)
sigma2n = 10^((-94-30)/10);

%% BS Parameters

% Number of BS antennas
%M = 64;

% Transmit power at the BS = 20 dBm (-26)
P = 10^((20-30)/10);

% RIS/BSantenna elements interdistance normalized wrt lambda
el_dist = 0.5;

%% HRIS Parameters

% Number of RIS elements on the x-axis
%Nx = 32;

% Number of RIS elements on the y-axis
Ny = 1;

% Number of RIS elements
N = Nx * Ny;

% HRIS absorption parameter
eta = 0.8;

%% UE Parameters

% Number of UEs
%K = 16;

% Transmit power at the UE = 20 dBm (-26)
P_ue = 10^((20-30)/10);

%% System Parameters

% Channel estimation relative length
%L = 64;

% Channel estimation length
tau_est = L * K;

% Coherence interval length
tau_c = 2 * tau_est;

% Communication length
tau_com = tau_est;

%% Simulation Parameters

% Range of probing relative lenght
C_vec = 2.^(1:6);     % iterate

% Range of probing length
tau_pro = C_vec * K;

% Range of reflection length
tau_ref = tau_c - tau_pro;

% Number of setups (simulation runs)
N_setup = 50;
N_channel_realizations = 50;

% Parameterized variable to define area of interest
% scenario_size = 50;

%% Simulation Setup

% Limits of the area of the interest
x_lim = [0, scenario_size];
y_lim = [-scenario_size, scenario_size];

% RISs
ris_0 = [0;y_lim(1)];                                        % Coordinates of RIS
ris_0_el = ris_0+lambda./2.*cat(1, 0:N-1, zeros(1, N));         % Coordinates of RIS elements

% BS coordinates
bs = [-scenario_size;0];

% Geometry (downlink)
phiB_0_a = atan2(bs(2)-ris_0(2), bs(1)-ris_0(1)); %135 deg Angle of Arrival BS-RIS
phiB_0_d = atan2(ris_0(2)-bs(2), ris_0(1)-bs(1)); %-45 deg Angle of Departure BS-RIS
d_0_G = norm(bs-ris_0);                          % Distance BS-RIS

% pilots
pilots = diag(ones(K, 1));

% false_alarm
false_alarm_prob_vec = [0.001, 0.01, 0.1];

%% Results
hat_detected_ue_sig = zeros(numel(C_vec), numel(false_alarm_prob_vec), N_setup, N_channel_realizations);
hat_detected_ue_pow = zeros(numel(C_vec), numel(false_alarm_prob_vec), N_setup, N_channel_realizations);

th_detected_ue_sig = zeros(numel(C_vec), numel(false_alarm_prob_vec), N_setup, N_channel_realizations);
th_detected_ue_pow = zeros(numel(C_vec), numel(false_alarm_prob_vec), N_setup, N_channel_realizations);

MSE_cha_est_sig = zeros(numel(C_vec), numel(false_alarm_prob_vec), N_setup, N_channel_realizations);
MSE_cha_est_pow = zeros(numel(C_vec), numel(false_alarm_prob_vec), N_setup, N_channel_realizations);

SINR_est_sig = zeros(K, numel(C_vec), numel(false_alarm_prob_vec), N_setup, N_channel_realizations);
SE_est_sig = zeros(K, numel(C_vec), numel(false_alarm_prob_vec), N_setup, N_channel_realizations);

SINR_est_pow = zeros(K, numel(C_vec), numel(false_alarm_prob_vec), N_setup, N_channel_realizations);
SE_est_pow = zeros(K, numel(C_vec), numel(false_alarm_prob_vec), N_setup, N_channel_realizations);

%% Simulation

rng(49)

for ind_setup = 1:N_setup
    
    ue_s = [(x_lim(2)-x_lim(1))*rand(1, K, N_setup) + x_lim(1); (y_lim(2)-y_lim(1))*rand(1, K, N_setup)+y_lim(1)]; % UEs coordinates --- drawn in x between [x_lim(1), x_lim(2)], y between [y_lim(1), y_lim(2)]
    
    %ue = [(x_lim(2)-x_lim(1))*rand(1,K)+x_lim(1); (y_lim(2)-y_lim(1))*rand(1,K)+y_lim(1)];
    ue = ue_s(:, 1:K, ind_setup);
    
    % Compute geometry wrt BS
    rot_angle = 90/180*pi;
    [bs_rot, ue_rot, ~] = rotate_geometry(rot_angle, bs, ue, ris_0_el, false());
    
    % Compute UE angles wrt BS and RIS
    phiUB = atan2(ue_rot(2, :) - bs_rot(2), ue_rot(1, :)-bs_rot(1))';
    phiU_0 = atan2(ue(2, :)-ris_0(2), ue(1, :)-ris_0(1))';
    
    % Compute UE distances wrt BS and RIS
    d_Bu = zeros(K, 1);
    d_u_0 = zeros(K, 1);
    
    for k=1:K

        % Distance BS-UE
        d_Bu(k) = norm(bs - ue(:, k));   

        % Distance RIS-UE
        d_u_0(k) = norm(ris_0-ue(:, k));   

    end
    
    for ind_ch =1:N_channel_realizations
        
        % Compute blockage probability of paths
        [block_u_0] = blockage_path(d_u_0, lambda_b, radius_B, height_B, height_U, hheight_A); % RIS block
        [block_Bu] = blockage_path(d_Bu, lambda_b, radius_B, height_B, height_U, hheight_A); % BS block
        
        % Channels
        [a_R_los, ~, G_los, h_los, h_D_los, a_BSU_los, a_BSR_los] = compute_channels(d_0_G, d_u_0, d_Bu, [beta,beta_nlos], N, M, phiB_0_a, phiB_0_d, phiU_0, phiUB,el_dist, block_u_0, block_Bu);
        
        % Channel realizations
        %h_los = h_los .* exp(1i.*rand(1, K).*2*pi);
        %h_D_los = h_D_los .* exp(1i.*rand(1, K).*2*pi);
       
        for ind_d = 1:length(C_vec)

            D = C_vec(ind_d);
            C = D;
            
            message = ['run', num2str(ind_setup), '/',num2str(N_setup),' Chan',num2str(ind_ch),'/',num2str(N_channel_realizations),' C',num2str(ind_d),'/',num2str(length(C_vec))];
            disp(message)
            
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
            
            false_alarm_prob = 0.001;
            for ind_prob = 1:numel(false_alarm_prob_vec)
                false_alarm_prob = false_alarm_prob_vec(ind_prob);

                [Theta_opt_sig, hat_det_rate_sig, th_det_rate_sig] = MARISA_EXTENSION(Y_r_pow, theta_in, false_alarm_prob, phiB_0_a, sigma2n, 'signal');
                [Theta_opt_pow, hat_det_rate_pow, th_det_rate_pow] = MARISA_EXTENSION(Y_r_pow, theta_in, false_alarm_prob, phiB_0_a, sigma2n, 'power');

                % Equivalent BS-UE channel
                Theta_over_blocks_pow = cat(3, Theta_prob_pow, repmat(Theta_opt_pow,[1,1,L-C])); % to change after hris optimization
                [SINR_pow, SE_pow]   = channel_estiamation_MMIMO(Theta_prob_pow, Theta_opt_pow, M, C, L, K, tau_est, tau_c, sigma2n,P_ue, G_los, h_los, h_D_los, eta);
                [SINR_sig, SE_sig]   = channel_estiamation_MMIMO(Theta_prob_sig, Theta_opt_sig, M, C, L, K, tau_est, tau_c, sigma2n,P_ue, G_los, h_los, h_D_los, eta);
                
                
                % Channel estimation mse
                [MSE_sig] = channel_estiamation_MSE(M,L,K,sigma2n,P_ue,G_los,Theta_prob_sig, Theta_opt_sig, h_los, eta);
                [MSE_pow] = channel_estiamation_MSE(M,L,K,sigma2n,P_ue,G_los,Theta_prob_pow, Theta_opt_pow, h_los, eta);
                
                % Save simulation results
                hat_detected_ue_sig(ind_d,ind_prob,ind_setup,ind_ch) = hat_det_rate_sig;
                hat_detected_ue_pow(ind_d,ind_prob,ind_setup,ind_ch) = hat_det_rate_pow;
                
                th_detected_ue_sig(ind_d,ind_prob,ind_setup,ind_ch) = th_det_rate_sig;
                th_detected_ue_pow(ind_d,ind_prob,ind_setup,ind_ch) = th_det_rate_pow;

                MSE_cha_est_sig(ind_d,ind_prob,ind_setup,ind_ch) = MSE_sig;
                MSE_cha_est_pow(ind_d,ind_prob,ind_setup,ind_ch) = MSE_pow;
                
                SINR_est_sig(:,ind_d,ind_prob,ind_setup,ind_ch) = SINR_sig;
                SINR_est_pow(:,ind_d,ind_prob,ind_setup,ind_ch) = SINR_pow;
                
                SE_est_sig(:,ind_d,ind_prob,ind_setup,ind_ch) = SE_sig;
                SE_est_pow(:,ind_d,ind_prob,ind_setup,ind_ch) = SE_pow;
                
                
            end
        end
    end
end


%% Save results
save(['RESULTS','_M', num2str(M),'_N', num2str(N),'_K', num2str(K),'_L', num2str(L),'_',num2str(scenario_size),'x',num2str(scenario_size)])



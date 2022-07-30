close all
clear all

addpath('functions')

plot_flag = false();

%% Physical Parameters
c = physconst('LightSpeed'); % speed of light
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
M = 64;                             

% Transmit power at the BS = 20 dBm (-26)
P = 10^((20-30)/10);   

% RIS/BSantenna elements interdistance normalized wrt lambda
el_dist = 0.5;                      

%% HRIS Parameters

% Number of RIS elements on the x-axis
Nx = 32;         

% Number of RIS elements on the y-axis
Ny = 1;   

% Number of RIS elements
N = Nx * Ny;   

% Nr codebook elements
M_cod = 64;          

% HRIS absorption parameter
eta = 0.8;                          

%% UE Parameters

% Number of UEs
K_vec = 16;

% Transmit power at the UE = 20 dBm (-26)
P_ue = 10^((20-30)/10);             

%% System Parameters

% Channel estimation relative length
L = 64; 

% Channel estimation length
tau_est = L * K_vec;

% Coherence interval length
tau_c = 2 * tau_est;    

% Communication length
tau_com = tau_est;

%% Simulation parameters

% Range of probing relative lenght
C_vec = 2.^(1:6);     % iterate

% Range of probing length
tau_pro = C_vec * K_vec;

% Range of reflection length
tau_ref = tau_c - tau_pro;

% Number of setups (simulation runs)
N_setup = 50;
N_channel_realizatios = 50;

% Don't know what is that
I = 1;  % risma iterations
Q_bit = 1;

%% Simulation
extr_scenario_vec = 50;

for ind_scen = 1:length(extr_scenario_vec)

    x_lim = [0, extr_scenario_vec(ind_scen)];
    y_lim = [-extr_scenario_vec(ind_scen), extr_scenario_vec(ind_scen)];
    
    % RISs
    ris_0 = [0;y_lim(1)];                                        % Coordinates of RIS
    ris_0_el = ris_0+lambda./2.*cat(1,0:N-1,zeros(1,N));         % Coordinates of RIS elements
    
    bs = [-extr_scenario_vec(ind_scen);0];
    
    % Geometry (downlink)
    phiB_0_a = atan2(bs(2)-ris_0(2),bs(1)-ris_0(1)); %135 deg Angle of Arrival BS-RIS
    phiB_0_d = atan2(ris_0(2)-bs(2),ris_0(1)-bs(1)); %-45 deg Angle of Departure BS-RIS
    d_0_G = norm(bs-ris_0);                          % Distance BS-RIS
    
    
    rng(49)
    for run_ind = 1:N_setup
    
    ue_s = [(x_lim(2)-x_lim(1))*rand(1,max(K_vec),N_setup)+x_lim(1); (y_lim(2)-y_lim(1))*rand(1,max(K_vec),N_setup)+y_lim(1)];
    
    %channel realization for 
    
    random_phase = exp(1i.*rand(max(K_vec),N_setup).*2*pi);
    


        for ind_k = 1:length(K_vec)

            for ind_d = 1:length(C_vec)

                D = C_vec(ind_d);
                    
                message = ['run', num2str(run_ind), '/',num2str(N_setup),' n_ue',num2str(ind_k),'/',num2str(length(K_vec)),' scenario',num2str(ind_scen),'/',num2str(length(extr_scenario_vec)),' D',num2str(ind_d),'/',num2str(length(D_vec))];
                disp(message)
                K = K_vec(ind_k);
                
                %ue = [(x_lim(2)-x_lim(1))*rand(1,K)+x_lim(1); (y_lim(2)-y_lim(1))*rand(1,K)+y_lim(1)];          % UEs coordinates --- drawn from [- x_max, x_max] x [0, y_max]
                ue = ue_s(:,1:K,run_ind);
                
                
                % Compute geometry wrt BS
                rot_angle = 90/180*pi;
                [bs_rot, ue_rot, ~] = rotate_geometry(rot_angle, bs, ue, ris_0_el, false());
                
                % Compute UE angles wrt BS and RIS
                phiUB = atan2(ue_rot(2,:)-bs_rot(2),ue_rot(1,:)-bs_rot(1))';
                phiU_0 = atan2(ue(2,:)-ris_0(2),ue(1,:)-ris_0(1))';
                
                % Compute UE distances wrt BS and RIS
                d_Bu = zeros(K,1);
                d_u_0 = zeros(K,1);
        
                for k=1:K
                    d_Bu(k) = norm(bs-ue(:,k));                                 % Distance BS-UE
                    d_u_0(k) = norm(ris_0-ue(:,k));                                % Distance RIS-UE
                end
                
                %% compute blockage probability of paths
                [block_u_0] = blockage_path(d_u_0, lambda_b, radius_B, height_B, height_U, hheight_A); % ris block
                [block_Bu] = blockage_path(d_Bu, lambda_b, radius_B, height_B, height_U, hheight_A);    % bs  block
                
                %block_u_0(:) = 0;       % remove blockages, to have channels all in LoS
                %block_Bu(:) = 0;
                
                %% Channels
                [a_R_0_los,~,G_0_los,h_0_los,h_D_los,a_BSU_los,a_BSR_0_los] = compute_channels(d_0_G,d_u_0,d_Bu,[beta,beta_nlos],N,M,phiB_0_a,phiB_0_d,phiU_0,phiUB,el_dist, block_u_0, block_Bu);
                
                
                %% pilots
                pilots = diag(ones(K,1));
                
                %% frame
                C = D;

                [~,theta_in, phi_vec] = generate_codebook(C,N);
                
                % Probing codebook definition
                Theta_prob_sig = repmat(eye(N), [1,1,C]);
                Theta_prob_pow = zeros(N,N,C);
                for t = 1:C
                    Theta_prob_pow(:,:,t) = diag(theta_in(:,t));
                end
                
                % Received signal at RIS
                [Y_r_pow] = received_signal_RIS(Theta_prob_pow, h_0_los, eta, pilots, P_ue, sigma2n);
                [Y_r_sig] = received_signal_RIS(Theta_prob_sig, h_0_los, eta, pilots, P_ue, sigma2n);
                
                threshold = 1*10^-7;
                [Theta_opt_sig] = MARISA_EXTENSION(Y_r_pow, theta_in, threshold, phiB_0_a, sigma2n, 'signal');
                [Theta_opt_pow] = MARISA_EXTENSION(Y_r_pow, theta_in, threshold, phiB_0_a, sigma2n, 'power');
                
                % Equivalent BS-UE channel
                Theta_over_blocks = cat(3, Theta_prob_pow,Theta_prob_sig); % to change after hris optimization
                [H_circ] = equivalent_BS_UE_channel(Theta_over_blocks, h_0_los, G_0_los, h_D_los, eta);
                [Y_b]    = received_signal_BS(H_circ, pilots, P_ue, sigma2n);
                
            end
        end
    end
end



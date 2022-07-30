close all
clear all

addpath('functions')

%% Multi-user scenario

%% Parameters
plot_flag = false();

% Settings
K_vec = (5)';                                                       % Nr UEs
extr_scenario_vec = 50;


% Physical constant
c = 3*10^8;                         % Light speed

%blockage
lambda_b = 1;%0.1, 0.3, 0.5, 1.0     % density of pedestrians
radius_B = 0.6;                     % [m] radius of pedestrians
height_B = 1.7;                     % [m] height of pedestrians
height_U = 1.5;                     % [m] height of UEs
hheight_A = 6;                      % [m] height of NR

% simulation parameters
I = 1;                              % risma iterations
N_run = 1;                         % simulation runs
eta = 0.8;                          %HRIS absorption parameter

% 
freq = 28*10^9;                     % [Hz] signal frequency
lambda = c/freq;                    % lambda
M_cod = 64;                        % Nr codebook elements
el_dist = 0.5;                      % RIS/BSantenna elements interdistance normalized wrt lambda


Q_bit = 1;

% Radio
M = 4;                              % Number of BS antennas
Nx = 32;                            % Number of RIS elements on the x axis
Ny = 1;                             % Number of RIS elements on the y axis
N = Nx*Ny;                          % Number of RIS elements
beta = 2;                           % Pathloss exponent
beta_nlos = 4;                      % Pathloss exponent
P = 10^((26-30)/10);                % Transmit power at the BS = 26 dBm
sigma2n = 10^((-80-30)/10);         % Noise power


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

direct_factor = 1;

rng(49)
ue_s = [(x_lim(2)-x_lim(1))*rand(1,max(K_vec),N_run)+x_lim(1); (y_lim(2)-y_lim(1))*rand(1,max(K_vec),N_run)+y_lim(1)];
random_phase = exp(1i.*rand(max(K_vec),N_run).*2*pi);

for run_ind = 1:N_run    
    for k_ind = 1:length(K_vec)
        
        message = ['run', num2str(run_ind), '/',num2str(N_run),' n_ue',num2str(k_ind),'/',num2str(length(K_vec)),' scenario',num2str(ind_scen),'/',num2str(length(extr_scenario_vec))];
        disp(message)
        K = K_vec(k_ind);
        
        %ue = [(x_lim(2)-x_lim(1))*rand(1,K)+x_lim(1); (y_lim(2)-y_lim(1))*rand(1,K)+y_lim(1)];          % UEs coordinates --- drawn from [- x_max, x_max] x [0, y_max]
        ue = ue_s(:,1:K,run_ind);
        
        
        %% Compute geometry wrt BS
        rot_angle = 90/180*pi;
        [bs_rot, ue_rot, ~] = rotate_geometry(rot_angle, bs, ue, ris_0_el, false());
        
        %% Compute UE angles wrt BS and RIS
        phiUB = atan2(ue_rot(2,:)-bs_rot(2),ue_rot(1,:)-bs_rot(1))';
        phiU_0 = atan2(ue(2,:)-ris_0(2),ue(1,:)-ris_0(1))';
        
        %% Compute UE distances wrt BS and RIS
        d_Bu = zeros(K,1);
        d_u_0 = zeros(K,1);

        for k=1:K
            d_Bu(k) = norm(bs-ue(:,k));                                 % Distance BS-UE
            d_u_0(k) = norm(ris_0-ue(:,k));                                % Distance RIS-UE
        end
        
        %% compute blockage probability of paths
        [block_u_0]  = blockage_path(d_u_0, lambda_b, radius_B, height_B, height_U, hheight_A); % ris block
        [block_Bu]   = blockage_path(d_Bu, lambda_b, radius_B, height_B, height_U, hheight_A);    % bs  block
        
        block_u_0(:) = 0;       % remove blockages, to have channels all in LoS
        block_Bu(:) = 0;
        
        %% Channels
        [a_R_0_los,~,G_0_los,h_0_los,h_D_los,a_BSU_los,a_BSR_0_los] = compute_channels(d_0_G,d_u_0,d_Bu,[beta,beta_nlos],N,M,phiB_0_a,phiB_0_d,phiU_0,phiUB,el_dist, block_u_0, block_Bu);
        
        h_D_los = h_D_los.*direct_factor;
        
        %% ORACLE MARISA estimation with LoS channels
        gain_ue_0_los = d_u_0.^-beta;
        gain_bs_los = d_0_G.^-beta;
        
        [~,C_in, phi_vec] = generate_codebook(M_cod,N);
                
        % ESTIMATION MARISA
        h_in = cat(1,h_0_los);
        G = cat(1,G_0_los);
        
        % Equivalent Channel        
        hd = h_D_los;
        Hb_marisa = zeros(N+1,M,K);
        for k=1:K
            Hb_marisa(:,:,k) = [sqrt(eta).*diag(h_in(:,k)')*G(:,:); hd(:,k)';];
        end

        H_MISO_eq = zeros(K, M, M_cod);
        H_MISO_eq_pow = zeros(K, M_cod);
        for c_ind = 1:M_cod
            Phim = diag(C_in(:,c_ind));
            H_MISO_eq(:,:,c_ind) = (h_in'*Phim*G + hd');
            for u = 1:K
                H_MISO_eq_pow(u,c_ind) = norm(H_MISO_eq(u,:,c_ind)');
            end 
        end

    end
end
end

figure
plot(phi_vec/pi*180, H_MISO_eq_pow')
legend('1', '2', '3', '4', '5')
xlabel('RIS beamformer direction $^{\circ}$','Interpreter', 'latex')
ylabel('$\|\mathbf{G}^{*} \mathbf{\Phi}^{*} \mathbf{h}  + \mathbf{h}_{d} \|_{F}$','Interpreter', 'latex')

plot_geometry(ris_0_el,ue,bs, true());


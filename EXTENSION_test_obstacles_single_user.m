close all
clear all

addpath('functions')

%% Multi-user scenario

%% Parameters
plot_flag = false();

% Settings
K_vec = (2:4:20)';                                                       % Nr UEs
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
N_run = 50;                         % simulation runs
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

% Virtual multipath component
delta_angle = 15/1800*pi; 


%% Results initialization
srate_marisa_oracle1=zeros(numel(K_vec),numel(extr_scenario_vec),N_run);
SMSE_marisa_oracle1 =zeros(numel(K_vec),numel(extr_scenario_vec),N_run);
MSE_marisa_oracle1  =zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
SINR_marisa_oracle1 =zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);

srate_marisa_oracle2=zeros(numel(K_vec),numel(extr_scenario_vec),N_run);
SMSE_marisa_oracle2 =zeros(numel(K_vec),numel(extr_scenario_vec),N_run);
MSE_marisa_oracle2  =zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
SINR_marisa_oracle2 =zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);

srate_risma  =zeros(numel(K_vec),numel(extr_scenario_vec),N_run);
SMSE_risma   =zeros(numel(K_vec),numel(extr_scenario_vec),N_run);
MSE_risma    =zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
SINR_risma   =zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);

srate_marisa_oracle1_nlos=zeros(numel(K_vec),numel(extr_scenario_vec),N_run);
SMSE_marisa_oracle1_nlos =zeros(numel(K_vec),numel(extr_scenario_vec),N_run);
MSE_marisa_oracle1_nlos  =zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
SINR_marisa_oracle1_nlos =zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);

srate_marisa_oracle2_nlos=zeros(numel(K_vec),numel(extr_scenario_vec),N_run);
SMSE_marisa_oracle2_nlos =zeros(numel(K_vec),numel(extr_scenario_vec),N_run);
MSE_marisa_oracle2_nlos  =zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
SINR_marisa_oracle2_nlos =zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);

srate_risma_nlos  =zeros(numel(K_vec),numel(extr_scenario_vec),N_run);
SMSE_risma_nlos   =zeros(numel(K_vec),numel(extr_scenario_vec),N_run);
MSE_risma_nlos    =zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
SINR_risma_nlos   =zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);


Direct_path_gain = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Reflected_path_gain_marisa_oracle1 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Reflected_path_gain_marisa_oracle2 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Reflected_path_gain_marisa_all1 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Reflected_path_gain_marisa_all2 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Reflected_path_gain_marisa_peak1 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Reflected_path_gain_marisa_peak2 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Reflected_path_gain_marisa_fpeak1 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Reflected_path_gain_marisa_fpeak2 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Reflected_path_gain_risma = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);


Reflected_path_power_marisa_oracle1 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Reflected_path_power_marisa_oracle2 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Reflected_path_power_marisa_all1 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Reflected_path_power_marisa_all2 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Reflected_path_power_marisa_peak1 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Reflected_path_power_marisa_peak2 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Reflected_path_power_marisa_fpeak1 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Reflected_path_power_marisa_fpeak2 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Reflected_path_power_risma = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);

Direct_path_power_marisa_oracle1 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Direct_path_power_marisa_oracle2 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Direct_path_power_marisa_all1 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Direct_path_power_marisa_all2 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Direct_path_power_marisa_peak1 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Direct_path_power_marisa_peak2 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Direct_path_power_marisa_fpeak1 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Direct_path_power_marisa_fpeak2 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Direct_path_power_risma = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);

Reflected_path_tx_power_marisa_oracle1 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Reflected_path_tx_power_marisa_oracle2 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Reflected_path_tx_power_marisa_all1 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Reflected_path_tx_power_marisa_all2 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Reflected_path_tx_power_marisa_peak1 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Reflected_path_tx_power_marisa_peak2 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Reflected_path_tx_power_marisa_fpeak1 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Reflected_path_tx_power_marisa_fpeak2 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Reflected_path_tx_power_risma = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);

Direct_path_power_tx_marisa_oracle1 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Direct_path_power_tx_marisa_oracle2 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Direct_path_power_tx_marisa_all1 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Direct_path_power_tx_marisa_all2 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Direct_path_power_tx_marisa_peak1 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Direct_path_power_tx_marisa_peak2 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Direct_path_power_tx_marisa_fpeak1 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Direct_path_power_tx_marisa_fpeak2 = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);
Direct_path_power_tx_risma = zeros(max(K_vec),numel(K_vec),numel(extr_scenario_vec),N_run);

ue_s_all = zeros(2,max(K_vec),N_run);

for ind_scen = 1:length(extr_scenario_vec)
x_lim = [0, extr_scenario_vec(ind_scen)];
y_lim = [-extr_scenario_vec(ind_scen), extr_scenario_vec(ind_scen)];

% RISs
ris_0 = [0;y_lim(1)];                                        % Coordinates of RIS
ris_0_el = ris_0+lambda./2.*cat(1,0:N-1,zeros(1,N));         % Coordinates of RIS elements
ris_0_orientation = 0/180*pi;                                % Orientation of the RIS (that is the direction of the line where ris element lies)

ris_1 = [0;y_lim(2)];                                        % Coordinates of RIS
ris_1_el = ris_1+lambda./2.*cat(1,0:N-1,zeros(1,N));         % Coordinates of RIS elements
ris_1_orientation = 180/180*pi;                              % Orientation of the RIS (that is the direction of the line where ris element lies)


bs = [-extr_scenario_vec(ind_scen);0];
% Geometry (downlink)
phiB_0_a = atan2(bs(2)-ris_0(2),bs(1)-ris_0(1))+ris_0_orientation; %135 deg Angle of Arrival BS-RIS (we are in DL)
phiB_0_d = atan2(ris_0(2)-bs(2),ris_0(1)-bs(1))+ris_0_orientation; %-45 deg Angle of Departure BS-RIS (we are in DL)
d_0_G = norm(bs-ris_0);                          % Distance BS-RIS

phiB_1_a = atan2(bs(2)-ris_1(2),bs(1)-ris_1(1))+ris_1_orientation; %-135 deg Angle of Arrival BS-RIS  (we are in DL)
phiB_1_d = atan2(ris_1(2)-bs(2),ris_1(1)-bs(1))+ris_1_orientation; %  45 deg Angle of Departure BS-RIS
d_1_G = norm(bs-ris_1);                          % Distance BS-RIS


direct_factor = 1;

rng(49)
ue_s = [(x_lim(2)-x_lim(1))*rand(1,max(K_vec),N_run)+x_lim(1); (y_lim(2)-y_lim(1))*rand(1,max(K_vec),N_run)+y_lim(1)];

for run_ind = 1:N_run    
    for k_ind = 1:length(K_vec)
        
        message = ['run', num2str(run_ind), '/',num2str(N_run),' n_ue',num2str(k_ind),'/',num2str(length(K_vec)),' scenario',num2str(ind_scen),'/',num2str(length(extr_scenario_vec))];
        disp(message)
        K = K_vec(k_ind);
        
        ue = ue_s(:,1:K,run_ind);      % UEs coordinates --- drawn from [- x_max, x_max] x [0, y_max]
        
        %% Compute geometry wrt BS
        rot_angle = 90/180*pi;
        [bs_rot, ue_rot, ~] = rotate_geometry(rot_angle, bs, ue, ris_0_el, false());
        
        
        %% Compute UE angles wrt BS and RIS
        phiUB = atan2(ue_rot(2,:)-bs_rot(2),ue_rot(1,:)-bs_rot(1))';
        phiU_00 = atan2(ue(2,:)-ris_0(2),ue(1,:)-ris_0(1))' + ris_0_orientation;
        phiU_10 = atan2(ue(2,:)-ris_1(2),ue(1,:)-ris_1(1))' + ris_1_orientation;
        
        phiU_01 = atan2(ue(2,:)-ris_0(2),ue(1,:)-ris_0(1))' + ris_0_orientation + delta_angle; % Add virtual multipath component at the RIS
        phiU_11 = atan2(ue(2,:)-ris_1(2),ue(1,:)-ris_1(1))' + ris_0_orientation - delta_angle; % Add virtual multipath component at the RIS
        
        %% Compute UE distances wrt BS and RIS
        d_Bu = zeros(K,1);
        d_u_0 = zeros(K,1);
        d_u_1 = zeros(K,1);

        for k=1:K
            d_Bu(k) = norm(bs-ue(:,k));                                 % Distance BS-UE
            d_u_0(k) = norm(ris_0-ue(:,k));                                % Distance RIS-UE
            d_u_1(k) = norm(ris_1-ue(:,k));                                % Distance RIS-UE
        end
        
        %% compute blockage probability of paths
        [block_u_0]  = blockage_path(d_u_0, lambda_b, radius_B, height_B, height_U, hheight_A); % ris block
        [block_u_1]  = blockage_path(d_u_1, lambda_b, radius_B, height_B, height_U, hheight_A); % ris  block
        [block_Bu]   = blockage_path(d_Bu, lambda_b, radius_B, height_B, height_U, hheight_A);    % bs  block
        
        block_u_0(:) = 0;       % remove blockages, to have channels all in LoS
        block_u_1(:) = 0;
        block_Bu(:) = 0;
        
        %% Channels
        % TODO: check why the steering angle of the RIS 1 and 0 are the same... That shouldn't be the case.
        [a_R_0_los,~,G_0_los,h_0_los,h_D_los,a_BSU_los,a_BSR_0_los] = compute_channels(d_0_G,d_u_0,d_Bu,[beta,beta_nlos],N,M,phiB_0_a,phiB_0_d,phiU_00,phiUB,el_dist, block_u_0, block_Bu);
        [a_R_1_los,~,G_1_los,h_1_los,~          ,~    ,a_BSR_1_los] = compute_channels(d_1_G,d_u_1,d_Bu,[beta,beta_nlos],N,M,phiB_1_a,phiB_1_d,phiU_10,phiUB,el_dist, block_u_1, block_Bu);
        
        
        block_u_0(:) = 1;       % remove blockages, to have channels all in LoS
        block_u_1(:) = 1;
        block_Bu(:) = 1;
        
        [a_R_0_nlos,~,G_0_nlos,h_0_nlos,h_D_nlos,a_BSU_nlos,a_BSR_0_nlos] = compute_channels(d_0_G,d_u_0,d_Bu,[beta,beta_nlos],N,M,phiB_0_a,phiB_0_d,phiU_00,phiUB,el_dist, block_u_0, block_Bu);
        [a_R_1_nlos,~,G_1_nlos,h_1_nlos,~       ,~         ,a_BSR_1_nlos] = compute_channels(d_1_G,d_u_1,d_Bu,[beta,beta_nlos],N,M,phiB_1_a,phiB_1_d,phiU_10,phiUB,el_dist, block_u_1, block_Bu);
        
        a_R_1_nlos(:) = 0 +1i*10^-25;
        G_1_nlos(:) = 0 +1i*10^-25;
        h_1_nlos(:) = 0 +1i*10^-25; 
        a_BSR_1_nlos(:) = 0 +1i*10^-25;
        
        h_D_los = h_D_los.*direct_factor;
        
        %% ORACLE MARISA estimation with LoS channels
        gain_ue_0_los = d_u_0.^-beta;
        gain_ue_1_los = d_u_1.^-beta;
        gain_bs_los = d_0_G.^-beta;
        
        gain_ue_0_nlos = d_u_0.^-beta_nlos;
        gain_ue_1_nlos = d_u_1.^-beta_nlos;
        gain_bs_nlos   = d_0_G.^-beta_nlos;        
        
        [oracle_v1_0] = MARISA_FAIR(N,phiB_0_a,gain_bs_los,phiU_00,gain_ue_0_los,el_dist,'sum');
        [oracle_v2_0] = MARISA_FAIR(N,phiB_0_a,gain_bs_los,phiU_00,gain_ue_0_los,el_dist,'wsum');
        
        [oracle_v1_1] = MARISA_FAIR(N,phiB_1_a,gain_bs_los, phiU_10, gain_ue_1_los,el_dist,'sum');
        [oracle_v2_1] = MARISA_FAIR(N,phiB_1_a,gain_bs_los, phiU_10, gain_ue_1_los,el_dist,'wsum');
        
        oracle_v1 = cat(1, oracle_v1_0, oracle_v1_1);
        oracle_v2 = cat(1, oracle_v2_0, oracle_v2_1);
        
        % ESTIMATION MARISA
        h_in = cat(1,h_0_los, h_1_los);
        G = cat(1,G_0_los, G_1_los);
        
        % Equivalent Channel        
        hd = h_D_los;
        Hb_marisa = zeros(2*N+1,M,K);
        for k=1:K
            Hb_marisa(:,:,k) = [sqrt(eta).*diag(h_in(:,k)')*G(:,:); hd(:,k)';];
        end
        Hb_risma = zeros(2*N+1,M,K);
        for k=1:K
            Hb_risma(:,:,k) = [diag(h_in(:,k)')*G(:,:); hd(:,k)';];
        end

%% LoS
[srate_marisa_oracle1(k_ind,ind_scen,run_ind),SMSE_marisa_oracle1(k_ind,ind_scen,run_ind),MSE_marisa_oracle1(1:K,k_ind,ind_scen,run_ind), ~,W_marisa_oracle1, SINR_marisa_oracle1(1:K,k_ind,ind_scen,run_ind)] = W_OPT_RISMA(M,P,sigma2n,K,G,h_in.*sqrt(eta),hd,Hb_marisa,oracle_v1);
[srate_marisa_oracle2(k_ind,ind_scen,run_ind),SMSE_marisa_oracle2(k_ind,ind_scen,run_ind),MSE_marisa_oracle2(1:K,k_ind,ind_scen,run_ind), ~,W_marisa_oracle2, SINR_marisa_oracle2(1:K,k_ind,ind_scen,run_ind)] = W_OPT_RISMA(M,P,sigma2n,K,G,h_in.*sqrt(eta),hd,Hb_marisa,oracle_v2);
[srate_risma(k_ind,ind_scen,run_ind)  ,SMSE_risma(k_ind,ind_scen,run_ind)  ,MSE_risma(1:K,k_ind,ind_scen,run_ind)  ,Phi_risma, W_risma, SINR_risma(1:K,k_ind,ind_scen,run_ind)] = RISMA(2*N,M,P,sigma2n,K,G,h_in,hd,Hb_risma,I);

%% RIS 1 in NLoS conditions
% ESTIMATION MARISA
h_in = cat(1,h_0_los, h_1_nlos);
G = cat(1,G_0_los, G_1_nlos);

% Equivalent Channel        
hd = h_D_los;
Hb_marisa = zeros(2*N+1,M,K);
for k=1:K
    Hb_marisa(:,:,k) = [sqrt(eta).*diag(h_in(:,k)')*G(:,:); hd(:,k)';];
end
Hb_risma = zeros(2*N+1,M,K);
for k=1:K
    Hb_risma(:,:,k) = [diag(h_in(:,k)')*G(:,:); hd(:,k)';];
end

[srate_marisa_oracle1_nlos(k_ind,ind_scen,run_ind),SMSE_marisa_oracle1_nlos(k_ind,ind_scen,run_ind),MSE_marisa_oracle1_nlos(1:K,k_ind,ind_scen,run_ind), ~,W_marisa_oracle1_nlos, SINR_marisa_oracle1_nlos(1:K,k_ind,ind_scen,run_ind)] = W_OPT_RISMA(M,P,sigma2n,K,G,h_in.*sqrt(eta),hd,Hb_marisa,oracle_v1);
[srate_marisa_oracle2_nlos(k_ind,ind_scen,run_ind),SMSE_marisa_oracle2_nlos(k_ind,ind_scen,run_ind),MSE_marisa_oracle2_nlos(1:K,k_ind,ind_scen,run_ind), ~,W_marisa_oracle2_nlos, SINR_marisa_oracle2_nlos(1:K,k_ind,ind_scen,run_ind)] = W_OPT_RISMA(M,P,sigma2n,K,G,h_in.*sqrt(eta),hd,Hb_marisa,oracle_v2);
[srate_risma_nlos(k_ind,ind_scen,run_ind)  ,SMSE_risma_nlos(k_ind,ind_scen,run_ind)  ,MSE_risma_nlos(1:K,k_ind,ind_scen,run_ind)  ,Phi_risma_nlos, W_risma_nlos, SINR_risma_nlos(1:K,k_ind,ind_scen,run_ind)] = RISMA(2*N,M,P,sigma2n,K,G,h_in,hd,Hb_risma,I);

        if plot_flag

            
            plot_geometry(cat(2,ris_0_el,ris_1_el),ue,bs)
            
            
            my_polarplot_v([oracle_v1_0, oracle_v2_0])
            legend('MARISA sum','MARISA wsum','MARISA maxmin', 'RISMA')
        end
    end
end
end

save('RESULTS_TEST_NLOS')

%% plot results
load('RESULTS_TEST_NLOS')
avg_srate_marisa_oracle1=mean(srate_marisa_oracle1,3);
avg_srate_marisa_oracle2=mean(srate_marisa_oracle2,3);
avg_srate_risma=mean(srate_risma,3);

avg_srate_marisa_oracle1_nlos=mean(srate_marisa_oracle1_nlos,3);
avg_srate_marisa_oracle2_nlos=mean(srate_marisa_oracle2_nlos,3);
avg_srate_risma_nlos=mean(srate_risma_nlos,3);


avg_SMSE_marisa_oracle2=mean(SMSE_marisa_oracle2 ,3);
avg_SMSE_marisa_oracle1=mean(SMSE_marisa_oracle1 ,3);
avg_SMSE_risma=mean(SMSE_risma,3);

avg_SMSE_marisa_oracle2_nlos=mean(SMSE_marisa_oracle2_nlos ,3);
avg_SMSE_marisa_oracle1_nlos=mean(SMSE_marisa_oracle1_nlos ,3);
avg_SMSE_risma_nlos=mean(SMSE_risma_nlos,3);

avg_sinr_marisa_oracle1= 10*log10(squeeze(mean(sum(SINR_marisa_oracle1,1),4))./K_vec');
avg_sinr_marisa_oracle2= 10*log10(squeeze(mean(sum(SINR_marisa_oracle2,1),4))./K_vec');
avg_sinr_risma= 10*log10(squeeze(mean(sum(SINR_risma,1),4))./K_vec');

avg_sinr_marisa_oracle1_nlos= 10*log10(squeeze(mean(sum(SINR_marisa_oracle1_nlos,1),4))./K_vec');
avg_sinr_marisa_oracle2_nlos= 10*log10(squeeze(mean(sum(SINR_marisa_oracle2_nlos,1),4))./K_vec');
avg_sinr_risma_nlos= 10*log10(squeeze(mean(sum(SINR_risma_nlos,1),4))./K_vec');

rate_marisa_oracle2=log2(1+SINR_marisa_oracle2);
rate_marisa_oracle1=log2(1+SINR_marisa_oracle1);
rate_risma=log2(1+SINR_risma);

rate_marisa_oracle2_nlos=log2(1+SINR_marisa_oracle2_nlos);
rate_marisa_oracle1_nlos=log2(1+SINR_marisa_oracle1_nlos);
rate_risma_nlos=log2(1+SINR_risma_nlos);


%% plots

figure('Name','LoS/NLos comparison','NumberTitle','off')

subplot(2,2,1)
hold on
plot(K_vec, avg_srate_marisa_oracle1)
plot(K_vec, avg_srate_marisa_oracle2)
plot(K_vec, avg_srate_risma)
title('sum rate LoS')
legend('MARISA sum','MARISA wsum', 'RISMA')
xlabel('Nr. users')
ylabel('b/s/Hz')

subplot(2,2,2)
hold on
plot(K_vec, avg_sinr_marisa_oracle1)
plot(K_vec, avg_sinr_marisa_oracle2)
plot(K_vec, avg_sinr_risma)
title('SINR LoS')
legend('MARISA sum','MARISA wsum', 'RISMA')
xlabel('Nr. users')
ylabel('dB')

subplot(2,2,3)
hold on
plot(K_vec, avg_srate_marisa_oracle1_nlos)
plot(K_vec, avg_srate_marisa_oracle2_nlos)
plot(K_vec, avg_srate_risma_nlos)
title('sum rate NLoS')
legend('MARISA sum','MARISA wsum', 'RISMA')
xlabel('Nr. users')
ylabel('b/s/Hz')

subplot(2,2,4)
hold on
plot(K_vec, avg_sinr_marisa_oracle1_nlos)
plot(K_vec, avg_sinr_marisa_oracle2_nlos)
plot(K_vec, avg_sinr_risma_nlos)
title('SINR NLoS')
legend('MARISA sum','MARISA wsum', 'RISMA')
xlabel('Nr. users')
ylabel('dB')


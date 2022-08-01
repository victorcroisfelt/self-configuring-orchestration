close all
clear all

addpath('functions')

plot_flag = false();

% M_vec = [32, 64];
% N_vec = [16, 32];
% K_vec = [8, 16];
% L_vec = [16, 32];
% scenario_size_vec = [250, 500];
% 
% for M = M_vec
%     for Nx = N_vec
%         for K = K_vec
%             for L = L_vec
%                 for scenario_size = scenario_size_vec
%                     main_marisa_extension
%                 end
%             end
%         end
%     end
% end

M = 64
Nx = 16
K = 16
L = 128
scenario_size = 250

main_marisa_extension_parfor
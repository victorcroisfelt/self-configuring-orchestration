function [SINR, SE] = channel_estiamation_MMIMO(Theta_prob, Theta_opt, M, C, L, K, tau_est, tau_c, sigma2n,P_ue, G, H, H_D, eta)

% evolution of theta over time during channel estimation phase
Theta_over_blocks = cat(3, Theta_prob, repmat(Theta_opt,[1,1,L-C]));

% equivalent channel during channel estimation
[H_circ] = equivalent_BS_UE_channel(Theta_over_blocks, H, G, H_D, eta);

% equivalent channel after channel estimation
[H_circ_star] = equivalent_BS_UE_channel(Theta_opt, H, G, H_D, eta);

% noise realization during channel estimation
N_b = (randn(M,K,L) + 1i*randn(M,K,L))*sqrt(sigma2n/2);

% desired estimate of the channel
H_bar = 1./L.*sum(H_circ,3);

% channel estimator input
Y_bar = sqrt(P_ue).*K.* H_bar + 1./L.*sum(N_b,3);

H_bar_hat = Y_bar ./ sqrt(P_ue)./K;

% noise realization during communication
%N_b_comm = (randn(M,K) + 1i*randn(M,K))*sqrt(sigma2n/2);

N_b_comm_test = (randn(M,K*1000) + 1i*randn(M,K*1000))*sqrt(sigma2n/2);
%A = mean(abs(H_bar_hat(:,2)'*N_b_comm_test).^2);


S = zeros(K,1);
I = zeros(K,1);
N = zeros(K,1);

for k =1:K
    S(k) = P_ue*abs(H_bar_hat(:,k)'*H_circ_star(:,k)).^2;
    I(k) = P_ue*sum(abs(H_bar_hat(:,1:K ~= k)'*H_circ_star(:,k)).^2);
    N(k) = mean(abs(H_bar_hat(:,k)'*N_b_comm_test).^2);
end

SINR = S./(I+N);
SE = (tau_c-tau_est)./tau_c.*log2(1+SINR);


end
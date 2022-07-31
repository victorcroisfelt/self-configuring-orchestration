function [SINR, SE, SINR_bound, SE_bound] = channel_estiamation_MMIMO(Theta_prob, Theta_opt, M, C, L, K, tau_est, tau_c, sigma2n,P_ue, G, H, H_D, eta)

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

N_b_comm_test = (randn(M,K*100) + 1i*randn(M,K*100))*sqrt(sigma2n/2);
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

V_bound = zeros(K,K);
E_bound = zeros(K,K);
for k =1:K
    for  i = 1:K
        V_bound(k,i) =  sum((abs(H_bar_hat(:,k)).^2+ 1./(L*K*K)*sigma2n/P_ue).*abs(H_circ(:,i)).^2);
        for m =1:M
            for m_prime = 1:M
                if m_prime ~= m
                    E_bound(k,i) = E_bound(k,i) + real(H_bar(m,k).*conj(H_bar(m_prime,k)).* conj(H_circ_star(m,i)).*H_circ(m_prime,k));
                end
            end
        end
    end
end


S_bound = zeros(K,1);
I_bound = zeros(K,1);
N_bound = zeros(K,1);
for k =1:K
    S_bound(k) = P_ue*abs(H_bar_hat(:,k)'*H_circ_star(:,k)).^2;
    I_bound(k) = P_ue*sum(V_bound(k,:)+E_bound(k,:))-S_bound(k);
    N_bound(k) = sigma2n*(abs(H_bar(:,k)'*H_bar(:,k)) + M/(L*K*K).*(sigma2n/P_ue));
end

SINR_bound = S_bound./(I_bound+N_bound);
SE_bound = (tau_c-tau_est)./tau_c.*log2(1+SINR_bound);


end
function [MSE] = channel_estimation_MSE(M,C,L,K,sigma2n,power, G, Theta_prob, Theta_opt, H, H_D, eta)

% evolution of theta over time during channel estimation phase
Theta_over_blocks = cat(3, Theta_prob, repmat(Theta_opt,[1,1,L-C]));

% equivalent channel during channel estimation
[H_circ] = equivalent_BS_UE_channel(Theta_over_blocks, H, G, H_D, eta);

% equivalent channel after channel estimation
[H_circ_star] = equivalent_BS_UE_channel(Theta_opt, H, G, H_D, eta);

% desired estimate of the channel
H_bar = 1./L * sum(H_circ, 3);

%intermediate_matrix = G'*(sum(Theta_prob, 3)- C*Theta_opt)*H;

intermediate_vector = zeros(K,1);

for k = 1:K
    intermediate_vector(k) = norm(H_bar(:, k) - H_circ_star(:, k))^2;  %the abs is to have real values in case of numerical errors
end

MSE = (M/(L*K))*(sigma2n/power) + mean(intermediate_vector);
end
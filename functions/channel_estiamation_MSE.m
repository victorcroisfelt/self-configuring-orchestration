function [MSE] = channel_estiamation_MSE(M,L,K,sigma2n,power, G, Theta_prob, Theta_opt, H, eta)

C = size(Theta_prob, 3);

intermediate_matrix = G'*(sum(Theta_prob, 3)- C*Theta_opt)*H;

intermediate_vector = zeros(K,1);
for k = 1:K
    intermediate_vector(k) = abs(intermediate_matrix(:,k)'*intermediate_matrix(:,k));  %the abs is to have real values in case of numerical errors
end

MSE =( M/(L*K^2))*(sigma2n/power)+ (eta/L^2)*(1/K)*sum(intermediate_vector);
end
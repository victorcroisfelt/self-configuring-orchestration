function [H_circ] = equivalent_BS_UE_channel(Theta, H, G, H_D, eta)
L = size(Theta,3);  % number of subblocks
M = size(G,2);      % number of BS elements
K = size(H_D,2);    % number of users

H_circ = zeros(M,K,L);

%H_circ_norm = zeros(K,L);

for l = 1:L
    
    H_circ(:,:,l) = H_D + sqrt(eta).*G'*Theta(:,:,l)'*H;
    %for k =1:K
    %    H_circ_norm(k,l) = H_circ(:,k,l)'*H_circ(:,k,l);
    %end
end
end
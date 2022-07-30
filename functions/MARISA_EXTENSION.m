function [Theta_out,detected_ue] = MARISA_EXTENSION(Y, Codebook, threshold, phi_B, sigma2n, mode)

N = size(Y,1);
K = size(Y,2);
C = size(Y,3);
D = size(Y,3);

% steering_vector to BS
g = steer(N,phi_B,0.5);
g_sigma = sum(g,2);
v_B = g_sigma./abs(g_sigma);
%v_B = diag(v_B);





if strcmpi(mode,'signal')
    Comb_matrix = zeros(N,D);
    
    if numel(size(Codebook)) == 3
        for d =1:D
            Comb_matrix(:,d) = diag(Codebook(:,:,d));
        end
    else
        Comb_matrix = Codebook;
    end
    
    Comb_matrix = (1/N).*Comb_matrix;
    
    %Y_comb = zeros(D,K,C);
    %for c = 1:C
    %    Y_comb(:,:,c) = Comb_matrix'*Y(:,:,c);
    %end
    Y_comb = reshape(Comb_matrix'*reshape(Y,[N,K*C]),[D,K,C]);
    
    Y_bar = (1/C).*sum(Y_comb,3);
    Y_bar_pow = abs(Y_bar).^2;
    
    max_pow = max(Y_bar_pow, [], 1);
    
    detected_ue = max_pow > (2*N*C/sigma2n)^-1*threshold;
    
    Theta_hat = Comb_matrix(:,detected_ue);
    A_hat = max_pow(detected_ue);
    
end

if strcmpi(mode,'power')
    theta=zeros(N,D);
    
    if numel(size(Codebook))==3
        for d =1:D
            theta(:,d) = diag(Codebook(:,:,d));
        end
    else
        theta = Codebook;
    end
    
    Y_d = zeros(D,K);
    for d = 1:D
        Y_d(d,:) = theta(:,d)'*Y(:,:,d);
    end
    
    
    Y_bar = Y_d;
    Y_bar_pow = abs(Y_bar).^2;
    
    max_pow = max(Y_bar_pow, [], 1);
    
    detected_ue = max_pow > threshold;
    
    Theta_hat = theta(:,detected_ue);
    A_hat = max_pow(detected_ue);
end


v_U = sum(Theta_hat.*sqrt(A_hat)./max(sqrt(A_hat)),2);

v_U = v_U./max(abs(v_U));
theta_out = v_U.*conj(v_B);
Theta_out = diag(theta_out);


end
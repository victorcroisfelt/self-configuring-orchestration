function [srate_ao,SMSE_ao,MSE_ao,Phim,W,SINR_ao] = W_OPT_RISMA(M,P,sigma2n,K,G,h,hd,Hb,Phim)

% Implements the Precoder optimized according to RISMA algorithm: created 19 July 2021
% Usage: [srate,SMSE,MSE] = RISMA(Nx,Ny,M,P,sigma2n,K,Rad,fig_map) where
% Nx = number of antenna elements of the RIS along x-axis
% Ny = number of antenna elements of the RIS along y-axis
% M = number of antenna elements of the BS
% P = Available power budget at the BS [W]
% sigma2n = Noise power at the UEs [W]
% K = Number of UEs
% Gammas = square root covariance matrix of the BS-RIS link
% Rs = square root covariance matrix of the RIS-UEs link
% Ts = square root covariance matrix of the BS-UEs link
% I = Number of Monte Carlo simulations

% Output:
% srate = Sum rate for all the I montecarlo runs
% SMSE = SMSE for all the I montecarlo runs
% MSE = MSE of each UE for all the I montecarlo runs


Phim = diag(Phim);


warning('off','all');

    v= diag(Phim');
    hb = zeros(size(Hb,2),size(Hb,3));
    for k=1:size(Hb,3)
        hb(:,k) = Hb(:,:,k)' * [v;1] ;
    end
    
    mu = K*sigma2n/P;
    W = ((hb*hb')+ mu*eye(M))\hb;
    W = W * sqrt(P)/norm(W,'fro');
    
    MSE_ao = sum(abs((h'*Phim*G + hd')*W).^2) - 2*(real(diag((h'*Phim*G + hd')*W)))' + 1+ sigma2n;
    MSE_ao = MSE_ao';
    SMSE_ao = sum(MSE_ao); % sum(sum(abs((h'*Phim*G + hd')*W).^2)) - 2*sum(real(diag((h'*Phim*G + hd')*W))) +  K*(1+ sigma2n);

    SINR_ao = zeros(1,K);
    for k=1:K
        jk = ones(1,K);
        jk(k) = 0;
        SINR_ao(k) = abs((h(:,k)'*Phim*G + hd(:,k)')*W(:,k))^2./...
            (sigma2n + sum(abs((h(:,k)'*Phim(:,:)*G + hd(:,k)')*W(:,logical(jk))).^2));
    end
    
    srate_ao = sum(log2(1+SINR_ao));


end
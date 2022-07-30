function [srate_ao,SMSE_ao,MSE_ao,Phi_min,W_min,SINR_ao] = RISMA(N,M,P,sigma2n,K,G_i,h_i,hd_i,Hb_i,I)

% Implements the RISMA algorithm: created 19 July 2021
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

%rng(123);


Maxiter = 102;%22;
Tol = 1e-5;

SMSE_ao = zeros(I,1); % With Alt. Opt.
MSE_ao = zeros(K,I);
srate_ao = zeros(I,1);

warning('off','all');

Phi_min = zeros(N,N,I);
W_min = zeros(M,K,I);

for i=1:I
    G = G_i(:,:,i);
    h = h_i(:,:,i);
    hd = hd_i(:,:,i);
    Hb = Hb_i(:,:,:,i);
    stop = 0;
    ll = 2;
    SMSE_ao_i = zeros(Maxiter,1);
    MSE_ao_i = zeros(Maxiter,K);
    SMSE_ao_i(1) = 200;
    SMSE_ao_min_i = 200;
    Phim = diag(exp(1i*rand(N,1)*pi).*ones(N,1));
    v= diag(Phim');
    hb = zeros(size(Hb,2),size(Hb,3));
    for k=1:size(Hb,3)
        hb(:,k) = Hb(:,:,k)' * [v;1] ;
    end
    
    mu = K*sigma2n/P;
    W = ((hb*hb')+ mu*eye(M))\hb;
    W = W * sqrt(P)/norm(W,'fro');
    
    In = eye(N+1);
    Hbb_k = zeros(size(Hb,1),size(Hb,1),K);
    hbb = zeros(size(Hb,1),K);
    while ~stop
        %                         v = cvx_min(Hb,W);
        WW = (W * W');
        for k=1:K
            Hbb_k(:,:,k) = Hb(:,:,k) * WW * Hb(:,:,k)';
            hbb(:,k) = Hb(:,:,k) * W(:,k);
        end
        Hbb = sum(Hbb_k,3);
        mu = sigma2n*ones(N+1,1);
        B = (Hbb + diag(mu));
        z = sum(hbb,2);
        nu = (In(:,N+1)'*(B\z)-1)/(B(N+1,N+1));
        vv = (B)\(z - nu*In(:,N+1));
        %v = vv(1:N)/norm(vv(1:N));
        %v = v/max(abs(v));
        v = vv(1:N)./abs(vv(1:N));
        
        hb = zeros(size(Hb,2),size(Hb,3));
        for k=1:size(Hb,3)
            hb(:,k) = Hb(:,:,k)' * [v;1] ;
        end
        
        mu = K*sigma2n/P;
        W = ((hb*hb')+ mu*eye(M))\hb;
        W = W * sqrt(P)/norm(W,'fro');
        
        
        Phim = diag(v');
        MSE_ao_i(ll,:) = sum(abs((h'*Phim*G + hd')*W).^2) - 2*(real(diag((h'*Phim*G + hd')*W)))' + 1+ sigma2n;
        SMSE_ao_i(ll) = sum(MSE_ao_i(ll,:)); % sum(sum(abs((h'*Phim*G + hd')*W).^2)) - 2*sum(real(diag((h'*Phim*G + hd')*W))) +  K*(1+ sigma2n);
        if SMSE_ao_i(ll)<SMSE_ao_min_i
            Phi_min(:,:,i) = Phim;
            W_min(:,:,i) = W;
            SMSE_ao_min_i = SMSE_ao_i(ll);
            MSE_ao_min_i = MSE_ao_i(ll,:);
        end
        
        if ll>Maxiter-2 || (abs(SMSE_ao_i(ll,:) - SMSE_ao_i(ll-1,:))./abs(SMSE_ao_i(ll-1,:)))<Tol %|| (SMSE_ao_i(ll) - SMSE_ao_i(ll-1))>0
            stop =1;
        end
        
        ll = ll+1;
    end
    SMSE_ao(i) = min(SMSE_ao_i(1:ll-1));
    MSE_ao(:,i) = MSE_ao_min_i;
    SINR_ao = zeros(1,K);
    for k=1:K
        jk = ones(1,K);
        jk(k) = 0;
        SINR_ao(k) = abs((h(:,k)'*Phi_min(:,:,i)*G + hd(:,k)')*W_min(:,k,i))^2./(sigma2n + sum(abs((h(:,k)'*Phi_min(:,:,i)*G + hd(:,k)')*W_min(:,logical(jk),i)).^2));
    end
    
    srate_ao(i) = sum(log2(1+SINR_ao));
end


end
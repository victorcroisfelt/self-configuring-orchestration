function [v, v_U, v_B] = MARISA_FAIR(N,phi_B_a,gain_G,phi_U,gain_U,d_lambda,mode)

h = steer(N,phi_U,d_lambda);
g = steer(N,phi_B_a,d_lambda);

if strcmpi(mode,'sum')
    h_sigma = sum(h,2);
    v_U = h_sigma./abs(h_sigma);
    
    g_sigma = sum(g,2);
    v_B = g_sigma./abs(g_sigma);
end
if strcmpi(mode,'wsum')
    h_sigma = sum(h.*sqrt(gain_U)',2);
    %v_U = h_sigma./norm(h_sigma)*sqrt(N);
    v_U = h_sigma./abs(h_sigma);
    
    g_sigma = sum(g.*sqrt(gain_G)',2);
    v_B = g_sigma./abs(g_sigma);
end
if strcmpi(mode,'maxmin')
    H=zeros(N,N,length(phi_U));
    for u =1:length(phi_U)
        H(:,:,u)=h(:,u)*h(:,u)';
    end
    v_U = run_max_min_solver(H,gain_U);
    
    
    G=zeros(N,N,length(phi_G));
    for u =1:length(phi_G)
        G(:,:,u)=g(:,u)*g(:,u)';
    end
    v_B = run_max_min_solver(G,gain_G);    
end

v = v_U.*conj(v_B);
end
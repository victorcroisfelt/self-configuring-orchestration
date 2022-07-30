function [C_out,C_in, phi_vec] = generate_codebook(M,N)
% Design codebook by spanning the entire straight angle [0, pi]
%OUTPUTS:
%     C 
%     phi_vec
% INPUTS:
%     M = codebook size 
%     N = number of RIS elements

phi_vec = linspace(0,pi,M);           % angular span

C_out = zeros(N,M);                             % Codebook
C_in = zeros(N,M);                             % Codebook

for i = 1:M
    a_in = steer(N, phi_vec(i),0.5);
    C_out(:,i) = conj(a_in);%./norm(a_in);
    C_in(:,i) = a_in;%./norm(a_in);
end




end

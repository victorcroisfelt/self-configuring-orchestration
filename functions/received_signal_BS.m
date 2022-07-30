function [Y] = received_signal_BS(H_circ, pilot, rho, sigma2n)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
K = size(pilot,1);  % nr users
M = size(H_circ,1);
L = size(H_circ,3);  % nr ris el

Y = zeros(M,K,L);

for l = 1:L
    Y(:,:,l) = H_circ(:,:,l)*pilot;
end
Y = Y * sqrt(rho);
Y = Y + (randn(M,K,L) + 1i*randn(M,K,L))*sqrt(sigma2n/2);
end


function [Y] = received_signal_RIS(Theta, H, eta, pilot, rho,sigma2n)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
K = size(pilot,1);  % nr users
N = size(Theta,1);  % nr ris el
C = size(Theta,3);  % probings

Y = zeros(N,K,C);

for c = 1:C
    Y(:,:,c) = Theta(:,:,c)*H;
end
Y = Y * sqrt(1-eta) * sqrt(rho);
Y = Y + (randn(N,K,C) + 1i*randn(N,K,C))*sqrt(sigma2n/2);
end


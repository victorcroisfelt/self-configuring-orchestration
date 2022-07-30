function [v] = steer(M,theta,d_lambda)
v = zeros(M,length(theta));
for j = 1:length(theta)
    v(:,j)=exp(1i*2*pi*d_lambda*kron((0:M-1)',cos(theta(j))));
end
end
function [a_R,g,G,h,h_D,a_BSU,a_BSR]= compute_channels(d_G,d_u,d_Bu,beta,N,M,phiB_a,phiB_d,phiU,phiUB,el_dist, varargin)
%a_R = ris response at direction of BS
%g = bs response at direction of RIS
%G = Overall BS_RIS channel
%h = RIS_UEs channel
%h_D = BS UE channel


% Nlos case
if numel(varargin) > 0
    nlos = true();
    block_u = varargin{1};
    block_Bu = varargin{2};
else
    nlos = false();
    block_u = false(size(d_u));
    block_Bu = false(size(d_Bu));
end

gamma_G = d_G.^-beta(1);                % Average pathloss of BS-RIS link
gamma_h = d_u.^-beta(1);                % Average pathloss of RIS-UE link
gamma_h_D = d_Bu.^-beta(1);             % Average pathloss of BS-UE link

if nlos
    gamma_h(block_u) = d_u(block_u).^-beta(2);                % Average pathloss of RIS-UE link
    gamma_h_D(block_Bu) = d_Bu(block_Bu).^-beta(2);             % Average pathloss of BS-UE link
end



G = sqrt(gamma_G)* steer(N,phiB_a,el_dist)*steer(M,phiB_d,el_dist)';            % BS-RIS channel
a_R = steer(N,phiB_a,el_dist);                                                  % BS-RIS channel (RIS steering)
g = sqrt(gamma_G)'.*steer(M,phiB_d,el_dist);                                    % BS-RIS channel (BS steering)
h = sqrt(gamma_h)'.* steer(N,phiU,el_dist);                                     % RIS-UE channel
h_D = sqrt(gamma_h_D)'.* steer(M,phiUB,el_dist);                                % BS-UE channel

a_BSU = steer(M,phiUB,el_dist);
a_BSR = steer(M,phiB_d,el_dist);

end
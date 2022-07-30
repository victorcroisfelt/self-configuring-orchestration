function [block] = blockage_path(d, lambda_b, r_B, h_B, h_U, h_A)
% blockage probability as per  Koucheryavy, Yevgeni, et al.
%"Quantifying the millimeter wave new radio base stations density
%for network slicing with prescribed SLAs." Computer Communications 174 (2021): 13-27.
% formula (5)

p = 1-exp(-2.*lambda_b.*r_B.*(sqrt(d).*(h_B-h_U)./(h_A-h_U)+r_B));

block = rand(size(d))<= p;


end
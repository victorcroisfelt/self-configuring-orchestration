function [outputArg1,outputArg2] = th_detection_performance( ...
   
    inputArg1, ...
    inputArg2)

%TH_DETECTION_PERFORMANCE Summary of this function goes here
%   Detailed explanation goes here

if strcmpi(mode,'signal')

    threshold = chi2inv(prob_false_alarm, 2);

end

if strcmpi(mode,'power')

    threshold = 2 * num_elements * noise_power * log(2 * num_elements * noise_power * prob_false_alarm);

end

end


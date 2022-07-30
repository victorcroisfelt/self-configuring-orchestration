function [threshold] = detection_threshold( ...
        mode, ...
        prob_false_alarm, ...
        num_elements, ...
        noise_power ...
        )

if strcmpi(mode,'signal')

    threshold = chi2inv(prob_false_alarm, 2);

end

if strcmpi(mode,'power')

    threshold = 2 * num_elements * noise_power * log(2 * num_elements * noise_power * prob_false_alarm);

end
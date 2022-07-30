function [threshold] = detection_threshold( ...
        mode, ...
        prob_false_alarm, ...
        num_elements, ...
        noise_power ...
        )

if strcmpi(mode,'signal')

    if prob_false_alarm == 0.1

        threshold = 0.2107;

    end

    if prob_false_alarm == 0.01

        threshold = 0.0201;

    end

    if prob_false_alarm == 0.001

        threshold = 0.0020;

    end
       
    %threshold = chi2inv(prob_false_alarm, 2);

end

if strcmpi(mode,'power')

    threshold = - (2 / (num_elements * noise_power)) * log(prob_false_alarm);

end
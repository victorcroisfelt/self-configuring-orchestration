function [threshold] = detection_threshold( ...
        mode, ...
        prob_false_alarm, ...
        num_elements, ...
        noise_power ...
        )

if strcmpi(mode,'signal')

    if prob_false_alarm == 0.1

        threshold = 4.6052;

    end

    if prob_false_alarm == 0.01

        threshold = 9.2103;

    end

    if prob_false_alarm == 0.001

        threshold = 13.8155;

    end
       
    %threshold = chi2inv(1 - prob_false_alarm, 2);

end

if strcmpi(mode,'power')

    threshold = - (2 * num_elements * noise_power) * log(prob_false_alarm);

end
function [steering]= array_steering_vector(lambda, pos1, pos2, pos_el)

if length(pos1) > 1


    steering = zeros()


    % Compute wave vector
    wave_vector = 2 * pi / lambda * (pos1 - pos2) / norm(pos2 - pos1)^2;
    
    % Compute inner product
    angles = dot(wave_vector, (pos_el - pos2));
    
    % Compute steering vector
    steering = exp(1j * angles);

    end

end
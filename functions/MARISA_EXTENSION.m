function [Theta_out, hat_prob_detection, true_prob_detection, distance] = MARISA_EXTENSION(Y, Codebook, prob_false_alarm, phi_B, sigma2n, mode)

% Extract sizes 
N = size(Y,1); % number of RIS elements
K = size(Y,2); % number of UEs
C = size(Y,3); % number of configurations
D = size(Y,3); % number of directions

% Compute steering vector to BS
g = steer(N, phi_B, 0.5);
g_sigma = sum(g, 2);

% Direction towards BS
v_B = g_sigma./abs(g_sigma);
%v_B = diag(v_B);

% Get threshold given a probability of false alarm
threshold = detection_threshold(mode, prob_false_alarm, N, sigma2n);

if strcmpi(mode,'signal')
    
    % Prepare to save combining matrix
    Comb_matrix = zeros(N, D);
    
    % Obtain combining matrix
    if numel(size(Codebook)) == 3
        for d =1:D
            Comb_matrix(:, d) = diag(Codebook(:, :, d));
        end
    else
        Comb_matrix = Codebook;
    end
    
    % Normalizing
    Comb_matrix = (1/N) .* Comb_matrix;
    
    Y_comb = zeros(D,K,C);
    for c = 1:C
        Y_comb(:,:,c) = Comb_matrix'*Y(:,:,c);
    end

    % Compute received signal according to Eq. (9)
    %Y_comb = reshape(Comb_matrix'*reshape(Y,[N, K*C]), [D, K, C]);
    
    % Compute average received signal according to Eq. (13)
    Y_bar = (1/C).*sum(Y_comb, 3);

    % Compute alphas as defined below Eq. (14) 
    Y_bar_pow = abs(Y_bar).^2;
    
    % Obtain the power received at the best direction according to Eq. (15)
    [max_pow, max_ind] = max(Y_bar_pow, [], 1);

    % Apply detection test in Eq. (22)
    detected_ue = max_pow > (2*N*C/sigma2n)^-1 * threshold;
    
    % Obtain estimated CSI for the detected UEs in Eqs. (16) and (17) 
    Theta_hat = Codebook(:, max_ind);
    Theta_hat_masked = Theta_hat(:, detected_ue);
    A_hat = max_pow(detected_ue);
    
    % Compute true probability of detection
    true_prob_detection = zeros(K, 1);

    for k = 1:K

        true_prob_detection(k) = Qchipr2(2, (2*N*C/sigma2n) * max_pow(k), threshold, 1e-5);
    
    end

    true_prob_detection(isinf(true_prob_detection)|isnan(true_prob_detection)) = 1;
    true_prob_detection = mean(true_prob_detection);

end

if strcmpi(mode,'power')

    % Prepare to save codebook
    theta = zeros(N, D);
    
    % Obtain codebook as in Eq. (18)
    if numel(size(Codebook))==3
        for d =1:D
            theta(:, d) = diag(Codebook(:, :, d));
        end
    else
        theta = Codebook;
    end
    
    % Compute received signal at each direction according to Eq. (19)
    Y_d = zeros(D,K);
    for d = 1:D
        Y_d(d,:) = theta(:, d)'*Y(:, :, d);
    end
    
    % Compute received power according to Eq. (20)
    Y_bar = Y_d;
    Y_bar_pow = abs(Y_bar).^2;
    
    % Obtain the power received at the best direction
    [max_pow, max_ind] = max(Y_bar_pow, [], 1);
    
    % Apply detection test in Eq. (30)
    detected_ue = max_pow > threshold;
    
    % Obtain estimated CSI for the detected UEs
    Theta_hat = theta(:, max_ind);
    Theta_hat_masked = Theta_hat(:, detected_ue);
    A_hat = max_pow(detected_ue);
    
    % Compute true probability of detection 
    true_prob_detection = exp(-(2 * N * sigma2n)^(-1) * (threshold - max_pow));

    pos = (threshold - max_pow < 0);
    true_prob_detection(pos) = nan;

    %true_prob_detection(isinf(true_prob_detection)|isnan(true_prob_detection)) = 1;

    true_prob_detection = mean(true_prob_detection(pos));

end

% Compute estimated probability of detection
hat_prob_detection = sum(detected_ue) / K;

% Compute desired case 
%v_d = sum(Theta_hat.*max_pow./sum(max_pow), 2);
v_d = mean(Theta_hat, 2);
v_d = sqrt(N) * (v_d/norm(v_d));
v_d(isnan(v_d)) = 0;

% Compute normalized directions according to Eqs. (32) and (33)
%v_U = sum(Theta_hat_masked.*A_hat./sum(A_hat), 2);
v_U = mean(Theta_hat_masked, 2);
v_U = sqrt(N) * (v_U/norm(v_U));
v_U(isnan(v_U)) = 0;

% Compute desired optimal reflection direction
theta_desired = v_d .* conj(v_B);

% Compute optimal reflection direction according to Eq. (32)
theta_out = v_U .* conj(v_B);
Theta_out = diag(theta_out);

% Compute distances between the two
distance = norm(theta_out - theta_desired);

end
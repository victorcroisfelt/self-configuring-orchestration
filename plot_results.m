clear all
close all

load('RESULTS_M64_N16_K16_L128_50x50.mat')

avg_hat_detected_ue_pow = mean(mean(hat_detected_ue_pow,4),3); 
avg_hat_detected_ue_sig = mean(mean(hat_detected_ue_sig,4),3);

avg_th_detected_ue_pow = mean(mean(th_detected_ue_pow,4),3); 
avg_th_detected_ue_sig = mean(mean(th_detected_ue_sig,4),3);

MSE_cha_est_pow = mean(mean(MSE_cha_est_pow,4),3);
MSE_cha_est_sig = mean(mean(MSE_cha_est_sig,4),3);

figure
hold on
p11 = plot(C_vec./L,MSE_cha_est_pow(:,1));
p12 = plot(C_vec./L,MSE_cha_est_sig(:,1));
% p21 = plot(C_vec./L,MSE_cha_est_pow(:,2));
% p22 = plot(C_vec./L,MSE_cha_est_sig(:,2));
% p31 = plot(C_vec./L,MSE_cha_est_pow(:,3));
% p32 = plot(C_vec./L,MSE_cha_est_sig(:,3));
set(gca, 'YScale', 'log')


ylabel('Channel estiation MSE', 'Interpreter', 'latex')
xlabel('$C/L$', 'Interpreter', 'latex')

legend('power-based', 'signal-based', 'Interpreter', 'latex')


figure
hold on
p11 = plot(C_vec./L,avg_hat_detected_ue_pow(:,1));
p12 = plot(C_vec./L,avg_hat_detected_ue_sig(:,1));
p21 = plot(C_vec./L,avg_hat_detected_ue_pow(:,2));
p22 = plot(C_vec./L,avg_hat_detected_ue_sig(:,2));
p31 = plot(C_vec./L,avg_hat_detected_ue_pow(:,3));
p32 = plot(C_vec./L,avg_hat_detected_ue_sig(:,3));


ylabel('Detecion probability', 'Interpreter', 'latex')
xlabel('$C/L$', 'Interpreter', 'latex')

legend(['power-based ', num2str(false_alarm_prob_vec(1))], ['signal-based'], 'Interpreter', 'latex')


function [theta] = channelEstimation(csi_sample)
% theta is a L*3 vector: [alpha, delta, tau] 
% [alpha, -, -]_l: the amplitude of the path l;
% [-, delta, -]_l: the azimuth angle of the path l;
% [-, -, tau]_l: the delay of the path l;

% labmda: wavelength of the signal;
% freq: frequency of the signal;
% M: antenna array size;
freq = 2.4e9;
lambda = 3e8/2.4e9;
M = 3;

% L is the number of propagation paths;
L = 8;

d = 12;
% d is spacing between adjacent antennas of the receiver antenna array;

iteration = 999;
parameter = cell(iteration+1);

theta = parameter_initialization(L);
parameter{1} = theta;


% Iterating 1000 times;
for a = 1:iteration
    
    % Expectation step (E-step)
    X = eStep(parameter{a},freq, L, M, csi_sample, lambda,d);
    
    % Maximization step (M-step)
    theta_next = mStep (X, freq, L, M, lambda,d);
    parameter{a+1} = theta_next;
end



end
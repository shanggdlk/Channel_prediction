function [theta] = channelEstimation(csi_sample)
% csi_sample: M*1 complex
% theta is a L*3 vector: [alpha, phi, tau] 
% [alpha, -, -]_l: the amplitude of the path l;
% [-, phi, -]_l: the azimuth angle of the path l;
% [-, -, tau]_l: the delay of the path l;
if nargin <= 0
    csi_sample = [0.5+0.1j; 0.4+0.2j; 0.3+0.3j];
end


globals_init();

global ITERATION L

% parameters: ITERATION+1 cell, inside each cell is a struct 
% (alpha, phi, tau)
parameters = cell(ITERATION+1, 1);
for i = 1:ITERATION+1
    parameters{i} = struct('alpha', zeros(1, L), ...
    'phi', zeros(1, L), 'tau', zeros(1, L));
end

parameters{1} = init(csi_sample, parameters{1});


% Iterating 1000 times;
for i = 2:ITERATION+1
    
    % Expectation step (E-step)
    X = eStep(parameters{i});
    
    % Maximization step (M-step)
    theta_next = mStep (X, freq, L, M, lambda,d);
    parameter{a+1} = theta_next;
end



end

function globals_init
    %% physical
    % labmda: wavelength of the signal;
    % freq: frequency of the signal;
    % M: antenna array size;
    % L is the number of propagation paths;
    % D is spacing between adjacent antennas of the receiver antenna array;
    % N: # of sample
    global FREQUENCY SPEED_OF_LIGHT LAMBDA M L D ITERATION N
    FREQUENCY = 5.2e9;  %unit hz
    SPEED_OF_LIGHT = 3e8;  %unit m/s
    LAMBDA = SPEED_OF_LIGHT/FREQUENCY;
    M = 3;
    L = 8;
    D = 12;
    ITERATION = 999;
    N = 1;
end
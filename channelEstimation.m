function [ret_parameter] = channelEstimation(csi_sample)
%% csi_sample: M*1 complex
% theta is a L*3 vector: [alpha, phi, tau] 
% [alpha, -, -]_l: the amplitude of the path l;
% [-, phi, -]_l: the azimuth angle of the path l;
% [-, -, tau]_l: the delay of the path l;

globals_init();

%  no input csi -> simulation
if nargin <= 0
    csi_sample = generate_simulation();
end


global ITERATION L

% parameters: ITERATION+1 cell, inside each cell is a struct 
% (alpha, phi, tau)
parameters = cell(ITERATION+1, 1);
for i = 1:ITERATION+1
    parameters{i} = struct('alpha', zeros(1, L), ...
    'phi', zeros(1, L), 'tau', zeros(1, L));
end

parameters{1} = init(csi_sample, parameters{1});


%% Iterating;
for i = 1:ITERATION
    
    % Expectation step (E-step)
    X = eStep(parameters{i}, csi_sample);
    
    % Maximization step (M-step)
    parameters{i+1}.tau = opt_tau(parameters{i}.phi, X);
    parameters{i+1}.phi = opt_phi(parameters{i+1}.tau, X);
    parameters{i+1}.alpha = compute_alpha(parameters{i+1}.tau,... 
        parameters{i+1}.phi, X);
end

ret_parameter = parameters{ITERATION+1};
end

function ret_tau = opt_tau(phi, X)

global DOMAIN_TAU L

tau_space = zeros(DOMAIN_TAU.length, L);
for i = 1:DOMAIN_TAU.length
    tau_space(i,:) = compute_Z(zeros(1, L)+DOMAIN_TAU.step*(i-1)...
        +DOMAIN_TAU.start, phi, X);
end
[~, I] = max(abs(tau_space));
ret_tau = DOMAIN_TAU.step*(I-1)+DOMAIN_TAU.start;

end


function ret_phi = opt_phi(tau, X)

global DOMAIN_PHI L

phi_space = zeros(DOMAIN_PHI.length, L);
for i = 1:DOMAIN_PHI.length
    phi_space(i,:) = compute_Z(tau, zeros(1, L)+DOMAIN_PHI.step*(i-1)...
        +DOMAIN_PHI.start, X);
end
[~, I] = max(abs(phi_space));
ret_phi = DOMAIN_PHI.step*(I-1)+DOMAIN_PHI.start;

end

function ret_alpha = compute_alpha(tau, phi, X)

global M N 

ret_alpha = 1/(M*N)*abs(compute_Z(tau, phi, X));

end


function csi_sample = generate_simulation()
    global SIMULATION_TAU SIMULATION_PHI SPEED_OF_LIGHT M FREQUENCY
    C = compute_C(SIMULATION_PHI);
    ALPHA = 1./(SPEED_OF_LIGHT*SIMULATION_TAU);
    csi_sample = repmat(ALPHA, M, 1) .* C .* ...
        repmat(exp(-1j*2*pi*SIMULATION_TAU*FREQUENCY), M, 1);
    csi_sample = sum(csi_sample, 2);
    
    


end




function globals_init
    %% physical
    % labmda: wavelength of the signal;
    % freq: frequency of the signal;
    % M: antenna array size;
    % L is the number of propagation paths;
    % D is spacing between adjacent antennas of the receiver antenna array;
    % N: # of sample
    global FREQUENCY SPEED_OF_LIGHT LAMBDA M L D ITERATION N DOMAIN_TAU ...
        DOMAIN_PHI
    FREQUENCY = 5.2e9;  %unit hz
    SPEED_OF_LIGHT = 3e8;  %unit m/s
    LAMBDA = SPEED_OF_LIGHT/FREQUENCY;
    M = 3;
    L = 8;
    D = 12;
    ITERATION = 999;
    N = 1;
    DOMAIN_TAU = struct('start', 10e-9, 'end', 30e-9, 'step', 1e-9); % unit: s
    DOMAIN_TAU.length = round((DOMAIN_TAU.end - DOMAIN_TAU.start) ...
        / DOMAIN_TAU.step + 1);
    
    DOMAIN_PHI = struct('start', 0, 'end', pi, 'step', pi/10); % unit: radius
    DOMAIN_PHI.length = round((DOMAIN_PHI.end - DOMAIN_PHI.start) ...
        / DOMAIN_PHI.step + 1);
    
    %% simulation
    global SIMULATION_TAU SIMULATION_PHI
    SIMULATION_TAU = [12 13 15 17 19 20 22 24]*1e-9;
    SIMULATION_PHI = [0.1 0.3 0.4 0.45 0.5 0.6 0.76 0.8]*pi;
   
end
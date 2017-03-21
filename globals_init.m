function globals_init
%% physical
% LAMBDA: wavelength of the signal;
% FREQUENCY: frequency of the signal;
% M: antenna array size;
% L: # propagation paths;
% D: spacing between adjacent rx antenna
% N: # sample
% F: # measured subcarrier
% DELTA_FREQUENCY: difference between adjacent subcarrier
    global CENTRAL_FREQUENCY SPEED_OF_LIGHT M L D ITERATION DOMAIN_TAU ...
        DOMAIN_PHI F DELTA_FREQUENCY FREQUENCIES LAMBDAS
    SPEED_OF_LIGHT = 3e8;  %unit m/s
    CENTRAL_FREQUENCY = 5.2e9;  %unit hz
    DELTA_FREQUENCY = 20e6/64; % 20Mhz, 64 slots
    F = 56;
    FREQUENCIES = ((1:F) - (F+1)/2) * DELTA_FREQUENCY+CENTRAL_FREQUENCY;
    LAMBDAS = SPEED_OF_LIGHT./FREQUENCIES;
    M = 10;
    L = 8;
    D = mean(LAMBDAS)/2;
    ITERATION = 800;
    DOMAIN_TAU = struct('start', 10e-9, 'end', 30e-9, 'step', 1e-9); % unit: s
    DOMAIN_TAU.length = round((DOMAIN_TAU.end - DOMAIN_TAU.start) ...
        / DOMAIN_TAU.step + 1);
    
    DOMAIN_PHI = struct('start', 0, 'end', pi, 'step', pi/100); % unit: radius
    DOMAIN_PHI.length = round((DOMAIN_PHI.end - DOMAIN_PHI.start) ...
        / DOMAIN_PHI.step + 1);
    
    %% simulation
    global SIMULATION_TAU SIMULATION_PHI
    SIMULATION_TAU = [12 13 15 17 19 20 22 24]*1e-9;
    SIMULATION_PHI = [0.1 0.3 0.4 0.45 0.5 0.6 0.76 0.8]*pi;
   
end
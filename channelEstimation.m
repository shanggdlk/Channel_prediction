function [parameters,ret_parameter] = channelEstimation(csi_sample)
%% csi_sample: M*1 complex
% theta is a L*3 vector: [alpha, phi, tau] 
% [alpha, -, -]_l: the amplitude of the path l;
% [-, phi, -]_l: the azimuth angle of the path l;
% [-, -, tau]_l: the delay of the path l;

globals_init();

%  no input csi -> simulation
if nargin <= 0
    [csi_sample, ~] = generate_simulation();
end

% disp(csi_sample);

global ITERATION L DOMAIN_TAU SIMULATION_PHI SIMULATION_TAU SPEED_OF_LIGHT...
    
% parameters: ITERATION+1 cell, inside each cell is a struct 
% (alpha, phi, tau)
parameters = cell(ITERATION+1, 1);

for i = 1:ITERATION+1
    parameters{i} = struct('alpha', zeros(1, L), ...
    'phi', zeros(1, L), 'tau', zeros(1, L));
end
% parameters{1}.tau = zeros(1, L)+...
%     DOMAIN_TAU.start+DOMAIN_TAU.step*round(DOMAIN_TAU.length/2);
parameters{1} = init(csi_sample, parameters{1});
% groundtruth as input
% parameters{1}.tau = SIMULATION_TAU;
% parameters{1}.phi = SIMULATION_PHI;
% parameters{1}.alpha = (1+1j)./(SPEED_OF_LIGHT*SIMULATION_TAU);

%% Iterating
for I = 1:ITERATION
    for K = 1:L
    % Expectation step (E-step)
    X = eStep(parameters{I}, csi_sample);

    % Maximization step (M-step)
    parameters{I+1}.tau(K) = opt_tau(parameters{I}.alpha(K),... 
        parameters{I}.phi(K), X(:,K));
    parameters{I+1}.phi(K) = opt_phi(X(:,K));
    parameters{I+1}.alpha(K) = compute_alpha(parameters{I+1}.tau(K),...
        parameters{I+1}.phi(K), X(:,K));
    end

    %% compute expectation
%     alpha = parameters{I+1}.alpha(K);
%     alpha = repmat(alpha,M,1);
% 
%     phi = parameters{I+1}.phi(K);
%     phi = repmat(phi,M,1);
% 
%     tau = parameters{I+1}.tau(K);
%     tau = repmat(tau,M,1);
% 
%     C_M = repmat(transpose(0:M-1), 1, L);
%     C = exp(1j*2*pi/LAMBDA*D*C_M.*cos(phi));
% 
%     S = sum(alpha.*C.*exp(-1j*2*pi*tau*FREQUENCY),2);
    %disp(csi_sample-S);
end

ret_parameter = parameters{ITERATION+1};

% H = [];
% for G = 1:ITERATION/L  
%     H = [H;parameters{G*L}.tau];
% end
% plot(H);

end




function [csi_sample,csi_hidden] = generate_simulation()
    global SIMULATION_TAU SIMULATION_PHI SPEED_OF_LIGHT M FREQUENCY
    C = compute_C(SIMULATION_PHI);
    ALPHA = (1+1j)./(SPEED_OF_LIGHT*SIMULATION_TAU);
    csi_sample = repmat(ALPHA, M, 1) .* C .* ...
        repmat(exp(-1j*2*pi*SIMULATION_TAU*FREQUENCY), M, 1);
    csi_hidden = csi_sample;
    csi_sample = sum(csi_sample, 2);
    
    
    %csi_sample = [-0.57 + 0.77; 0.11 - 1.65i; 1.43 - 0.65i];

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
    FREQUENCY = 5.24e9;  %unit hz
    SPEED_OF_LIGHT = 3e8;  %unit m/s
    LAMBDA = SPEED_OF_LIGHT/FREQUENCY;
    M = 30;
    L = 8;
    D = LAMBDA/2;
    ITERATION = 800;
    N = 1;
    DOMAIN_TAU = struct('start', 10e-9, 'end', 30e-9, 'step', 1e-9); % unit: s
    DOMAIN_TAU.length = round((DOMAIN_TAU.end - DOMAIN_TAU.start) ...
        / DOMAIN_TAU.step + 1);
    
    DOMAIN_PHI = struct('start', 0, 'end', pi, 'step', pi/100); % unit: radius
    DOMAIN_PHI.length = round((DOMAIN_PHI.end - DOMAIN_PHI.start) ...
        / DOMAIN_PHI.step + 1);
    
    %% simulation
    global SIMULATION_TAU SIMULATION_PHI
    %SIMULATION_TAU = [0, 1.1]*1e-9;
    %SIMULATION_PHI = [1/18,1/9]*pi;
    SIMULATION_TAU = [12 13 15 17 19 20 22 24]*1e-9;
    SIMULATION_PHI = [0.1 0.3 0.4 0.45 0.5 0.6 0.76 0.8]*pi;
   
end
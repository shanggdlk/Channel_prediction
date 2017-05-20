function channelEstimation(csi_sample)
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

global ITERATION L SIMULATION_PHI SIMULATION_TAU ALPHA ...
 EPSILON
    
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
    X = eStep(parameters{I}, csi_sample, K);

    % Maximization step (M-step)
    parameters{I+1}.tau(K) = opt_tau(parameters{I}.phi(K), X);
    parameters{I+1}.phi(K) = opt_phi(parameters{I+1}.tau(K), X);
    parameters{I+1}.alpha(K) = compute_alpha(parameters{I+1}.tau(K),...
        parameters{I+1}.phi(K), X);
    
    parameters{I}.tau(K) = parameters{I+1}.tau(K);
    parameters{I}.phi(K) = parameters{I+1}.phi(K);
    parameters{I}.alpha(K) = parameters{I+1}.alpha(K);
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

ret_parameter = parameters{ITERATION};
if all((ret_parameter.tau - SIMULATION_TAU).^2 < EPSILON)
    fprintf('predicting tau success\n');
else
    error('predicting tau fail\n')
end

if all((ret_parameter.phi - SIMULATION_PHI).^2 < EPSILON)
    fprintf('predicting phi success\n');
else
    error('predicting phi fail\n')
end

if all((ret_parameter.alpha - ALPHA).^2 < EPSILON)
    fprintf('predicting alpha success\n');
else
    error('predicting alpha fail\n')
end

% h = zeros(4,ITERATION);
% for i = 1:ITERATION
%     h(:,i) = parameters{i}.phi;
% end

% plot(h(1,:));
% hold on;
% plot(h(2,:));
% plot(h(3,:));
% plot(h(4,:));
% hold off;


% H = [];
% for G = 1:ITERATION/L  
%     H = [H;parameters{G*L}.tau];
% end
% plot(H);

end




function [csi_sample,csi_hidden] = generate_simulation()
    global SIMULATION_TAU SIMULATION_PHI L SPEED_OF_LIGHT FREQUENCIES F ALPHA
    C = compute_C(SIMULATION_PHI);
    ALPHA = reshape((1+1j)./(SPEED_OF_LIGHT*SIMULATION_TAU), [1, L, 1]);
    CSI_TAU = reshape(SIMULATION_TAU, [1, L, 1]);
    CSI_FREQUENCY = reshape(FREQUENCIES, [1,1,F]);
    
    
    csi_sample = bsxfun(@times, ALPHA, bsxfun(@times, C, ...
        exp(-1j*2*pi.*bsxfun(@times, CSI_TAU, CSI_FREQUENCY))));
    
    csi_hidden = csi_sample;
    csi_sample = squeeze(sum(csi_sample, 2));
    %disp(size(csi_sample));
end




function globals_init
%% physical
% LAMBDA: wavelength of the signal;
% FREQUENCY: frequency of the signal;
% M: antenna array size;
% L: # propagation paths;
% D: spacing between adjacent rx antenna
% N: # sample
% F: # measured subcarrier
% matrix order: [M, L, F]
% DELTA_FREQUENCY: difference between adjacent subcarrier
    global CENTRAL_FREQUENCY SPEED_OF_LIGHT SHAPE_M M L D ITERATION DOMAIN_TAU ...
        DOMAIN_PHI F DELTA_FREQUENCY FREQUENCIES LAMBDAS EPSILON
    CENTRAL_FREQUENCY = 2.462e9;  %unit hz
 
    F = 64;   % # of frequency 
    DELTA_FREQUENCY = 20e6/F; % 20Mhz, 64 slots
  
    SPEED_OF_LIGHT = 3e8;  %unit m/s
    
    FREQUENCIES = ((1:F) - (F+1)/2) * DELTA_FREQUENCY+CENTRAL_FREQUENCY;
    %disp(FREQUENCIES(1,2)-FREQUENCIES(1,1));
    LAMBDAS = SPEED_OF_LIGHT./FREQUENCIES;
    
    %shape of the 2d antenna array, in this case, 2 row * 6 column
    SHAPE_M = [2 6];
    M = prod(SHAPE_M);
    L = 2;
    %D = mean(LAMBDAS)/2;
    D = 0.031;
    ITERATION = 50;
    EPSILON = 0.1;  %required accuracy
    % search space
    DOMAIN_TAU = struct('start', 10e-9, 'end', 30e-9, 'step', 1e-9); % unit: s
    DOMAIN_TAU.length = round((DOMAIN_TAU.end - DOMAIN_TAU.start) ...
        / DOMAIN_TAU.step + 1);
    
    DOMAIN_PHI = struct('start', 0, 'end', pi, 'step', pi/100); % unit: radius
    DOMAIN_PHI.length = round((DOMAIN_PHI.end - DOMAIN_PHI.start) ...
        / DOMAIN_PHI.step + 1);
    
    %% simulation
    global SIMULATION_TAU SIMULATION_PHI
    %SIMULATION_TAU = [12 13 15 17 19 20 22 24]*1e-9;
    %SIMULATION_PHI = [0.1 0.3 0.4 0.45 0.5 0.6 0.76 0.8]*pi;
    %SIMULATION_TAU = [12 13 15 17]*1e-9;
    %SIMULATION_PHI = [0.1 0.3 0.4 0.6]*pi;
    SIMULATION_PHI = [0.1 0.3]*pi;
    SIMULATION_TAU = [12 13]*1e-9;
end
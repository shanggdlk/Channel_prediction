function [parameters,ret_parameter] = channelEstimation(csi_sample)
%% csi_sample: M*1 complex
% theta is a L*3 vector: [alpha, phi, tau] 
% [alpha, -, -]_l: the amplitude of the path l;
% [-, phi, -]_l: the azimuth angle of the path l;
% [-, -, tau]_l: the delay of the path l;

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
    X = eStep(parameters{I}, csi_sample, K);

    % Maximization step (M-step)
    parameters{I+1}.tau(K) = opt_tau(parameters{I}.phi(K), X);
    parameters{I+1}.phi(K) = opt_phi(parameters{I+1}.tau(K), X);
    parameters{I+1}.alpha(K) = compute_alpha(parameters{I+1}.tau(K),...
        parameters{I+1}.phi(K), X);
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
    global SIMULATION_TAU SIMULATION_PHI L SPEED_OF_LIGHT M FREQUENCIES F
    C = compute_C(SIMULATION_PHI);
    ALPHA = repmat((1+1j)./(SPEED_OF_LIGHT*SIMULATION_TAU), M, 1, F);
    CSI_TAU = repmat(SIMULATION_TAU, M, 1, F);
    CSI_FREQUENCY = repmat(reshape(FREQUENCIES, [1,1,F]), M, L);
    csi_sample = ALPHA .* C .* exp(-1j*2*pi.*CSI_TAU.*CSI_FREQUENCY);
    csi_hidden = csi_sample;
    csi_sample = squeeze(sum(csi_sample, 2));
end

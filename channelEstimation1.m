function [parameters,ret_parameter] = channelEstimation1(csi_sample)
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

disp(csi_sample);

global ITERATION L SIMULATION_PHI SIMULATION_TAU SPEED_OF_LIGHT

% parameters: ITERATION+1 cell, inside each cell is a struct 
% (alpha, phi, tau)
parameters = cell(ITERATION+1, 1);

for i = 1:ITERATION+1
    parameters{i} = struct('alpha', zeros(1, L), ...
    'phi', zeros(1, L), 'tau', zeros(1, L));
end

% parameters{1} = init1(csi_sample, parameters{1});
parameters{1}.tau = SIMULATION_TAU;
parameters{1}.phi = SIMULATION_PHI;
parameters{1}.alpha = (1+1j)./(SPEED_OF_LIGHT*SIMULATION_TAU);

%% Iterating;
for I = 1:ITERATION
    
    K = mod(I-1,L)+1;
    % Expectation step (E-step)
    X = eStep(parameters{I}, csi_sample);

    % Maximization step (M-step)
    tau = opt_tau(parameters{I}.phi(K), X(:,K));
    phi = opt_phi(tau, X(:,K));
    alpha = compute_alpha(tau, phi, X(:,K));
    
    parameters{I+1} = parameters{I};
    parameters{I+1}.tau(K) = tau;
    parameters{I+1}.phi(K) = phi;
    parameters{I+1}.alpha(K) = alpha;
%     
%     
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

function ret_tau = opt_tau(phi, X)

global DOMAIN_TAU M LAMBDA D FREQUENCY
% 
C_M = transpose(0:M-1);
C_L = repmat(cos(phi), M, 1);
C = exp(1j*2*pi/LAMBDA*D*C_M.*C_L);
C = transpose(ctranspose(C));
C = repmat(C,1,DOMAIN_TAU.length);


X = repmat(X,1,DOMAIN_TAU.length);

tau_space = repmat((DOMAIN_TAU.start:DOMAIN_TAU.step:DOMAIN_TAU.end),M,1);

Z_abs = real(sum(C .* X .* exp(1j*2*pi*FREQUENCY.*tau_space),1));
[~, I] = max(Z_abs);

ret_tau = DOMAIN_TAU.step*(I-1)+DOMAIN_TAU.start;

% %
% tau_space = zeros(DOMAIN_TAU.length, 1);
% for i = 1:DOMAIN_TAU.length
%     tau_space(i,:) = compute_Z1(zeros(1, 1)+DOMAIN_TAU.step*(i-1)...
%         +DOMAIN_TAU.start, phi, X);
% end
% 
% [~, I] = max(abs(tau_space));
% ret_tau = DOMAIN_TAU.step*(I-1)+DOMAIN_TAU.start;

end


function ret_phi = opt_phi(tau, X)

global DOMAIN_PHI M LAMBDA D FREQUENCY

C_M = repmat(transpose(0:M-1),1,DOMAIN_PHI.length);
phi_space = repmat((DOMAIN_PHI.start:DOMAIN_PHI.step:DOMAIN_PHI.end),M,1);
C_L = cos(phi_space);
C = exp(1j*2*pi/LAMBDA*D*C_M.*C_L);
C = transpose(ctranspose(C));

X = repmat(X,1,DOMAIN_PHI.length);

tau_space = repmat(tau,M,DOMAIN_PHI.length);

Z_abs = abs(sum(C .* X .* exp(1j*2*pi*FREQUENCY.*tau_space),1));
[~, I] = max(Z_abs);

ret_phi = DOMAIN_PHI.step*(I-1)+DOMAIN_PHI.start;

%global DOMAIN_PHI

% phi_space = zeros(DOMAIN_PHI.length, 1);
% for i = 1:DOMAIN_PHI.length
%     phi_space(i,:) = compute_Z1(tau, zeros(1, 1)+DOMAIN_PHI.step*(i-1)...
%         +DOMAIN_PHI.start, X);
% end
% [~, I] = max(abs(phi_space));
% ret_phi = DOMAIN_PHI.step*(I-1)+DOMAIN_PHI.start;
end

function ret_alpha = compute_alpha(tau, phi, X)

global M N LAMBDA D FREQUENCY

C_M = transpose(0:M-1);
phi_space = repmat(phi,M,1);
C_L = cos(phi_space);
C = exp(1j*2*pi/LAMBDA*D*C_M.*C_L);
C = transpose(ctranspose(C));

tau_space = repmat(tau,M,1);

Z = sum(C .* X .* exp(1j*2*pi*FREQUENCY.*tau_space),1);
ret_alpha = 1/(M*N)*Z;
% 
% ret_alpha = 1/(M*N)*(compute_Z1(tau, phi, X));

end


function csi_sample = generate_simulation()
    global SIMULATION_TAU SIMULATION_PHI SPEED_OF_LIGHT M FREQUENCY
    C = compute_C(SIMULATION_PHI);
    ALPHA = (1+1j)./(SPEED_OF_LIGHT*SIMULATION_TAU);
    csi_sample = repmat(ALPHA, M, 1) .* C .* ...
        repmat(exp(-1j*2*pi*SIMULATION_TAU*FREQUENCY), M, 1);
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
    M = 10;
    L = 8;
    D = LAMBDA/2;
    ITERATION = 800;
    N = 2;
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
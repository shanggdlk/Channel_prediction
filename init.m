function [parameter] = init(csi, parameter)

global N L M FREQUENCY DOMAIN_TAU DOMAIN_PHI

for K = 1:L
%% init tau
    X = eStep(parameter, csi);
    if K ~= L
        X(:,K+1:end) = 0;
    end

    tau_space = repmat(X, 1, DOMAIN_TAU.length);
    var = reshape(repmat(DOMAIN_TAU.start:DOMAIN_TAU.step:DOMAIN_TAU.end, ...
    L*M, 1), [M, DOMAIN_TAU.length*L]);
    %disp(tau_space .* exp(-1j*2*pi*FREQUENCY*var));
    
    tau_space = abs(tau_space .* exp(-1j*2*pi*FREQUENCY*var));
   
    %tau_space = tau_space .* exp(-1j*2*pi*FREQUENCY*var);

    opt_tau = transpose(reshape(sum(tau_space, 1), [L, DOMAIN_TAU.length]));

    
    [~, I] = max(opt_tau);

    parameter.tau = DOMAIN_TAU.step*(I-1)+DOMAIN_TAU.start;

    %% init phi
    DOMAIN_PHI.start = 0; DOMAIN_PHI.end = pi; DOMAIN_PHI.step = pi/10;  
    DOMAIN_PHI.length = (DOMAIN_PHI.end - DOMAIN_PHI.start) / DOMAIN_PHI.step + 1;

    phi_space = zeros(DOMAIN_PHI.length, L);
    for i = 1:DOMAIN_PHI.length
    phi_space(i,:) = compute_Z(parameter.tau,... 
        zeros(1, L)+DOMAIN_PHI.step*(i-1)+DOMAIN_PHI.start, X);
    end
    [~, I] = max(abs(phi_space));
    parameter.phi = DOMAIN_PHI.step*(I-1)+DOMAIN_PHI.start;
    

    %% init alpha
    parameter.alpha = 1/(M*N)*abs(compute_Z(parameter.tau, parameter.phi, X));
    
    if K ~= L
        parameter.alpha(K+1:end) = 0;
        parameter.phi(K+1:end) = 0;
        parameter.tau(K+1:end) = 0;
    end
    
    if K == 1
        old_parameter = parameter;
    else
        old_parameter.alpha(K) = parameter.alpha(K);
        old_parameter.tau(K) = parameter.tau(K);
        old_parameter.phi(K) = parameter.phi(K);
    end
    parameter = old_parameter;
    %disp(old_parameter);
end

% 
% %% init tau
% X = eStep(parameter, csi);
% %disp(X);
% tau_space = repmat(X, 1, DOMAIN_TAU.length);
% 
% var = reshape(repmat(DOMAIN_TAU.start:DOMAIN_TAU.step:DOMAIN_TAU.end, ...
%     L*M, 1), [M, DOMAIN_TAU.length*L]);
% 
% tau_space = abs(tau_space .* exp(-1j*2*pi*FREQUENCY*var));
% %tau_space = tau_space .* exp(-1j*2*pi*FREQUENCY*var);
% 
% opt_tau = transpose(reshape(sum(tau_space, 1), [L, DOMAIN_TAU.length]));
% 
% 
% 
% [~, I] = max(opt_tau);
% 
% parameter.tau = DOMAIN_TAU.step*(I-1)+DOMAIN_TAU.start;
% 
% %% init phi
% DOMAIN_PHI.start = 0; DOMAIN_PHI.end = pi; DOMAIN_PHI.step = pi/10;  
% DOMAIN_PHI.length = (DOMAIN_PHI.end - DOMAIN_PHI.start) / DOMAIN_PHI.step + 1;
% 
% phi_space = zeros(DOMAIN_PHI.length, L);
% for i = 1:DOMAIN_PHI.length
%     phi_space(i,:) = compute_Z(parameter.tau,... 
%         zeros(1, L)+DOMAIN_PHI.step*(i-1)+DOMAIN_PHI.start, X);
% end
% [~, I] = max(abs(phi_space));
% parameter.phi = DOMAIN_PHI.step*(I-1)+DOMAIN_PHI.start;
% 
% 
% %% init alpha
% parameter.alpha = 1/(M*N)*abs(compute_Z(parameter.tau, parameter.phi, X));

<<<<<<< HEAD
=======
%% init alpha
parameter.alpha = 1/(M*N)*compute_Z(parameter.tau, parameter.phi, X);
        
>>>>>>> 4507a2069eeedf937fbc64f2742b032f65c6c57d
end
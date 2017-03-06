function [parameter] = init1(csi, parameter)

global L M D FREQUENCY DOMAIN_TAU DOMAIN_PHI LAMBDA

for K = 1:L
    
    %% compute the X_k;
    phi = parameter.phi;
    phi = repmat(phi,M,1);
    
    tau = parameter.tau;
    tau = repmat(tau,M,1);
    
    alpha = parameter.alpha;
    alpha = repmat(alpha,M,1);
    
    C_M = repmat(transpose(0:M-1), 1, L);
    
    C_L = cos(phi);
    
    C = exp(1j*2*pi/LAMBDA*D*C_M.*C_L);    
    
    S_matrix = alpha.* C .* exp(-1j*2*pi*FREQUENCY*tau);
    
    X_k = csi - sum(S_matrix,2) + S_matrix(:,K);
    % TODO init tau
    
    % init phi
    parameter.phi(K) = opt_phi(X_k);
    
    % init alpha
    parameter.alpha(K) = compute_alpha(parameter.tau(K), parameter.phi(K), X_k);
end

%disp(parameter.alpha);
%disp(parameter.phi);
%disp(parameter.tau);

%% testing whether the initialization works correctly.
alpha = parameter.alpha;
alpha = repmat(alpha,M,1);

phi = parameter.phi;
phi = repmat(phi,M,1);

tau = parameter.tau;
tau = repmat(tau,M,1);

C_M = repmat(transpose(0:M-1), 1, L);
C = exp(1j*2*pi/LAMBDA*D*C_M.*cos(phi));

S = sum(alpha.*C.*exp(-1j*2*pi*tau*FREQUENCY),2);
disp(S);

end
function [parameter] = init1(csi, parameter)

global N L M D FREQUENCY DOMAIN_TAU DOMAIN_PHI LAMBDA

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
    %% init tau
    tau_space = repmat(X_k, 1, DOMAIN_TAU.length);
    var = (DOMAIN_TAU.start:DOMAIN_TAU.step:DOMAIN_TAU.end);
    var = repmat(var,M,1);
    
    tau_space = abs(tau_space .* exp(1j*2*pi*FREQUENCY*var));
    %disp(tau_space);
    res = sum(tau_space,1);
    %disp(res);
    I=find(res==max(res));
    I = I(1);
    opt_tau = DOMAIN_TAU.step*(I-1)+DOMAIN_TAU.start;
    %disp(parameter.tau);
    parameter.tau(K) = opt_tau;
    
    %% init phi
    
    phi_matrix = repmat((DOMAIN_PHI.start:DOMAIN_PHI.step:DOMAIN_PHI.end),M,1);
    
    C_M = repmat(transpose(0:M-1), 1, DOMAIN_PHI.length);
    C = exp(1j*2*pi/LAMBDA*D*C_M.*cos(phi_matrix));
    C = transpose(ctranspose(C));
    
    opt_tau_matrix = repmat(opt_tau, M, DOMAIN_PHI.length);
    Z = sum(C .* repmat(X_k, 1, DOMAIN_PHI.length) .* exp(1j*2*pi*FREQUENCY*opt_tau_matrix),1);
    Z_val = abs(Z)
    I=find(Z_val==max(Z_val));
    I = I(1);
    opt_phi = DOMAIN_PHI.step*(I-1)+DOMAIN_PHI.start;
    parameter.phi(K) = opt_phi;
    
    %% init alpha
    opt_tau = repmat(opt_tau,M,1);
    opt_phi = repmat(opt_phi,M,1);
    C_M = transpose(0:M-1);
    C = exp(1j*2*pi/LAMBDA*D*C_M.*cos(opt_phi));
    C = transpose(ctranspose(C));
    %disp(size(C));
    Z = sum(C .* X_k .* exp(1j*2*pi*FREQUENCY*opt_tau),1);
    alpha = 1/(M*N)*(Z);
    parameter.alpha(K) = alpha;
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
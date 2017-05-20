function [parameter] = init(csi, parameter)

global L FREQUENCIES F

S_FREQUENCY = reshape(FREQUENCIES, [1,1,F]);
for K = 1:L
    
    %% compute the X_k;
    alpha = reshape(parameter.alpha, [1, L, 1]);
    
    C = compute_C(parameter.phi);
    
    S_TAU = reshape(parameter.tau, [1, L, 1]);
    
    S_matrix = bsxfun(@times, alpha, bsxfun(@times, C, ...
        exp(-1j*2*pi.*bsxfun(@times, S_FREQUENCY, S_TAU))));
    
    X_k = csi - squeeze(sum(S_matrix,2)) + squeeze(S_matrix(:,K,:));
    
    % init tau
    parameter.tau(K) = init_tau(X_k);

    % init phi
    parameter.phi(K) = opt_phi(parameter.tau(K), X_k);
   
    % init alpha
    parameter.alpha(K) = compute_alpha(parameter.tau(K), parameter.phi(K), X_k);
end

% disp(parameter.alpha);
% disp(parameter.phi);

%% testing whether the initialization works correctly.
% alpha = parameter.alpha;
% alpha = repmat(alpha,M,1);
% 
% phi = parameter.phi;
% phi = repmat(phi,M,1);
% 
% tau = parameter.tau;
% tau = repmat(tau,M,1);
% 
% C_M = repmat(transpose(0:M-1), 1, L);
% C = exp(1j*2*pi/LAMBDA*D*C_M.*cos(phi));
% 
% S = sum(alpha.*C.*exp(-1j*2*pi*tau*FREQUENCY),2);
% disp(S);

end
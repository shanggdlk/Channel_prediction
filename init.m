function [parameter] = init(csi, parameter)

global N L M FREQUENCY DOMAIN_TAU DOMAIN_PHI


%% init tau
X = eStep(parameter, csi);

tau_space = repmat(X, 1, DOMAIN_TAU.length);

var = reshape(repmat(DOMAIN_TAU.start:DOMAIN_TAU.step:DOMAIN_TAU.end, ...
    L*M, 1), [M, DOMAIN_TAU.length*L]);

tau_space = abs(tau_space .* exp(-1j*2*pi*FREQUENCY*var));

opt_tau = transpose(reshape(sum(tau_space, 1), [L, DOMAIN_TAU.length]));

[~, I] = max(opt_tau);

parameter.tau = DOMAIN_TAU.step*(I-1)+DOMAIN_TAU.start;

%% init phi
phi_space = zeros(DOMAIN_PHI.length, L);
for i = 1:DOMAIN_PHI.length
    phi_space(i,:) = compute_Z(parameter.tau,... 
        zeros(1, L)+DOMAIN_PHI.step*(i-1)+DOMAIN_PHI.start, X);
end
[~, I] = max(abs(phi_space));
parameter.phi = DOMAIN_PHI.step*(I-1)+DOMAIN_PHI.start;


%% init alpha
parameter.alpha = 1/(M*N)*abs(compute_Z(parameter.tau, parameter.phi, X));
        
end
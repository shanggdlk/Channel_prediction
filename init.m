function [parameter] = init(csi, parameter)

global N L M FREQUENCY


%% init tau
X = eStep(parameter,csi);

start_tau = 10; end_tau = 30;step_tau = 1;  % unit: ns
n_tau = (end_tau - start_tau) / step_tau + 1;

tau_space = repmat(X, 1, n_tau);

var = reshape(repmat(start_tau:step_tau:end_tau, L*M, 1), [M, n_tau*L]);

tau_space = abs(tau_space .* exp(-1j*2*pi*FREQUENCY*var));

opt_tau = transpose(reshape(sum(tau_space, 1), [L, n_tau]));

[~, I] = max(opt_tau);

parameter.tau = step_tau*(I-1)+start_tau;

%% init phi
start_phi = 0; end_phi = pi; step_phi = pi/10;  % unit: radius
n_phi = (end_phi - start_phi) / step_phi + 1;

phi_space = zeros(n_phi, L);
for i = 1:n_phi
    phi_space(i,:) = compute_Z(parameter.tau,... 
        zeros(1, L)+step_phi*(i-1)+start_phi, X);
end
[~, I] = max(abs(phi_space));
parameter.phi = step_phi*(I-1)+start_phi;


%% init alpha
parameter.alpha = 1/(M*N)*abs(compute_Z(parameter.tau, parameter.phi, X));
        
end
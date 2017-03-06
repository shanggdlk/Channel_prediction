function ret_phi = opt_phi(X)

global DOMAIN_PHI M LAMBDA D

C_M = repmat(transpose(0:M-1),1,DOMAIN_PHI.length);
phi_space = repmat((DOMAIN_PHI.start:DOMAIN_PHI.step:DOMAIN_PHI.end),M,1);
C_L = cos(phi_space);
C = exp(1j*2*pi/LAMBDA*D*C_M.*C_L);
C = transpose(ctranspose(C));

X = repmat(X,1,DOMAIN_PHI.length);

% no need to compute tau
% tau_space = repmat(tau,M,DOMAIN_PHI.length);
% Z_abs = abs(sum(C .* X .* exp(1j*2*pi*FREQUENCY.*tau_space),1));
Z_abs = abs(sum(C .* X,1));
[~, I] = max(Z_abs);

ret_phi = DOMAIN_PHI.step*(I-1)+DOMAIN_PHI.start;

end
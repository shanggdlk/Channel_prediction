function ret_phi = opt_phi(tau, X)

global DOMAIN_PHI M LAMBDAS D F FREQUENCIES

C_M = repmat(transpose(0:M-1),1,DOMAIN_PHI.length, F);
phi_space = repmat((DOMAIN_PHI.start:DOMAIN_PHI.step:DOMAIN_PHI.end),M,1, F);
C_L = cos(phi_space);
C_LAMBDA = repmat(reshape(LAMBDAS, [1, 1, F]), M, DOMAIN_PHI.length, 1);
C = exp(1j*2*pi./C_LAMBDA*D.*C_M.*C_L); 

C = conj(C);

X = repmat(reshape(X,[M,1,F]),1,DOMAIN_PHI.length, 1);

Z_FREQUENCY = repmat(reshape(FREQUENCIES, [1,1,F]), M, DOMAIN_PHI.length);
Z_TAU = repmat(tau,M,DOMAIN_PHI.length, F);
Z_abs = abs(squeeze(sum(sum(C .* X .* ...
    exp(1j*2*pi.*Z_FREQUENCY.*Z_TAU),1), 3)));
[~, I] = max(Z_abs);

ret_phi = DOMAIN_PHI.step*(I-1)+DOMAIN_PHI.start;

end
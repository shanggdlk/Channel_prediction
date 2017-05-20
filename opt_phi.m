function ret_phi = opt_phi(tau, X)

global DOMAIN_PHI M F FREQUENCIES


phi_space = DOMAIN_PHI.start:DOMAIN_PHI.step:DOMAIN_PHI.end;
C = conj(compute_C(phi_space));

X = reshape(X,[M,1,F]);

Z_FREQUENCY = reshape(FREQUENCIES, [1,1,F]);

Z_abs = abs(squeeze(sum(sum(bsxfun(@times, C, bsxfun(@times, X, ...
    exp(1j*2*pi.*Z_FREQUENCY.*tau))),1), 3)));
[~, I] = max(Z_abs);

ret_phi = DOMAIN_PHI.step*(I-1)+DOMAIN_PHI.start;

end
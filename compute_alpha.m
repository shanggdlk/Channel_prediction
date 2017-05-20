function ret_alpha = compute_alpha(tau, phi, X)

global M F FREQUENCIES

C = compute_C(phi);
C = squeeze(conj(C));

Z = (bsxfun(@times, C.*X, exp(1j*2*pi.*reshape(FREQUENCIES, [1 F]).*tau)));
ret_alpha = 1/(M*F)*sum(Z(:));

end
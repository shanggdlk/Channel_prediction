function ret_tau = opt_tau(phi, X)

global DOMAIN_TAU M FREQUENCIES L F

% only need a slice of C
C = conj(compute_C(phi));

X = reshape(X, [M, 1, F]);

Z_FREQUENCY = reshape(FREQUENCIES, [1,1,F]);
Z_TAU = DOMAIN_TAU.start:DOMAIN_TAU.step:DOMAIN_TAU.end;

Z_abs = bsxfun(@times, C, bsxfun(@times, X, exp(1j*2*pi.*bsxfun(@times, Z_FREQUENCY, Z_TAU))));
Z_abs = squeeze(sum(sum(Z_abs, 3), 1));
[~, I] = max(Z_abs);

ret_tau = DOMAIN_TAU.step*(I-1)+DOMAIN_TAU.start;

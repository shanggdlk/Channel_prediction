function ret = compute_Z(tau, phi, X)
% tau: 1*L
% phi: 1*L
% X: M*L
global FREQUENCY M


C = compute_C(phi);
C = transpose(ctranspose(C));

ret = C .* X .* repmat(exp(1j*2*pi*FREQUENCY.*tau), M, 1);
ret = sum(ret, 1);

end
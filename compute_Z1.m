function ret = compute_Z1(tau, phi, X)
% tau: 1*L
% phi: 1*L
% X: M*L
global FREQUENCY M


C = compute_C1(phi);
C = transpose(ctranspose(C));
Z = C .* X .* repmat(exp(1j*2*pi*FREQUENCY.*tau), M, 1);
disp(size(Z));
end
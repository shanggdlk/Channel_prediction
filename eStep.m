function X = eStep(parameter, csi)

global M L LAMBDA D FREQUENCY


C_M = repmat(transpose(0:M-1), 1, L);
C_L = repmat(cos(parameter.phi), M, 1);
C = exp(-1j*2*pi/LAMBDA*D*C_M.*C_L);

S = repmat(parameter.alpha, M, 1) .* C .* ...
    exp(-1j*2*pi*FREQUENCY*repmat(parameter.tau, M, 1));

X = repmat(csi - sum(S, 2), 1,L) + S;

end

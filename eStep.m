function X = eStep(parameter, csi)

global M L FREQUENCY

C = compute_C(parameter.phi);

S = C .* ...
    exp(-1j*2*pi*FREQUENCY*repmat(parameter.tau, M, 1));
X = repmat(csi - sum(S, 2), 1,L) + S;
end


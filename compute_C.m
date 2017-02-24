function ret = compute_C(phi)

global M L LAMBDA D

C_M = repmat(transpose(0:M-1), 1, L);
C_L = repmat(cos(phi), M, 1);
ret = exp(1j*2*pi/LAMBDA*D*C_M.*C_L);

end
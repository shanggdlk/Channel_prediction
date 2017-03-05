function ret = compute_C1(phi)

global M LAMBDA D

C_M = transpose(0:M-1);
C_L = repmat(cos(phi), M, 1);
ret = exp(1j*2*pi/LAMBDA*D*C_M.*C_L);
end
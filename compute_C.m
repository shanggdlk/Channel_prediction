function ret = compute_C(phi)
%M * L * F
global M L LAMBDAS D F

C_M = repmat(transpose(0:M-1), 1, L, F);
C_L = repmat(cos(phi), M, 1, F);
C_LAMBDA = repmat(reshape(LAMBDAS, [1,1,F]), M, L, 1);
ret = exp(1j*2*pi./C_LAMBDA*D.*C_M.*C_L);
end
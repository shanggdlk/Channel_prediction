function ret_alpha = compute_alpha(tau, phi, X)

global M N LAMBDA D FREQUENCY

C_M = transpose(0:M-1);
phi_space = repmat(phi,M,1);
C_L = cos(phi_space);
C = exp(1j*2*pi/LAMBDA*D*C_M.*C_L);
C = ctranspose(C);


Z = sum(C * X * exp(1j*2*pi*FREQUENCY*tau),1);
ret_alpha = 1/(M*N)*Z;

end
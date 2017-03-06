function ret_tau = opt_tau(alpha, phi, X)

global DOMAIN_TAU M LAMBDA D FREQUENCY
% 
C_M = transpose(0:M-1);
C_L = repmat(cos(phi), M, 1);
C = exp(1j*2*pi/LAMBDA*D*C_M.*C_L);
C = transpose(ctranspose(C));
C = repmat(C,1,DOMAIN_TAU.length);


X = repmat(X,1,DOMAIN_TAU.length);

tau_space = repmat((DOMAIN_TAU.start:DOMAIN_TAU.step:DOMAIN_TAU.end),M,1);

Z_abs = real(sum(C .* X .* exp(1j*2*pi*FREQUENCY.*tau_space)*...
    ctranspose(alpha),1));
[~, I] = max(Z_abs);

ret_tau = DOMAIN_TAU.step*(I-1)+DOMAIN_TAU.start;

% %
% tau_space = zeros(DOMAIN_TAU.length, 1);
% for i = 1:DOMAIN_TAU.length
%     tau_space(i,:) = compute_Z1(zeros(1, 1)+DOMAIN_TAU.step*(i-1)...
%         +DOMAIN_TAU.start, phi, X);
% end
% 
% [~, I] = max(abs(tau_space));
% ret_tau = DOMAIN_TAU.step*(I-1)+DOMAIN_TAU.start;

end
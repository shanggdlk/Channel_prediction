function X = eStep(theta,freq, L, M, csi, lambda,d)
csi = csi';
% csi is a vector;
% theta = [alpha, azimuthal_angle, tau];
% C = zeros(M,L);

% S: MxL matrix;
S = zeros(M,L);

for a = 1:M
    
    C = zeros(1,L);
    G = zeros(1,L);
    alpha = zeros(1,L);
    % X = alpha*C*G;
    
    for b = 1:L
        
        alpha(b) = theta(b,1);
        
        a_angle = theta(b,2);
        C(b) = exp((1i)*(2*pi/lambda)*(a-1)*d*cos(a_angle));
        
        tau = theta(b,3);
        G(b) = exp(-(1i)*2*pi*tau*freq);
    end
    
    alpha = repmat(alpha,L,1);
    C = repmat(C,L,1);
    G = repmat(G,L,1);
    G = G.';
    C = C-C.*eye(size(C));
    alpha = alpha - alpha.*eye(size(alpha));
    % has the following code or do not have has the same effect;
    %G = G - G.*eye(size(G));
    S_temp = alpha.*C*G;
    S(a,:) = diag(S_temp).';
end

H = repmat(csi,L,1);
H = H.';
X = H - S;

end

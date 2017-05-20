function [csi_sample] = compute_CSI(PHI,TAU,ALPHA)
    M = 7;
    L = 2;
    D = 0.031;
    CENTRAL_FREQUENCY = 2.462e9;  %unit hz
    DELTA_FREQUENCY = 20e6/64; % 20Mhz, 64 slots
    F = 64;    
    SPEED_OF_LIGHT = 3e8;  %unit m/s    
    FREQUENCIES = ((1:F) - (F+1)/2) * DELTA_FREQUENCY+CENTRAL_FREQUENCY;
    LAMBDAS = SPEED_OF_LIGHT./FREQUENCIES;
    
    C_M = repmat(transpose(0:M-1), 1, L, F);
    C_L = repmat(cos(PHI), M, 1, F);
    C_LAMBDA = repmat(reshape(LAMBDAS, [1,1,F]), M, L, 1);
    C = exp(1j*2*pi./C_LAMBDA*D.*C_M.*C_L);
    
    ALPHA = repmat(ALPHA, M, 1, F);
    CSI_TAU = repmat(TAU, M, 1, F);
    CSI_FREQUENCY = repmat(reshape(FREQUENCIES, [1,1,F]), M, L);
    %size(ALPHA)
    csi_sample = ALPHA .* C .* exp(-1j*2*pi.*CSI_TAU.*CSI_FREQUENCY);
    %csi_sample = C .* exp(-1j*2*pi.*CSI_TAU.*CSI_FREQUENCY);
    csi_sample = squeeze(sum(csi_sample, 2));
    %disp(size(csi_sample));
    csi_sample(:,28:38) = [];
    csi_sample(:,1) = [];
end
function ret_tau = init_tau(X)

global DOMAIN_TAU M FREQUENCIES F

tau_X = reshape(X, [M, 1, F]);

tau_FREQUENCY = reshape(FREQUENCIES, [1,1,F]);
tau_count = reshape(DOMAIN_TAU.start:DOMAIN_TAU.step:DOMAIN_TAU.end,...
    [1,DOMAIN_TAU.length,1]);

Tau_val = squeeze(sum(abs(sum(bsxfun(@times, tau_X,... 
    exp(1j*2*pi.*bsxfun(@times, tau_FREQUENCY, tau_count))),3)),1));
[~, I] = max(Tau_val);

ret_tau = DOMAIN_TAU.step*(I-1)+DOMAIN_TAU.start;

end
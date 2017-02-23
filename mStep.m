function alpha_next = mStep (X, freq, L, M, lambda,d)

% alpha_next: store the new channel parameter computed via EM algorithm;
alpha_next = zeros(L,3);

for a = 1:L
    X_l = X(:,a);
    parameter = MLE (X_l, freq, M, lambda,d);
    alpha_next(a,:) = parameter;
end

end
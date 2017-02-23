function parameter = sage(X, freq, M, L, lambda,d,theta)

% theta_new: store the next theta estimation: (amplitude, angle, delay);
% theta_new: Lx3;
theta_new = zeros(L,3);

for a = 1:L
    % X_l: Mx1 vector;
    X_l = X(:,a);
    % theta_in: 1X3 vector;
    theta_old = theta(a,:);
    parameter = MLE1 (X_l, theta_old, freq, M, lambda,d);
    theta_new(a,:) = parameter;
end

end
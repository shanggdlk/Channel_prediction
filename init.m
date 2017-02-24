function [parameter] = init(csi, parameter)

global N L M FREQUENCY

% X: M*L

X = eStep(parameter,csi);

start_tau = 10; end_tau = 30;step_tau = 1;
n_tau = (end_tau - start_tau) / step_tau + 1;

space = repmat(X, 1, n_tau);

var = reshape(repmat(start_tau:step_tau:end_tau, L*M, 1), [M, n_tau*L]);

space = space .* exp(-1j*2*pi*FREQUENCY*var);

opt_tau = transpose(reshape(sum(space, 1), [L, n_tau]));

[~, I] = max(opt_tau);

parameter.tau = step_tau*(I-1)+start_tau;

for a = 1:L

    % X_in: 1XM vector;
    X_in = X(:,a).';
    max_val = 0;
    
    % find the optimal delay;
    for b = 1:50
        dnm = b/(1e-9);
        c_val = exp((1i)*2*pi*freq*dnm);
        
        % c_in: 1XM vector;
        c_in = (c_val*ones(M,1))';
        value = sum(abs((c_in.*X_in)));
        
        % record the optimal delay;
        if value > max_val
            max_val = value;
            parameter(a,3) = dnm;
        end
        
    end
    
    % angle_out: 1XL;
    % c: Mx1;
    
    c = zeros(M,1);
    max_z = 0;
    
    % find the optimal azimuthal angle;
    for angle = 1:pi/360:pi
    
        for k = 1:M
            c(k) = exp((1i)*(2*pi/lambda)*cos(angle)*d*(k-1)); 
        end
        
        % c = c^H;
        c = c';
        % c: 1xM;
        % g: 1x1;
        g = exp(-(1i)*(2*pi)*parameter(a,3)*freq);
        
        % z: 1x1 = (1XM)*(M*1)*1;
        z = c*(X_in.')*g;
        
        if (abs(z) > parameter(a,2))
            max_z = abs(z);
            parameter(a,2) = abs(z);            
        end
        
    end 
    
    
% compute the amplitude;
c = zeros(M,1);

for k = 1:M
    c(k) = exp((1i)*(2*pi/lambda)*cos(parameter(a,2))*d*(k-1)); 
end

% c=c^H;
c = c';

g = exp(-(1i)*(2*pi)*parameter(a,3)*freq);
z = c*X(:,a)*g;

end

% parameter: Lx3
parameter(a,1) = (1/(M*N))*z;

end
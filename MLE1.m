function parameter = MLE1 (X, theta, freq, M, lambda,d)
% X: Mx1 vector;
% theta: 1x3 vector;

N = 1;
% N is the number of samples in the freq domain;

% parameter stores the estimated channel parameter (amplitude, angle, delay);
parameter = zeros(1,3);
% max_z: denote the maximum z;
max_z = 0;
% an_angle: denote the azimuthal angle;
a_angle = theta(1,2);


c = zeros(1,M);

for k = 1:M
    c(k) = exp((1i)*(2*pi/lambda)*cos(a_angle)*d*(k-1)); 
end

for delay = 1:50
    % searching space of delay range from 1ns to 50ns;
    
    dnm = delay/(1e-9);
    
    g = exp(-(1i)*(2*pi)*dnm*freq);
    
    % c(1xM) * X(Mx1) * g(1);
    % z = 1x1 (a+bj);
    z = c*X*g;
    
    if (abs(z) > max_z)
        max_z = abs(z);
        % store the maximum delay;
        parameter(1,3) = dnm;
    end 
end

% argmax azimuth angle;

% archive the optimal delay;
dnm = parameter(1,3);

g = exp(-(1i)*(2*pi)*dnm*freq);

% store the maximum z;
max_z = 0;

% c: 1xM;
c = zeros(1,M);

for angle = 0:pi/360:pi
    
    % c: 1xM vector;
    for k = 1:M
        c(k) = exp((1i)*(2*pi/lambda)*cos(angle)*d*(k-1)); 
    end
    
    % z (1x1) = c(1xM) * X(Mx1) * g(1x1);
    z = c*X*g;
    
    if (abs(z) > max_z)
        max_z = abs(z);
        % store the optimal theta;
        parameter(1,2) = theta;
    end    
end

% compute the amplitude;

c = zeros(1,M);

for k = 1:M
    c(k) = exp((1i)*(2*pi/lambda)*cos(parameter(1,2))*d*(k-1)); 
end

g = exp(-(1i)*(2*pi)*parameter(1,3)*freq);
z = c*X*g;

% store the optimal alpha;
parameter(1,1) = (1/(M*N))*z;


end
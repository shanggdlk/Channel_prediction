function parameter = MLE (X, freq, M, lambda,d)

N = 10;
%N is the number of samples in the freq domain;
parameter = zeros(1,3);
max_z = 0;

for theta = 0:pi/360:pi
    
    c = zeros(1,M);
    for k = 1:M
        c(k) = exp((1i)*(2*pi/lambda)*cos(theta)*d*(k-1)); 
    end
    
    for delay = 1:50
        dnm = delay/(1e-9);
        g = exp(-(1i)*(2*pi)*dnm*freq);
        z = c*X*g;
        if (abs(z) > max_z)
            max_z = abs(z);
            parameter(2) = theta;
            parameter(3) = dnm;
        end
    end
end

parameter(1) = (1/(M*N))*z;


end
function value = Zernike(n,m,r,phi)
% m = +1 for y-coma, m = -1 for x-coma
    
    value = 0;
    if mod(n-m,2) == 0
        for k = 0:(n-abs(m))/2
            value = value + (-1)^k * factorial(n-k) / factorial(k) / factorial((n+abs(m))/2-k) / factorial((n-abs(m))/2-k) * r^(n-2*k);
        end
        
        if m >= 0
            value = value * cos(m*phi);
        else
            value = value * sin(abs(m)*phi);
        end
    end

end
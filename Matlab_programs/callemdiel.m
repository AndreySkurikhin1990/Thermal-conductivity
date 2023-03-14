%       Call Emdiel
repeat = true;
while repeat == true
    n = input( 'input n   \n');
    
    if (n<=1.e0)
        fprintf('*** n must be >1!! ***\n');
        repeat = true;
    else
        emnml = 1.e0 - ((n-1)/(n+1))^2;
        fprintf([' n = %f, \n normal eps = %f,',  ...
                '\n hemispherical eps = %f \n'], ...
                 n, emnml, emdiel(n));
        repeat = false;
    end
end



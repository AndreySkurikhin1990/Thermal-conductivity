%       call emmet

fprintf('input n and k   \n');
Inputs = input('Example: [1.1, 1.3]   ');
n = Inputs(1);
k = Inputs(2);

emnml = 1.e0-((n-1)*(n-1)+k*k)/((n+1)*(n+1)+k*k);
fprintf( [' n = %f, k = %f, \n normnal eps = %f,', ...
          '\n hemispherical eps = %f \n'], ...
          n, k, emnml, emmet(n,k) );
      

function inteexp = integroexpon(n,x)
N=250;
E=0;
s=0;
if (x==0)
    E=1/(n-1);
else
for m=0:N 
    if (m==(n-1))
        continue; 
    end
        t=1;
        for k=1:m
            t=t*k;
        end
        s=s+(-x)^m/((m-n+1)*t);
end
e=0.577215665057043;
t=1;
for k=1:(n-1)
    t=t*k;
end;
E=(-x)^(n-1);
E=E/t;
t=0;
for k=1:(n-1)
    t=t+1/k;
end
psi=-e+t;
E=E*(-log(x)+psi);
E=E-s;
end
inteexp = PrAd(E);
end

function mk = PrAd(m)
if (isnan(m))
    m=0; 
end; 
if (isinf(m)) 
    m=0;  
end; 
mk=abs(real(m));
end
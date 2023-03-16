function t = tmp38()
format long g;
t=PoisKorn()
t=0;
end

function pk = PoisKorn()
a=1e-5; b=1e5; ep=1e-6; Nit=1e3;
k=0; ffv=1e-4; ra=abs(a-b); %ffv = 0,01 % - такая доля частиц не попадает в диапазон фракции
while (ra>ep)
    c=(a+b)/2;
    fa=erfc(a)-ffv;
    fb=erfc(b)-ffv;
    fc=erfc(c)-ffv;
    if ((fc*fb)>0) 
        if ((fa*fc)<0) 
            b=c; 
        end
    end
    if ((fc*fa)>0) 
        if ((fb*fc)<0) 
            a=c; 
        end
    end
    k=k+1;
    if (k>Nit) 
        break; 
    end
    ra=abs(a-b);
end
pk=c;
end
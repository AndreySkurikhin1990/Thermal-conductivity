%определение показателя преломления по коэффициенту отражения
lam = [464,560,656,750];
arRef464=[44,51,47,45,50,59,69,71,72,76];
arRef560=[34,39,33,33,38,48,59,63,66,71];
arRef656=[23,27,25,24,26,34,47,53,58,64];
arRef750=[13,15,15,13,15,23,35,40,44,50];
n=length(arRef464);
arRef1000 = [arRef750(n),arRef656(n),arRef560(n),arRef464(n)];
arRef900 = [arRef750(n-1),arRef656(n-1),arRef560(n-1),arRef464(n-1)];
for j=1:length(arRef1000)
    %Ref = arRef1000(j)/1e2;
    Ref = arRef900(j)/1e2;
    eps=1e-6;
    ae=1+eps;
    be=1e2;
    nit=0;
    ra=abs(ae-be);
    while (ra>eps)
ce=(ae+be)/2;
ma=dependRefln(ae)-Ref;
mc=dependRefln(ce)-Ref;
mb=dependRefln(be)-Ref;
if (ma*mc>0)        
ae=ce;    
else
if (mc*mb>0)
be=ce;
else disp(ce);
end
end
        nit=nit+1;
        if (nit>1e2)
            break; 
        end
    end
    nc(j)=ce;
end
nc=nc
arRef900=arRef900'
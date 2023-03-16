function [ prib ] = nachPribTemN(T0,Th,l,y0,n,ktp,kap)
sigma=5.668e-8; 
y=0:l:y0;
ksim=y/y0;
p=length(ksim);
cons=4*(n^2)*sigma/(3*ktp*kap);
for k=1:p
    ksi=ksim(k);
prch(k)=(T0*(1-ksi)+Th*ksi)+cons*((T0^4)*(1-ksi)+(Th^4)*ksi);
end
tem=0; ep=1e-7;
for q=1:p
    ksi=ksim(q);
    ta=0;
    tb=1e4;
    dl=abs(ta-tb);
    tc=(ta+tb)/2;
    while (dl>ep)
        fa=(ta+cons*ta^4-prch(q));
        fb=(tb+cons*tb^4-prch(q));
        fc=(tc+cons*tc^4-prch(q));
        if (fb*fc<0)
            ta=tc;
        end
        if (fa*fc<0)
            tb=tc;
        end
        dl=abs(ta-tb);
        tc=(ta+tb)/2;
    end
    tem(q)=tc;
end
prib=tem;
end
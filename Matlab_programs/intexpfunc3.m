function [ integ3 ] = intexpfunc3(tau0,ep,taum,Nto)
intt3=0;
if (tau0==0) 
    r=-1;
    intt3(1)=1/2;
else
    r=1;
end
for k=1:Nto+1
taum(k)=tau0-r*taum(k);
end
q=1; 
for k=1:Nto+1 
    a=taum(k);  
    y=abs(exp(-a)/a); 
    if (y<ep)          
        tm=0; 
    else          
        tm=integroexpon(3,a); 
    end; 
    intt3(q)=tm; 
    q=q+1; 
end;%int4tau=0; int4tau(1)=1/3; q=1; for k=1:Nto+1      a=taum(k); y=abs(exp(-a)/a); if (y<ep)         tm=0; else         tm=(exp(-a)-a*int3taus(q))/3; end; int4tau(q)=tm; q=q+1;end; int4tau=0; int4tau(1)=1/3; q=1; for k=1:Nto+1      a=taum(k); y=abs(exp(-a)/a); if (y<ep)         tm=0; else tm=integroexpon(4,a); end; int4tau(q)=tm; q=q+1;end; %disp(int4tau);%int4tau0m=0; int4tau0m(Nto+1)=1/3; q=1; for k=1:Nto+1      a=taum(k); a=tau0-a; y=abs(exp(-a)/a); if (y<ep)         tm=0;  else tm=(exp(-a)-int3tau0mins(q)*a)/3; end; int4tau0m(q)=tm; q=q+1; end; %int4tau0m=0; int4tau0m(Nto+1)=1/3; q=1; for k=1:Nto+1  a=taum(k); a=tau0-a; y=abs(exp(-a)/a); if (y<ep)  tm=0;  else tm=integroexpon(4,a); end; int4tau0m(q)=tm; q=q+1; end; disp(int4tau0m);
integ3=intt3;
end
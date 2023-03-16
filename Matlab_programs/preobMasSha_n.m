function [ pre ] = preobMasSha_n()
format longg;
alp2=1e6*SredGrafSha();
dlv2=(dlvoSham1()+dlvoSham2())/2;
ledlv2=length(dlv2);
dlv1=RasshDiapDlinVoln();
ledlv1=length(dlv1);
dlv1=preobDlinVoln(ledlv1,ledlv2,dlv1,dlv2,1);
%dlv1=dlvoVer53101(); ledlv1=length(dlv1);
t=1;
for j=1:ledlv1
    a=dlv1(j);
    fa=-1;
for k=2:ledlv2
        if ((dlv2(k-1)>=a) && (dlv2(k)<a))
            a1a=dlv2(k-1);
            b1a=dlv2(k);
            f1a=alp2(k-1);
            f2a=alp2(k);
            fa=1;
            break;
        end;
end;
if (fa>0)
    tm(t)=f1a+(a-a1a)*(f2a-f1a)/(b1a-a1a);
    dvtm(t)=a;
    t=t+1;
end
end;
a1=dlv2(1);            
b1=dlv2(2);
f1=alp2(1);
f2=alp2(2);
tm(1)=alp2(1)+(f2-f1)*(dlv1(1)-a1)/(b1-a1);
n=ledlv1;
m=ledlv2;
a1=dlv2(m-1);            
b1=dlv2(m);
f1=alp2(m-1);
f2=alp2(m);
tm(n)=f2+(f2-f1)*(dlv1(n)-b1)/(b1-a1);
tm(n-1)=(tm(n)+tm(n-2))/2;
for k=1:n
    tm(k)=ProvAd(tm(k));
end
pre=tm;
end

function [ w ] = preobDlinVoln(ledlv1,ledlv2,dlv1,dlv2,v)
for k=1:ledlv1
    dlv1(k)=1e-2/dlv1(k);
end
for k=1:ledlv2
    dlv2(k)=1e-2/dlv2(k);
end
if (v==1)
    w=dlv1;
else
    w=dlv2;
end
end

function mk = ProvAd(m)
if (isnan(m))   
    m=0; 
end; 
if (isinf(m))  
    m=0;  
end; 
mk=real(abs(m));
end
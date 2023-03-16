function [ pre ] = preobMas()
alp2=1e6*SredGrafSha();
dlv1=dlvoVer53101();
dlv2=dlvoSham1();
ledlv1=length(dlv1);
ledlv2=length(dlv2);
%disp(length(alp2));
tm=0;
for k=1:length(dlv1)
    if (dlv2(k)==dlv1(1)) 
        break; end;
end
n=k; t=1; a=0; b=0; a1=0; b1=0; f1=0; f2=0;
for j=1:ledlv1-1
    a=dlv1(j);
    b=dlv1(j+1);
    for k=n-1:ledlv2-1
        if (dlv2(k)<a)
            if (dlv2(k-1)>a)
            a1=dlv2(k-1);
            b1=dlv2(k);
            f1=alp2(k-1);
            f2=alp2(k);
            break; 
            end; 
        end;
    end;
    tm(t)=f2+(a-b1)*(f1-f2)/(a1-b1);
    t=t+1;
    if (b1<b)
        tm(t)=f2+(b-b1)*(f1-f2)/(a1-b1);
    t=t+1;
    end;
end;
tm(1)=alp2(n); 
tm(2)=(tm(3)+tm(1))/2;
for k=1:ledlv1
    if (isnan(tm(k)))
        disp(k);
    end
end
pre=tm;
end
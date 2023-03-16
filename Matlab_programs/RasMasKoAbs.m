function [ rasmaskopo ] = RasMasKoAbs()
alsr=skleiv2Mas(soglMas(1e6*SredGraf(),1e6*SredGraSt()),1e6*SredGraf());
rasmaskopo=alsr;
end

function [ obma ] = skleiv2Mas(ar1,ar2)
q=1;
obm=0;
for k=1:(length(ar1)+length(ar2)-1)
    if (k<length(ar1))
        obm(k)=ar1(k);
    else
        obm(k)=ar2(q);
        q=q+1;
    end 
end
obma=obm;
end

function [ prked ] = privedkEdiPropus(ar)
arr=ar;
p=length(arr);
for k=1:p
    if (arr(k)>1)
        arr(k)=1;
    end
    if (arr(k)<0)
        arr(k)=0;
    end
end
prked=arr;
end

function [ dlv ] = soglMas(ar1,ar2)
dl=dlvoi(); 
dlvv=dlvoVer53101();
kon=dlvv(1);
p=length(dl); 
q=0;
arr=0;
for k=1:p  
    if (dl(k)>=kon)
        q=q+1;
    else
        break;
    end
end
t2=ar2(q);
t1=ar1(1)-t2;
ko=(ar1(2)-ar1(1))/(dlvv(2)-dlvv(1));
t1=ar1(1)+(dl(q)-dlvv(1))*ko;
t1=(t1-t2);
for k=1:q+p
    if (k<=q)
    arr(k)=ar2(k)+t1;
    else
        break;
    end
end
dlv = arr;
end

function [ sred ] = SredGraSt()
Tal=privedkEdiPropus(Trko2()); Tal2=privedkEdiPropus(Trko4());
mkbr=232.188; mv=1.041; tol=0.67; rokbr=2.75; rov=0.112; 
mkbr2=228.831; mv2=0.979; tol2=0.64; rov2=0.092; 
vkbr=mkbr/(1e3*rokbr); vv=mv/(1e3*rov); 
vkbr2=mkbr2/(1e3*rokbr); vv2=mv2/(1e3*rov2);
xv=(vv/(vv+vkbr))*tol*1e3; 
xv2=(vv2/(vv2+vkbr2))*tol2*1e3; 
nTal=length(Tal);
for k=1:nTal 
    Tal(k)=-log(Tal(k))/xv; 
    Tal2(k)=-log(Tal2(k))/xv2;
end
alsr1=0;
for k=1:nTal 
    alsr1(k)=Tal(k)+Tal2(k);
end
sred=alsr1/2;
end
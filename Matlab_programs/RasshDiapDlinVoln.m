function [ ndv ] = RasshDiapDlinVoln()
Sp=skleiv2Mas(dlinvolndoMas(),dlinvoln());
ndv=Sp;
end

function [ dlv ] = dlinvoln()
dl=0; 
dl=dlvoVer53101(); 
p=length(dl); 
for k=1:p  
    dl(k)=1e-2/dl(k);
end
dlv = dl;
end

function [ dlv ] = dlinvolndoMas()
dl=0; 
dl=dlvoi();
wl=dlvoVer53101();
kon=wl(1);
dli=0;
p=length(dl); 
q=1;
for k=1:p  
    if (dl(k)>kon)
    dli(q)=1e-2/dl(k);
    q=q+1;
    else 
        break;
    end
end
dlv = dli;
end

function [ obma ] = skleiv2Mas(ar1,ar2)
q=1;
obm=0;
p=(length(ar1)+length(ar2));
r=length(ar1);
for k=1:p
    if (k<=r)
        obm(k)=ar1(k);
    else
        obm(k)=ar2(q);
        q=q+1;
    end 
end
obma=obm;
end
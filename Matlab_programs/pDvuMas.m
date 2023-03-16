function [ tetadm ] = pDvuMas(tetan,Tk,Tz)
teta=0;
q=1; 
for j=1:Ntt
    for k=1:Nto 
        teta(j,k)=tetan(q);    
        q=q+1; 
    end; 
end;
teta(Ntt,Nto+1)=Tk/Tz;
for k=1:Ntt-1
    teta(k,Nto+1)=teta(k+1,1);
end;
tetadm=teta;
end
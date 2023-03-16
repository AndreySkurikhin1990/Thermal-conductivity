function [koor] = vychKoor(y0,y0l,l,Nto,Ntt)
koo=0; ko=0; teta=0; tem=0;
for k=1:Ntt 
    koo(k,1)=ko; ko=ko+y0l; end;
ko=y0l;
for k=1:Ntt-1 
    koo(k,Nto+1)=ko; ko=ko+y0l; end; 
ko=l;
for k=1:Ntt 
    for j=2:Nto 
        koo(k,j)=ko; ko=ko+l; end; 
    ko=koo(k,1); end; 
koo(Ntt,Nto+1)=y0; 
koor=koo;
end
function [] = checkDerivTemFie()
kn=0; fl=0; 
for k=1:Nto 
    if (tetam(k)<tetam(Nto-n0)) 
        fl=1; break; end; end;
for k=1:Nto+1 
    kn(k)=(tetam(k)-tetam(1))/(koom(k)-koom(1)); end;
for k=round(Nto/2)-1:Nto 
    if (kn(k+1)>kn(k)) 
        fl=1; break; end; end;
n1=Nto-1; 
if (fl==1) 
    koop=koom; tetap=0; 
    for g=1:Nto 
        if (g>Nto) 
            break; 
        else
            tetap(g)=(tetam(g+1)-tetam(g))/(koop(g+1)-koop(g)); 
        end; end; 
    tetapp=tetap;
    for g=1:Nto 
        if (g==Nto) 
            break; else 
            tetapp(g)=(tetap(g+1)-tetap(g))/(koop(g+1)-koop(g)); end; end;
mpp=0; 
for k=1:Nto 
    if (mpp>tetapp(k)) 
        mpp=tetapp(k); 
        n1=k; end; 
end; end;
end
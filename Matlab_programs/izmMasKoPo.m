function [ alfom ] = izmMasKoPo(teta,Tz)
alfm=0;
for j=1:Ntt
    Tsr=teta(j,1)*Tz;
    %alf=koefpoglPatch(koo(j,1),Tsr);
    %alf=1e6/sprko(Tsr);
    alf=knusreddv(alfs,npp,dl,Tz);
    alfm(j)=alf;
    Tsr=0;
end
alfom=alfm;
end
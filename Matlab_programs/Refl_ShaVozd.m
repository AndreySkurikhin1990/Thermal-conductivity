function kootr = Refl_ShaVozd()
te=2e2:1e2:12e2; te0=273.15;
for j=1:length(te)
    tem=te(j)
knuSha=0; knuSha=preobMasSha_n(); 
%knuSha=0; knuSha=preobMas(); 
dl=0; dl=dlivoln();
ronu=0; nsha=0; 
nsh=0; nsh=Kramers_n_Sha(); 
hiSha=0;
p=length(nsh);
for k=1:p  
    hiSha(k)=dl(k)*knuSha(k)/(4*pi);
    nsha(k)=sqrt(nsh(k)^2+hiSha(k)^2);
    ronu(k)=formulaDunkla(nsha(k));
    ronup(k)=1-emmet(nsh(k),hiSha(k));
    ronud(k)=1-emdiel(nsha(k));
end
inpopo=(Refsred(hiSha,nsh,dl,tem+te0) + knusreddvvsr(hiSha,nsh,dl,tem+te0)) / 2;
inpopr=(Refsred(nsh,nsh,dl,tem+te0) + knusreddvvsr(nsh,nsh,dl,tem+te0)) / 2;
if (inpopo<inpopr)
Refl1 = (Refsred(ronu,nsha,dl,tem+te0) + knusreddvvsr(ronu,nsha,dl,tem+te0)) / 2;
Refl2 = (Refsred(ronu,nsh,dl,tem+te0) + knusreddvvsr(ronu,nsh,dl,tem+te0)) / 2;
Refl = (Refl1 + Refl2) / 2;
Refl1 = (Refsred(ronud,nsha,dl,tem+te0) + knusreddvvsr(ronud,nsha,dl,tem+te0)) / 2;
Refl2 = (Refsred(ronud,nsh,dl,tem+te0) + knusreddvvsr(ronud,nsh,dl,tem+te0)) / 2;
Refl = Refl/2 + (Refl1 + Refl2) / 4
else
Refl1 = (Refsred(ronup,nsha,dl,tem+te0) + knusreddvvsr(ronup,nsha,dl,tem+te0)) / 2;
Refl2 = (Refsred(ronup,nsh,dl,tem+te0) + knusreddvvsr(ronup,nsh,dl,tem+te0)) / 2;
Refl = (Refl1 + Refl2) / 2
end
%Refl = Refsred(ronud,nsha,dl,tem)
%Refl=knusreddvvsr(ronud,nsha,dl,tem)
end
kootr=0;
end
function t = formulaDunkla(nsha)
    %ronun=abs(((nsha-1)^2+hiSha^2)/((nsha+1)^2+hiSha^2)); %нормальное падение
    ronu=(1/2)+((nsha-1)*(3*nsha+1))/(6*((nsha+1)^2));
    ronu=ronu-(2*(nsha^3)*(nsha^2+2*nsha-1))/((nsha^2+1)*(nsha^4-1)); 
    podln=abs((nsha-1)/(nsha+1));
    ronu=ronu+(8*(nsha^4)*((nsha^4)+1)*log(nsha))/((nsha^2+1)*((nsha^4-1)^2));
    ronu=ronu+(nsha^2)*((nsha^2-1)^2)*log(podln)/((nsha^2+1)^3); %формула Данкла
    t = ronu;
end
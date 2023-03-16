%определеяет коэффициент отражения на границе вермикулит - шамот
function kootr = Refl_Sha(tem)
knuSha=0; 
knuSha=preobMas_n(); 
alfs=0; alfs=1e6*SredGraf; 
npp=Kramers_n();
dl=dlivoln();
ronu=0; nsh=0; nsha=0; 
nsh=Kramers_n_Sha(); 
hi=0; hiSha=0; podln=0;
p=length(npp);
for k=1:p  
    hiSha(k)=dl(k)*knuSha(k)/(4*pi); 
    hi(k)=dl(k)*alfs(k)/(4*pi); 
    nsha(k)=sqrt((nsh(k)^2+hiSha(k)^2)/(npp(k)^2+hi(k)^2)); 
    nsha(k)=abs(nsha(k));
    ronu(k)=0.5+((nsha(k)-1)*(3*nsha(k)+1))/(6*(nsha(k)+1)^2);
    ronu(k)=ronu(k)-(2*(nsha(k)^3)*(nsha(k)^2+2*nsha(k)-1))/((nsha(k)^2+1)*(nsha(k)^4-1)); 
    podln=(nsha(k)-1)/(nsha(k)+1);
    podln=abs(podln);
    ronu(k)=ronu(k)+(8*(nsha(k)^4)*((nsha(k)^4)+1)*log(nsha(k)))/((nsha(k)^2+1)*((nsha(k)^4-1)^2));
    ronu(k)=ronu(k)+(nsha(k)^2)*((nsha(k)^2-1)^2)*log(podln)/((nsha(k)^2+1)^3);
end;
%su=0; 
%for j=1:p-1    
%    su=su+(dl(j+1)-dl(j))*(nsha(j)+nsha(j+1))/2; 
%end; 
%su=0; 
%for j=1:p-1 
    %su=su+(dl(j+1)-dl(j))*(ronu(j)+ronu(j+1))/2; 
%end; 
%Refl=su/(dl(p)-dl(1)); 
%Refl=abs(real(Refl)); 
Refl = Refsred(ronu,npp,dl,tem)
Refl=knusreddvvsr(ronu,npp,dl,tem);
%su=0; for j=1:p-1   
%su=su+(dl(j+1)-dl(j))*(npp(j)+npp(j+1))/2; end; 
%su=0; for j=1:p-1 su=su+(dl(j+1)-dl(j))*(alfs(j)+alfs(j+1))/2; end; 
%alf=su/(dl(p)-dl(1)); alf=real(alf);
kootr=Refl;
end
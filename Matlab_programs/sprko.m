function sredpoRosdloslab = sprko(T)
alfs=1e6*SredGraf; 
npp=0; npp=Kramers_n(); 
dl=dlivoln(); p=length(dl);
betar=0;
fz = funizlvtorrod(npp,dl,T);
for k=1:p-1
betar=betar+(fz(k+1)-fz(k))/alfs(k);
end
betar=(1/betar);
sredpoRosdloslab=1/(betar*1e-6);
end
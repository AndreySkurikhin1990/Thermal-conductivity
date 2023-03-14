function [ dlv ] = dlivoln()
%vl=299792458; 
format long g;
%dl=0; dl=dlvoVer53101(); 
dl=RasshDiapDlinVoln();
p=length(dl); 
fid = fopen('Dlina_volny.txt','w');
for k=1:p  
    %dl(k)=1e-2/dl(k); %fr(k)=vl/dl(k); 
    fprintf(fid,'%0.20f\n',dl(k));
end
fclose(fid);
dlv = dl;
end
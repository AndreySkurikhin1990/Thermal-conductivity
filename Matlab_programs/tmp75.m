function t = tmp75()
Te=(1e2:1e2:1e3)+273.15;
koef=KoefEffecKTP(salosha,vybsha,vystsha);
ektpvi=EffecKTPSha(koef,Te-te0);
ektpvi=perevodvK(ektpvi,Te);
t=0;
end

function [ klaefm ] = KoefEffecKTP(salosha,vybsha,vystsha)
if ((salosha>=28e-2) && (salosha<38e-2))
if (vybsha==0) 
        if (vystsha==0) 
            kektp(4)=-0.435e-9; 
            kektp(3)=0.685e-6; 
            kektp(2)=0.134e-3; 
            kektp(1)=0.725; 
            porsha=(20+24)*1e-2/2;
        elseif (vystsha==1)
            kektp(4)=-0.867e-9; 
            kektp(3)=1.77e-6; 
            kektp(2)=-0.523e-3; 
            kektp(1)=0.806; 
            porsha=(24+30)*1e-2/2;
        end
end %задание коэффициентов - шамот средней пористости
if (vybsha==1)
    kektp(4)=-0.397e-9; 
    kektp(3)=0.71e-6; 
    kektp(2)=0.011e-3; 
    kektp(1)=0.851; 
    porsha=(16+20)*1e-2/2;
end %уплотненный шамот
if (vybsha==2) 
    kektp(4)=-0.377e-9; 
    kektp(3)=0.918e-6; 
    kektp(2)=-0.338e-3; 
    kektp(1)=0.77; 
    porsha=(30+33)*1e-2/2; 
end %низкоплотный шамот
if (vybsha==3) 
    kektp(4)=0; 
    kektp(3)=-0.607e-6; 
    kektp(2)=1.14e-3; 
    kektp(1)=0.641; 
    porsha=(10+16)*1e-2/2;
end
end %повышенной плотности шамот
if ((salosha>=38e-2) && (salosha<=45e-2))
if (vybsha==0) 
        if (vystsha==0) 
            kektp(4)=-0.124e-9; 
            kektp(3)=0.215e-6; 
            kektp(2)=0.125e-3; 
            kektp(1)=1.01; 
            porsha=(20+24)*1e-2/2;   
        elseif (vystsha==1) 
            kektp(4)=-0.333e-9; 
            kektp(3)=0.805e-6; 
            kektp(2)=-0.289e-3; 
            kektp(1)=0.903; 
            porsha=(24+30)*1e-2/2;
        end
end %задание коэффициентов - шамот средней пористости
if (vybsha==1) 
    kektp(4)=0; 
    kektp(3)=-0.154e-6; 
    kektp(2)=0.369e-3; 
    kektp(1)=1.03; 
    porsha=(16+20)*1e-2/2;
end %уплотненный шамот
if (vybsha==2) 
    kektp(4)=-0.377e-9; 
    kektp(3)=0.918e-6; 
    kektp(2)=-0.338e-3; 
    kektp(1)=0.77; 
    porsha=(30+33)*1e-2/2;
end %низкоплотный шамот
if (vybsha==3) 
    kektp(4)=0; 
    kektp(3)=-0.141e-6; 
    kektp(2)=0.437e-3; 
    kektp(1)=1.32; 
    porsha=(10+16)*1e-2/2;
end
end %повышенной плотности шамот
klaefm=kektp';
end

function [ t ] = EffecKTPSha(koef,Te)
ktp=0;
for k=1:length(Te)
    s=0;
    for j=1:length(koef)
        s=s+koef(j)*(Te(k)^(j-1));
    end
    ktp(k)=s;
end
t=ktp';
end

function [ ktp ] = perevodvK(ektp,Te)
koef=koefPribSha(ektp,Te);
ktp=EffecKTPSha(koef,Te);
end
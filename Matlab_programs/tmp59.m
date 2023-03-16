function t = tmp59(salosha,tem,vybsha,vystsha)
%format short;
%laefm=ObrMatr();
t=EffecKTP(salosha,tem,vybsha,vystsha);
end

function laefm = EffecKTP(salosha,tem,vybsha,vystsha)
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
    ktp=0;
for k=1:length(kektp)
    ktp=ktp+kektp(k)*(tem^(k-1));
end
laefm=ktp;
end

function t = ObrMatr()
n=3;
%md1=randi(10,n,n)
md1=[1,6,5;4,3,1;2,5,3];
for k=1:n
    for j=1:n
        md2(k,j)=0;
        if (k==j)
            md2(k,j)=1;
        end
    end
end
md2=preobMatr(md1,md2,n);
t=0;
end

function [ t ] = preobMatr(md1,md2,n)
md0=md1;
disp('Pervoe');
for k=1:n-1
    perv=md1(k,k);
    if (perv==0)
    else
    for j=1:n
        md1(k,j)=md1(k,j)/perv;
        md2(k,j)=md2(k,j)/perv;
    end
    end
    perv=1;
    for j=k+1:n
        vtor=md1(j,k);
        if (vtor==0)
        else
        for p=1:n
            md1(j,p)=md1(j,p)*perv/vtor-md1(k,p);
            md2(j,p)=md2(j,p)*perv/vtor-md2(k,p);
        end
        end
    end
    md2=md2;
end
disp('Vtoroe');
k=n;
while (k>=2)
    perv=md1(k,k);
    if (perv==0)
    else
    j=n;
    while (j>=1)
        md1(k,j)=md1(k,j)/perv;
        md2(k,j)=md2(k,j)/perv;
        j=j-1;
    end
    end
    md1=md1;
    md2=md2;
    j=k-1;
    while (j>=1)
        vtor=md1(j,k);
        if (vtor==0)
        else
            p=n;
        while (p>=1)
            md1(j,p)=md1(j,p)-md1(k,p)*vtor;
            md2(j,p)=md2(j,p)-md2(k,p)*vtor;
            p=p-1;
        end
        end
        j=j-1;
    end
    md1=md1;
    md2=md2;
    k=k-1;
end
md=0;
for k=1:n
    for j=1:n
        s=0;
        for p=1:n
            s=s+md0(k,p)*md2(p,j);
        end
        md(k,j)=s;
    end
end
md=md
md=md0*md2
t=md;
end
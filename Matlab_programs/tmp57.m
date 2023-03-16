function [ t ] = tmp57()
%format shortg;
format long g;
%p1=(16+20)/2; p2=(30+33)/2; p3=(10+12)/2;
%m=poisrasprpor0();
%t=poisrasprpor1(p1);
%t=poisrasprpor2(p2);
%t=poisrasprpor22(p3);
%t=poisrasprpor12(p3);
%n=4; %Ï‡ÍË  ¬»: 0 - 400, 1 - 500, 2 - 600, 3 - 700, 4 - 800, 5 - 900, 6 - 1000
%m=poisrasprporkvi(n);
vvi=0;
m=rasPorpoRazitom(vvi);
t=0;
end
%≈ÒÎË ·ÓÎ¸¯Â
function [ p ] = rasppor2(p2,p0,legr01,prgr01,ras01)
legr21=legr01*p2/p0; prgr21=prgr01*p2/p0; n=length(prgr01);
ras21=ras01; q=1; xp=0;
for k=1:n
    xt=prgr01(k);
    for  j=1:length(prgr21)-1
    if ((legr21(j)<xt) && (prgr21(j)>xt))
        kot=ras21(j)/(prgr21(j)-legr21(j));
        ras22(q)=kot*(xt-legr21(j));
        if (xp>0)
            if (xp>legr21(j)) 
                if (xp<prgr21(j))
                dob=kot*(xt-xp);
                ras22(q)=dob;
                end
            else
                if (j>1)
                kop=ras21(j-1)/(prgr21(j-1)-legr21(j-1));
                ras22(q)=ras22(q)+(legr21(j)-xp)*kop;
                end
            end
        end
        q=q+1;
    end
    end
    xp=xt;
    kop=kot;
    srra(k)=(legr21(k)+prgr21(k))/2;
end
n=length(ras22);
s=0;
for k=1:n
    s=s+ras22(k)*srra(k);
end
s=1e-2*1e-6*s'
p=ras22';
end
%≈ÒÎË ·ÓÎ¸¯Â
function t = poisrasprpor2(p2)
rasp0=rasprPorpoRazm0();
p0=raps0(1);
q=1; n=length(rasp0);
for k=2:n
    raspr0(q)=abs(rasp0(k)-rasp0(k-1))*1e2/p0; q=q+1;
end
legr0=LegrRasPorpoRazm0(); 
n=length(legr0);
for k=1:n-1
    prgr0(k)=legr0(k+1);
end
prgr0(n)=350; s=0; prgr0=prgr0';
for k=1:n
    if (prgr0(k)<=1)
        p=k;
        s=s+raspr0(k);
    end
end
ras01(1)=s; legr01(1)=0; prgr01(1)=1; q=2; n=length(legr0);
for k=p:n
    r=prgr0(k)-legr0(k);
    if ((prgr0(k)>1) && (r<=1))
        ras01(q)=raspr0(k);
        prgr01(q)=prgr0(k);
        q=q+1;
    end
    if (r>1)
        m=raspr0(k)/r;
        for j=1:r
        ras01(q)=m;
        prgr01(q)=prgr01(q-1)+1;
        q=q+1;
        end
    end
end
n=length(prgr01);
for k=2:n
    legr01(k)=prgr01(k-1);
end
ras01=ras01';
legr01=legr01';
prgr01=prgr01';
ras22=rasppor2(p2,p0,legr01,prgr01,ras01);
ras22=ras22';
t=0;
end
%≈ÒÎË ÏÂÌ¸¯Â
function t = poisrasprpor1(p2)
rasp0=rasprPorpoRazm0();
p0=rasp0(1);
q=1; n=length(rasp0); raspr0(1)=0;
for k=2:length(rasp0)
    raspr0(k)=abs(rasp0(k)-rasp0(k-1))*1e2/p0;
end
raspr0=raspr0';
s=0; n=length(raspr0);
for k=1:n
    s=s+raspr0(k);
end
legr0=LegrRasPorpoRazm0(); 
n=length(legr0);
for k=1:n-1
    prgr0(k)=legr0(k+1);
end
prgr0(n)=13e1; s=0; prgr0=prgr0';
for k=1:n
    if (prgr0(k)<=1)
        p=k;
        s=s+raspr0(k);
    end
end
ras01(1)=s; legr01(1)=0; prgr01(1)=1; q=2; 
n=(length(legr0)+length(prgr0)+length(raspr0))/3;
for k=p:n
    r=prgr0(k)-legr0(k);
    if ((prgr0(k)>1) && (r<=1))
        ras01(q)=raspr0(k);
        prgr01(q)=prgr0(k);
        q=q+1;
    end
    if (r>1)
        m=raspr0(k)/r;
        for j=1:r
        ras01(q)=m;
        prgr01(q)=prgr01(q-1)+1;
        q=q+1;
        end
    end
end
n=length(prgr01);
for k=2:n
    legr01(k)=prgr01(k-1);
end
ras01=ras01';
legr01=legr01';
prgr01=prgr01';
ras11=rasppor1(p2,p0,legr01,prgr01,ras01);
ras11=ras11';
t=0;
end
%≈ÒÎË ÏÂÌ¸¯Â
function [ p ] = rasppor1(p1,p0,legr01,prgr01,ras01)
legr11=legr01*p1/p0; prgr11=prgr01*p1/p0; n=length(prgr01); legr11=legr11'; prgr11=prgr11';
ras11=ras01'; legrm=0; prgrm=0; srra21=0;
q=1; xp=0;
for k=1:n
    xt=prgr01(k);
    for  j=1:n
        if ((legr11(j)<xt) && (prgr11(j)>xt))
        kot=ras11(j)/(prgr11(j)-legr11(j));
        dob1=kot*(xt-legr11(j));
        ras21(q)=dob1;
                if (xp<legr11(j))
                    if (j>1)
                kop=ras11(j-1)/(prgr11(j-1)-legr11(j-1));
                dob2=(legr11(j)-xp)*kop;
                ras21(q)=ras21(q)+dob2;
                    end
                end
                %if ((xp<prgr11(q)) && (q<2) && (j<3))
                    %ras21(q)=ras01(q)+ras21(q);
                %end
                %if ((xp>legr11(j)) && (xp<prgr11(j)))
                %dob=kot*(xt-xp);
                %ras21(q)=dob;
                %end
        end
    end
    legrm(q)=q-1; prgrm(q)=q; srra21(q)=(legrm(q)+prgrm(q))/2;
    q=q+1;
    xp=xt;
end
ras21=ras21';
n=length(ras21);
s=0;
for k=1:n
    s=s+ras21(k);
end
ras21=ras21*1e2/s;
s=0;
for k=1:n
    s=s+ras21(k);
end
s=s';
s=0;
for k=1:n
    s=s+srra21(k)*ras21(k);
end
s=1e-2*1e-6*s'
p=ras21;
end

function t = poisrasprpor0()
rasp0=rasprPorpoRazm0();
p0=rasp0(1);
q=1; n=length(rasp0);
for k=2:length(rasp0)
    raspr0(q)=abs(rasp0(k)-rasp0(k-1))*1e2/p0; 
    q=q+1;
end
raspr0=raspr0';
s=0; n=length(raspr0);
for k=1:n
    s=s+raspr0(k);
end
legr0=LegrRasPorpoRazm0(); 
n=length(legr0);
for k=1:n-1
    prgr0(k)=legr0(k+1);
end
prgr0(n)=legr0(n); s=0; prgr0=prgr0';
for k=1:n
    if (prgr0(k)<=1)
        p=k;
        s=s+raspr0(k);
    end
end
ras01(1)=s; legr01(1)=0; prgr01(1)=1; q=2; n=length(legr0);
for k=p:n
    r=prgr0(k)-legr0(k);
    if ((prgr0(k)>1) && (r<=1))
        ras01(q)=raspr0(k);
        prgr01(q)=prgr0(k);
        q=q+1;
    end
    if (r>1)
        m=raspr0(k)/r;
        for j=1:r
        ras01(q)=m;
        prgr01(q)=prgr01(q-1)+1;
        q=q+1;
        end
    end
end
n=length(prgr01);
for k=2:n
    legr01(k)=prgr01(k-1);
end
legr01=legr01';
prgr01=prgr01';
n=length(ras01); s=0;
for k=1:n
    s=s+ras01(k);
end
ras01=ras01';
t=0;
end

function [ r ] = rasprPorpoRazm0()
r=[21.8,21.8,21.75,21.625,21.25,20.75,20.25,19.75,19,18.25,17.5,12,9.875,8.75,8.25,7.25,6.25,5,4.25,3.25,0.875,0.325,0,0]';
end

function [ legr ] = LegrRasPorpoRazm0()
legr=[0,0.05,0.1,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50]';
end
%ÿœƒ-41
function [ r ] = rasprPorpoRazm1()
t=[0,0,0,0.23,0.29,0.4,0.46,0.75,2.48,2.48,7.86,10.8,23.92,41.88,1.79,1.04,3.64,1.96];
r=fliplr(t);
end

function [ legr ] = LegrRasPorpoRazm1()
t=[120,110,100,90,80,70,60,50,40,30,20,15,10,5,3,1,0.1,0]; %‡ÁÏÂ ‚ ÏÍÏ
legr=fliplr(t);
end

function t = poisrasprpor21(p2)
raspr0=rasprPorpoRazm1();
p0=11.0144;
legr0=LegrRasPorpoRazm1(); 
n=length(legr0);
for k=1:n-1
    prgr0(k)=legr0(k+1);
end
prgr0(n)=350; s=0; prgr0=prgr0';
for k=1:n
    if (prgr0(k)<=1)
        p=k;
        s=s+raspr0(k);
    end
end
ras01(1)=s; legr01(1)=0; prgr01(1)=1; 
q=2; n=length(legr0);
for k=p:n
    r=prgr0(k)-legr0(k);
    if ((prgr0(k)>1) && (r<=1))
        ras01(q)=raspr0(k);
        prgr01(q)=prgr0(k);
        q=q+1;
    end
    if (r>1)
        m=raspr0(k)/r;
        for j=1:r
        ras01(q)=m;
        prgr01(q)=prgr01(q-1)+1;
        q=q+1;
        end
    end
end
n=length(prgr01);
for k=2:n
    legr01(k)=prgr01(k-1);
end
ras01=ras01';
legr01=legr01';
prgr01=prgr01';
ras22=rasppor2(p2,p0,legr01,prgr01,ras01);
ras22=ras22';
t=0;
end

function t = poisrasprpor11(p2)
raspr0=rasprPorpoRazm1();
p0=11.0144;
legr0=LegrRasPorpoRazm1(); 
n=length(legr0);
for k=1:n-1
    prgr0(k)=legr0(k+1);
end
prgr0(n)=35e1; s=0; 
for k=1:n
    if (prgr0(k)<=1)
        p=k;
        s=s+raspr0(k);
    end
end
ras01(1)=s; legr01(1)=0; prgr01(1)=1; q=2; 
n=(length(legr0)+length(prgr0)+length(raspr0))/3;
for k=p:n
    r=prgr0(k)-legr0(k);
    if ((prgr0(k)>1) && (r<=1))
        ras01(q)=raspr0(k);
        prgr01(q)=prgr0(k);
        q=q+1;
    end
    if (r>1)
        m=raspr0(k)/r;
        for j=1:r
        ras01(q)=m;
        prgr01(q)=prgr01(q-1)+1;
        q=q+1;
        end
    end
end
n=length(prgr01);
for k=2:n
    legr01(k)=prgr01(k-1);
end
ras01=ras01';
legr01=legr01';
prgr01=prgr01';
ras11=rasppor1(p2,p0,legr01,prgr01,ras01);
ras11=ras11';
t=0;
end

function t = poisrasprpor01()
raspr0=rasprPorpoRazm1();
p0=11.0144;
legr0=LegrRasPorpoRazm1(); 
n=length(legr0);
for k=1:n-1
    prgr0(k)=legr0(k+1);
end
prgr0(n)=legr0(n); s=0; prgr0=prgr0';
for k=1:n
    if (prgr0(k)<=1)
        p=k;
        s=s+raspr0(k);
    end
end
ras01(1)=s; legr01(1)=0; prgr01(1)=1; q=2; n=length(legr0);
for k=p:n
    r=prgr0(k)-legr0(k);
    if ((prgr0(k)>1) && (r<=1))
        ras01(q)=raspr0(k);
        prgr01(q)=prgr0(k);
        q=q+1;
    end
    if (r>1)
        m=raspr0(k)/r;
        for j=1:r
        ras01(q)=m;
        prgr01(q)=prgr01(q-1)+1;
        q=q+1;
        end
    end
end
n=length(prgr01);
for k=2:n
    legr01(k)=prgr01(k-1);
end
legr01=legr01';
prgr01=prgr01';
ras01=ras01';
t=0;
end
%ÿ¡-1 2-1
function [ r ] = rasprPorpoRazm2()
t=[0, 10.6294, 14.0694, 17.8574, 23.8485, 33.3956, 39.3094, 44.8367, 50.0161 ...
53.8813, 59.6405, 61.4827, 63.8741, 68.1047, 71.3586, 75.9268, 79.3882, 83.2566 ...
86.9122, 89.9145, 92.5579, 95.1927, 96.8355, 97.998, 98.8324, 99.3257, 99.6412 ...
99.9201, 99.9442, 99.9575, 100]';
r=fliplr(t);
end

function [ legr ] = LegrRasPorpoRazm2()
t=[340.838, 200.63933, 45.34257, 26.45502, 20.39776, 16.04155, 12.64329 ...
10.19442, 8.15041, 6.63919, 5.27988, 4.16287, 3.46842, 2.80983, 2.2606 ...
1.81613, 1.44981, 1.17497, 0.93713, 0.75399, 0.61365, 0.49093, 0.39032 ...
0.31608, 0.2547, 0.20469, 0.16692, 0.1359, 0.10807, 0.0862, 0.06979]';
legr=fliplr(t);
end

function t = poisrasprpor22(p2)
rasp0=rasprPorpoRazm2();
p0=25.2;
rasp0=rasp0*p0/1e2;
n=length(rasp0);
raspr0(n)=0;
for k=2:n
    raspr0(k-1)=abs(rasp0(k)-rasp0(k-1))*1e2/p0; 
end
legr0=LegrRasPorpoRazm2();
n=length(legr0);
for k=1:n-1
    prgr0(k)=legr0(k+1);
end
prgr0(n)=35e1; s=0; prgr0=prgr0';
for k=1:n
    if (prgr0(k)<=1)
        p=k;
        s=s+raspr0(k);
    end
end
ras01(1)=s; legr01(1)=0; prgr01(1)=1; 
q=2; n=length(legr0);
for k=p:n
    r=prgr0(k)-legr0(k);
    if ((prgr0(k)>1) && (r<=1))
        ras01(q)=raspr0(k);
        prgr01(q)=prgr0(k);
        q=q+1;
    end
    if (r>1)
        m=raspr0(k)/r;
        for j=1:r
        ras01(q)=m;
        prgr01(q)=prgr01(q-1)+1;
        q=q+1;
        end
    end
end
n=length(prgr01);
for k=2:n
    legr01(k)=prgr01(k-1);
end
ras01=ras01';
legr01=legr01';
prgr01=prgr01';
ras22=rasppor2(p2,p0,legr01,prgr01,ras01);
ras22=ras22';
t=0;
end

function t = poisrasprpor12(p2)
rasp0=rasprPorpoRazm2();
p0=25.2;
rasp0=rasp0*p0/1e2;
n=length(rasp0);
for k=2:n
    raspr0(k-1)=abs(rasp0(k)-rasp0(k-1))*1e2/p0; 
end
raspr0(n)=0;
legr0=LegrRasPorpoRazm2(); 
n=length(legr0);
for k=1:n-1
    prgr0(k)=legr0(k+1);
end
prgr0(n)=35e1; s=0; 
for k=1:n
    if (prgr0(k)<=1)
        p=k;
        s=s+raspr0(k);
    end
end
ras01(1)=s; legr01(1)=0; prgr01(1)=1; q=2; 
n=(length(legr0)+length(prgr0)+length(raspr0))/3;
for k=p:n
    r=prgr0(k)-legr0(k);
    if ((prgr0(k)>1) && (r<=1))
        ras01(q)=raspr0(k);
        prgr01(q)=prgr0(k);
        q=q+1;
    end
    if (r>1)
        m=raspr0(k)/r;
        for j=1:r
        ras01(q)=m;
        prgr01(q)=prgr01(q-1)+1;
        q=q+1;
        end
    end
end
n=length(prgr01);
for k=2:n
    legr01(k)=prgr01(k-1);
end
ras01=ras01';
legr01=legr01';
prgr01=prgr01';
ras11=rasppor1(p2,p0,legr01,prgr01,ras01);
ras11=ras11';
t=0;
end

function t = poisrasprpor02()
rasp0=rasprPorpoRazm2();
p0=25.2;
rasp0=rasp0*p0/1e2;
n=length(rasp0);
for k=2:n
    raspr0(k-1)=abs(rasp0(k)-rasp0(k-1))*1e2/p0; 
    q=q+1;
end
raspr0(n)=0;
legr0=LegrRasPorpoRazm2(); 
n=length(legr0);
for k=1:n-1
    prgr0(k)=legr0(k+1);
end
prgr0(n)=legr0(n); s=0; prgr0=prgr0';
for k=1:n
    if (prgr0(k)<=1)
        p=k;
        s=s+raspr0(k);
    end
end
ras01(1)=s; legr01(1)=0; prgr01(1)=1; q=2; n=length(legr0);
for k=p:n
    r=prgr0(k)-legr0(k);
    if ((prgr0(k)>1) && (r<=1))
        ras01(q)=raspr0(k);
        prgr01(q)=prgr0(k);
        q=q+1;
    end
    if (r>1)
        m=raspr0(k)/r;
        for j=1:r
        ras01(q)=m;
        prgr01(q)=prgr01(q-1)+1;
        q=q+1;
        end
    end
end
n=length(prgr01);
for k=2:n
    legr01(k)=prgr01(k-1);
end
legr01=legr01';
prgr01=prgr01';
ras01=ras01';
t=0;
end
%ÿ¬-1 1-1
function [ r ] = rasprPorpoRazm3()
t=[0, 11.5833, 16.7564, 23.6914, 33.213, 42.772, 49.4821, 55.7423, 60.9905, 65.1889 ...
70.437, 71.3838, 72.3753, 73.7401, 76.0907, 79.529, 82.1748, 85.175, 88.2587, 90.8563 ...
93.1275, 95.443, 96.8905, 97.874, 98.6652, 99.0294, 99.3704, 99.6077, 99.7959, 99.9616 ...
99.9757, 100]';
r=fliplr(t);
end

function [ legr ] = LegrRasPorpoRazm3()
t=[341.24803, 200.87079, 45.36211, 26.43548, 20.38867, 16.04315, 12.64066 ...
10.19261, 8.14948, 6.63862, 5.27717, 4.144, 3.41414, 2.76709, 2.27202, 1.82495...
1.44478, 1.18063, 0.93948, 0.75206, 0.61279, 0.49187, 0.39041, 0.31641, 0.25516 ...
0.20479, 0.16705, 0.13588, 0.10807, 0.08619, 0.06981, 0.05648]';
legr=fliplr(t);
end

function t = poisrasprpor23(p2)
raspr0=rasprPorpoRazm3();
p0=26.5;
raspr0=raspr0*p0/1e2;
n=length(raspr0);
for k=2:n
    raspr0(k-1)=abs(rasp0(k)-rasp0(k-1))*1e2/p0; 
end
raspr0(n)=0;
legr0=LegrRasPorpoRazm3(); 
n=length(legr0);
for k=1:n-1
    prgr0(k)=legr0(k+1);
end
prgr0(n)=35e1; s=0; prgr0=prgr0';
for k=1:n
    if (prgr0(k)<=1)
        p=k;
        s=s+raspr0(k);
    end
end
ras01(1)=s; legr01(1)=0; prgr01(1)=1; q=2; n=length(legr0);
for k=p:n
    r=prgr0(k)-legr0(k);
    if ((prgr0(k)>1) && (r<=1))
        ras01(q)=raspr0(k);
        prgr01(q)=prgr0(k);
        q=q+1;
    end
    if (r>1)
        m=raspr0(k)/r;
        for j=1:r
        ras01(q)=m;
        prgr01(q)=prgr01(q-1)+1;
        q=q+1;
        end
    end
end
n=length(prgr01);
for k=2:n
    legr01(k)=prgr01(k-1);
end
ras01=ras01';
legr01=legr01';
prgr01=prgr01';
ras22=rasppor2(p2,p0,legr01,prgr01,ras01);
ras22=ras22'
t=0;
end

function t = poisrasprpor13(p2)
rasp0=rasprPorpoRazm3();
p0=26.5;
rasp0=rasp0*p0/1e2;
n=length(rasp0);
for k=2:n
    raspr0(k-1)=abs(rasp0(k)-rasp0(k-1))*1e2/p0;
end
raspr0(n)=0;
legr0=LegrRasPorpoRazm3(); 
n=length(legr0);
for k=1:n-1
    prgr0(k)=legr0(k+1);
end
prgr0(n)=35e1; s=0; 
for k=1:n
    if (prgr0(k)<=1)
        p=k;
        s=s+raspr0(k);
    end
end
ras01(1)=s; legr01(1)=0; prgr01(1)=1; q=2; 
n=length(legr0); n=length(prgr0); n=length(raspr0);
for k=p:n
    r=prgr0(k)-legr0(k);
    if ((prgr0(k)>1) && (r<=1))
        ras01(q)=raspr0(k);
        prgr01(q)=prgr0(k);
        q=q+1;
    end
    if (r>1)
        m=raspr0(k)/r;
        for j=1:r
        ras01(q)=m;
        prgr01(q)=prgr01(q-1)+1;
        q=q+1;
        end
    end
end
n=length(prgr01);
for k=2:n
    legr01(k)=prgr01(k-1);
end
ras01=ras01';
legr01=legr01';
prgr01=prgr01';
ras11=rasppor1(p2,p0,legr01,prgr01,ras01);
ras11=ras11';
t=0;
end

function t = poisrasprpor03()
rasp0=rasprPorpoRazm3();
p0=26.5;
rasp0=rasp0*p0/1e2;
n=length(rasp0);
for k=2:n
    raspr0(k-1)=abs(rasp0(k)-rasp0(k-1))*1e2/p0;
end
raspr0(n)=0;
legr0=LegrRasPorpoRazm3(); 
n=length(legr0);
for k=1:n-1
    prgr0(k)=legr0(k+1);
end
prgr0(n)=legr0(n); s=0; prgr0=prgr0';
for k=1:n
    if (prgr0(k)<=1)
        p=k;
        s=s+raspr0(k);
    end
end
ras01(1)=s; legr01(1)=0; prgr01(1)=1; q=2; n=length(legr0);
for k=p:n
    r=prgr0(k)-legr0(k);
    if ((prgr0(k)>1) && (r<=1))
        ras01(q)=raspr0(k);
        prgr01(q)=prgr0(k);
        q=q+1;
    end
    if (r>1)
        m=raspr0(k)/r;
        for j=1:r
        ras01(q)=m;
        prgr01(q)=prgr01(q-1)+1;
        q=q+1;
        end
    end
end
n=length(prgr01);
for k=2:n
    legr01(k)=prgr01(k-1);
end
legr01=legr01';
prgr01=prgr01';
ras01=ras01'
t=0;
end
%ÿœƒ
function [ r ] = rasprPorpoRazm4()
t=[0, 10.1713, 18.4934, 26.538, 36.3394, 48.4526, 56.4047, 65.0966, 71.6618, 76.5625, 83.4975, 84.9653 ...
87.1196, 90.9463, 93.8334, 96.7929, 97.6675, 98.5212, 99.0279, 99.396, 99.7339, 99.8082, 99.9293, 99.9293, 100]';
r=fliplr(t);
end

function [ legr ] = LegrRasPorpoRazm4()
t=[343.22645, 201.84496, 45.34218, 26.43562, 20.39193, 16.04629, 12.64759, 10.19305, 8.14791, 6.6392, 5.27874 ... 
4.12243, 3.40431, 2.77495, 2.25745, 1.82966, 1.45675, 1.1781, 0.94135, 0.75467, 0.61214, 0.49046, 0.3903, 0.31611, 0.2549]';
legr=fliplr(t);
end

function t = poisrasprpor24(p2)
raspr0=rasprPorpoRazm4();
p0=16.5;
raspr0=raspr0*p0/1e2;
n=length(raspr0);
for k=2:n
    raspr0(k-1)=abs(rasp0(k)-rasp0(k-1))*1e2/p0; 
end
raspr0(n)=0;
legr0=LegrRasPorpoRazm4(); 
n=length(legr0);
for k=1:n-1
    prgr0(k)=legr0(k+1);
end
prgr0(n)=35e1; s=0; prgr0=prgr0';
for k=1:n
    if (prgr0(k)<=1)
        p=k;
        s=s+raspr0(k);
    end
end
ras01(1)=s; legr01(1)=0; prgr01(1)=1; q=2; n=length(legr0);
for k=p:n
    r=prgr0(k)-legr0(k);
    if ((prgr0(k)>1) && (r<=1))
        ras01(q)=raspr0(k);
        prgr01(q)=prgr0(k);
        q=q+1;
    end
    if (r>1)
        m=raspr0(k)/r;
        for j=1:r
        ras01(q)=m;
        prgr01(q)=prgr01(q-1)+1;
        q=q+1;
        end
    end
end
n=length(prgr01);
for k=2:n
    legr01(k)=prgr01(k-1);
end
ras01=ras01';
legr01=legr01';
prgr01=prgr01';
ras22=rasppor2(p2,p0,legr01,prgr01,ras01);
ras22=ras22'
t=0;
end

function t = poisrasprpor14(p2)
rasp0=rasprPorpoRazm4();
p0=16.5;
rasp0=rasp0*p0/1e2;
n=length(rasp0);
for k=2:n
    raspr0(k-1)=abs(rasp0(k)-rasp0(k-1))*1e2/p0; 
end
raspr0(n)=0;
s=0; n=length(raspr0);
for k=1:n
    s=s+raspr0(k);
end
legr0=LegrRasPorpoRazm4(); 
n=length(legr0);
for k=1:n-1
    prgr0(k)=legr0(k+1);
end
prgr0(n)=35e1; s=0; 
for k=1:n
    if (prgr0(k)<=1)
        p=k;
        s=s+raspr0(k);
    end
end
ras01(1)=s; legr01(1)=0; prgr01(1)=1; q=2; 
n=length(legr0); n=length(prgr0); n=length(raspr0);
for k=p:n
    r=prgr0(k)-legr0(k);
    if ((prgr0(k)>1) && (r<=1))
        ras01(q)=raspr0(k);
        prgr01(q)=prgr0(k);
        q=q+1;
    end
    if (r>1)
        m=raspr0(k)/r;
        for j=1:r
        ras01(q)=m;
        prgr01(q)=prgr01(q-1)+1;
        q=q+1;
        end
    end
end
n=length(prgr01);
for k=2:n
    legr01(k)=prgr01(k-1);
end
ras01=ras01';
legr01=legr01';
prgr01=prgr01';
ras11=rasppor1(p2,p0,legr01,prgr01,ras01);
ras11=ras11';
t=0;
end

function t = poisrasprpor04()
rasp0=rasprPorpoRazm4();
p0=16.5;
rasp0=rasp0*p0/1e2;
n=length(rasp0);
for k=2:n
    raspr0(k-1)=abs(rasp0(k)-rasp0(k-1))*1e2/p0; 
end
raspr0(n)=0;
s=0; n=length(raspr0);
for k=1:n
    s=s+raspr0(k);
end
legr0=LegrRasPorpoRazm4(); 
n=length(legr0);
for k=1:n-1
    prgr0(k)=legr0(k+1);
end
prgr0(n)=legr0(n); s=0; prgr0=prgr0';
for k=1:n
    if (prgr0(k)<=1)
        p=k;
        s=s+raspr0(k);
    end
end
ras01(1)=s; legr01(1)=0; prgr01(1)=1; q=2; n=length(legr0);
for k=p:n
    r=prgr0(k)-legr0(k);
    if ((prgr0(k)>1) && (r<=1))
        ras01(q)=raspr0(k);
        prgr01(q)=prgr0(k);
        q=q+1;
    end
    if (r>1)
        m=raspr0(k)/r;
        for j=1:r
        ras01(q)=m;
        prgr01(q)=prgr01(q-1)+1;
        q=q+1;
        end
    end
end
n=length(prgr01);
for k=2:n
    legr01(k)=prgr01(k-1);
end
legr01=legr01';
prgr01=prgr01';
ras01=ras01';
t=0;
end

function [ t ] = rasprPorpoRazmAbskvi400()
t=[12.936, 12.851, 12.667, 12.500, 12.000, 10.320, 10.000, 9.125 ...
8.000, 6.980, 6.000, 4.917, 4.000, 3.125, 2.000, 1.423, 1.115, 0.846 ...
0.615, 0.471, 0.392, 0.314, 0.275]; %t=fliplr(t);
end

function [ t ] = LegrRasPorpoRazmkvi400()
t=[0.053, 0.631, 0.774, 1.000, 1.442, 5.374, 6.610, 10.000 ...
13.092, 16.408, 20.347, 26.640, 32.407, 43.481, 59.785 ...
71.225, 84.243, 100.000, 128.792, 149.908, 193.070 ... 
213.634, 261.567]'; %t=fliplr(t);
end

function [ r ] = rasprPorpoRazmOtnkvi400()
r=[0.66, 1.43, 1.29, 3.87, 12.99, 2.47, 6.76, 8.70, 7.89, 7.57, 8.37, 7.09, 6.76, 8.70, 4.46, 2.38 ...
2.08, 1.78, 1.12, 0.61, 0.61, 0.30, 2.12]; %r=fliplr(r);
end

function [ t ] = rasprPorpoRazmAbskvi500()
t=[10.414, 10.339, 10.207, 10.000, 9.822, 9.663, 9.465, 9.212 ...
9.059, 8.870, 8.664, 8.458, 8.335, 8.190, 8.000, 7.695, 7.477 ...
7.172, 6.859, 6.547, 6.105, 6.000, 4.929, 4.469, 4.000, 3.716 ...
3.406, 2.916, 2.699, 2.408, 2.168, 2.000, 1.814, 1.616, 1.414 ...
1.201, 0.935, 0.755, 0.667, 0.561, 0.464, 0.360, 0.298, 0.286 ...
0.280, 0.273]; %t=fliplr(t);
end

function [ t ] = LegrRasPorpoRazmkvi500()
t=[0.055, 0.445, 0.550, 0.640, 0.710, 0.802, 1.000, 1.482, 1.944, 2.348 ...
2.861, 3.521, 3.863, 4.665, 4.945, 6.015, 6.700, 7.786, 8.873, 10.000 ...
11.080, 11.256, 11.981, 13.127, 14.429, 15.924, 18.065, 22.529, 25.618 ...
30.105, 33.878, 37.381, 41.736, 46.477, 53.062, 61.330, 71.807, 84.073 ...
91.316, 100.000, 117.438, 137.916, 162.240, 190.638, 223.378, 252.815]'; %t=fliplr(t);
end

function [ r ] = rasprPorpoRazmOtnkvi500()
r=[0.72, 1.27, 1.99, 1.71, 1.52, 1.90, 2.44, 1.46, 1.82, 1.97, 1.97, 1.19, 1.38, 1.83, 2.93 ...
2.10, 2.93, 3.00, 3.00, 4.25, 1.00, 10.28, 4.42, 4.50, 2.73, 2.97, 4.71, 2.08, 2.80, 2.30 ...
1.62, 1.78, 1.90, 1.94, 2.04, 2.56, 1.72, 0.85, 1.02, 0.93, 1.00, 0.60, 0.12, 0.05, 0.07, 2.62];
r=fliplr(r);
end

function [ t ] = poisrasprporkvi(n)
prmi=1;
switch (n)
        case 4 %400
        raspr0=rasprPorpoRazmAbskvi400(); 
        %p0=51.74; 
        prgr0=LegrRasPorpoRazmkvi400();
        case 5 %500
        raspr0=rasprPorpoRazmAbskvi500(); 
        %p0=52.07; 
        prgr0=LegrRasPorpoRazmkvi500();
        case 6 %600
        raspr0=rasprPorpoRazmAbskvi600(); 
        %p0=51.51; 
        prgr0=LegrRasPorpoRazmkvi600();
        case 7 %700
        raspr0=rasprPorpoRazmAbskvi700(); 
        %p0=51.51; 
        prgr0=LegrRasPorpoRazmkvi700();
end
t=obrabMaskvi(raspr0,prmi,prgr0);
end

function [ t ] = obrabMaskvi(raspr0,prmi,prgr0)
v=1; prk=3e2; n1=length(prgr0); 
for k=v:prk
    prgr01(k)=k;
    legr01(k)=k-1;
    srra(k)=(legr01(k)+prgr01(k))/2e0;
    ras01(k)=0;
    ras02(k)=0;
end
n=length(prgr01);
r=1;
for k=1:n
    if ((prgr0(k)>=prmi) && (r>0))
        p=k; r=-1; break;
    end
end
legr0(1)=0;
for k=2:n1
    legr0(k)=prgr0(k-1);
end
mpg0=max(prgr0);
mrp0=min(raspr0);
mpg01=max(prgr01);
r=(raspr0(p)-raspr0(p-1))/(prgr0(p)-prgr0(p-1)); %ÓÚ 0 ‰Ó 1 ÏÍÏ
r=raspr0(p-1)+r*(prmi-prgr0(p-1));
w0=raspr0(v); pr=w0;
ras01(v)=w0; v=v+1; %2
ras01(v)=r; pr=r; v=v+1; v0=v; %3
u=1; 
ras02(u)=(w0-r)*1e2/w0; u=u+1; %2
for k=v0-1:n %Ë‰ÂÏ ÔÓ prgr01
    f=1;  x=0;
        for q=p:n1 %Ë‰ÂÏ ÔÓ prgr0
        if ((legr0(q)<=prgr01(k)) && (f>0) && (prgr0(q)>=prgr01(k)))
                x=q; f=-1; break;
        end
        end
        if ((x>1) && (f<0))
                s=raspr0(x); r=raspr0(x-1); s=s-r; 
                r=prgr0(x); w=legr0(x); r=r-w;
                ko=s/r;
                r=prgr01(k); w=legr0(x); r=r-w;
                w=raspr0(x-1); s=w+ko*r; 
        else
            r=mpg01-mpg0;
            w=-mrp0;
            r=w/r;
            w=prgr01(k)-mpg0;
            s=mrp0+r*w;
        end
        w=pr; w=(w-s)*1e2/w0;
        ras01(v)=s; pr=s; v=v+1;
        ras02(u)=w; u=u+1;
end
ko=240; v=v; u=u;
n=length(ras01); n=length(ras02); qg=n;
%ras010=ras01
s=0; k0=1e-2; k00=1e-6;
srra=srra*k00;
ras02=ras02*k0;
for k=1:n
    s=s+ras02(k)*srra(k);
end
s=s'
%s=sum(ras02);
%ras02=ras02*1e2/s;
%ras01=ras01'
s=sum(ras02);
ras02=ras02*1e2/s;
s=sum(ras02);
ras02=ras02(1:15);
%prgr01=prgr01';
t=0;
end

function [ t ] = rasprPorpoRazmAbskvi600()
t=[8.585, 8.558, 8.400, 8.246,  8.173, 8.000, 7.849, 7.652, 7.410, 7.115, 6.938, 6.721 ...
6.477, 6.222, 6.000, 5.886, 5.549, 5.338, 5.014, 4.593, 4.283, 4.000, 3.715, 3.178 ...
2.822, 2.362, 2.233, 2.000, 1.747, 1.531, 1.302, 1.093, 0.883, 0.739, 0.617, 0.528 ...
0.456, 0.313, 0.306, 0.275, 0.269]; %t=fliplr(t);
end

function [ t ] = LegrRasPorpoRazmkvi600()
t=[0.061, 0.606, 0.784, 0.921, 1.000, 1.238, 1.648, 2.296, 3.209 ...
4.477, 5.288, 6.247, 7.352, 8.682, 9.633, 10.321, 11.711, 12.674 ...
14.260, 16.050, 17.383, 18.013, 18.738, 19.882, 21.687, 24.608 ...
25.803, 27.377, 31.436, 35.390, 41.447, 48.541, 56.848, 66.578 ...
77.972, 84.381, 100.000, 123.240, 149.563, 190.423, 252.405]'; %t=fliplr(t);
end

function [ t ] = rasprPorpoRazmAbskvi700()
t=[5.678, 5.626, 5.579, 5.385, 5.142, 4.865, 4.709, 4.615 ...
4.561, 4.493, 4.446, 4.361, 4.242, 4.074, 4.000, 3.883 ...
3.440, 3.149, 2.864, 2.550, 2.259, 2.000, 1.904, 1.741 ...
1.617, 1.486, 1.381, 1.292, 1.205, 1.143, 1.043, 0.957 ...
0.876, 0.807, 0.725, 0.663, 0.594, 0.544, 0.469, 0.406 ...
0.384, 0.347, 0.285, 0.173, 0.149]; %t=fliplr(t);
end

function [ t ] = LegrRasPorpoRazmkvi700()
t=[0.060, 1.000, 1.105, 1.452, 1.865, 2.394, 2.947, 3.479 ...
4.109, 4.852, 5.730, 6.766, 7.990, 9.435, 10.000, 11.081 ...
13.716, 15.441, 17.367, 19.548, 22.002, 24.340, 25.761 ...
27.923, 30.219, 32.702, 35.390, 38.299, 41.447, 44.854 ...
48.541, 52.531, 56.848, 61.521, 66.578, 72.050, 77.972 ...
84.381, 91.316, 100.000, 107.542, 117.537, 138.150, 176.041, 233.572]'; %t=fliplr(t);
end
%»“ŒÃ-440
function [ t ] = rasPorpoRazmitom440()
rapo=[0.57, 0.3, 0.6, 0.87, 1.35, 2.07, 3.72, 3.81, 5.38 ...
    7.6, 9.67, 10.87, 34.68, 10.78, 6.25, 1.47, 0.0];
s=sum(rapo); 
rapo=rapo*1e2*1e-2/s;
rapo=fliplr(rapo);
t=rapo;
end
%»“ŒÃ-620
function [ t ] = rasPorpoRazmitom620()
rapo=[0.26, 0.15, 0.26, 0.22, 0.81, 0.81, 1.88 ...
    3.95, 5.54, 7.35, 7.09, 9.01, 34.9, 13.59, 7.5, 2.14, 4.54];
s=sum(rapo); 
rapo=rapo*1e2*1e-2/s;
rapo=fliplr(rapo);
t=rapo;
end
%»“ŒÃ-860
function [ t ] = rasPorpoRazmitom860()
rapo=[0.4, 0.09, 0.44, 0.22, 0.66, 1.02, 1.33 ...
    2.66, 4.07, 10.71, 12.17, 11.29, 35.06, 11.24, 7.13, 1.51, 0.0];
s=sum(rapo); 
rapo=rapo*1e2*1e-2/s;
rapo=fliplr(rapo);
t=rapo;
end
%»“ŒÃ-1000
function [ t ] = rasPorpoRazmitom1000()
rapo=[0.23, 0.19, 0.04, 0.61, 0.23, 1.03, 0.8 ...
    2.47, 5.66, 10.87, 14.18, 12.5, 32.61, 11.59, 5.25, 1.75, 0.0];
s=sum(rapo); 
rapo=rapo*1e2*1e-2/s;
rapo=fliplr(rapo);
t=rapo;
end
function [ t ] = PravGranPoromitom()
rapo=[15e1, 95, 85, 75, 65, 55, 45, 35, 25, 15, 7.5, 4, 2, 0.75, 0.3, 0.055, 0.0055]; 
rapo=(1e-6)*rapo; 
rapo=fliplr(rapo);
t=rapo;
end
function [ t ] = LevGranPoromitom(pravgr, levgr)
n=length(pravgr);
levgr(1)=0; 
for k=2:n
    levgr(k)=pravgr(k-1); 
end
t=levgr;
end
function [ t ] = SreRazPoritom(rpr, srra)
srpo=0; 
n=length(srra);
for k=1:n
    if (rpr(k)>0) 
        srpo=srpo+rpr(k)*srra(k)/1e2; 
    end
end
t=srpo;
end
function [ t ] = rasPorpoRazitom(vvi)
rpn=1e0; dp=1e0; 
s=0; l=0; h=1e-6; p=h; qg=2e2;
prgr0=PravGranPoromitom(); 
legr0=LevGranPoromitom(prgr0);
rpr0=NapMasRaspitom(vvi);
n=length(rpr0);
raspr=rasrasporporaz0(prgr0, legr0, rpr0, rpn, 0, n, dp);
qg=length(raspr);
for k=1:qg
    srra(k)=(p+l)/2e0; 
    s=s+srra(k)*raspr(k); 
    legr(k)=l; 
    prgr(k)=p; 
    p=p+h; 
    l=l+h;
end
s=s'
t=0;
end
function [ r ] = NapMasRaspitom(vvi)
if (vvi==0) 
    r=rasPorpoRazmitom440(); %»“ŒÃ-440
elseif (vvi==1) 
    r=rasPorpoRazmitom620(); %»“ŒÃ-620
elseif (vvi==2) 
    r=rasPorpoRazmitom860(); %»“ŒÃ-860
elseif (vvi==3) 
    r=rasPorpoRazmitom1000(); %»“ŒÃ-1000
end
end
function [ t ] = rasrasporporaz0(prgr0, legr0, raspr0, rpn, vyb, n, de)
e=1e-1; ht=1e0; e1=ht+e; koef=1e6;
prgr00=prgr0*koef;
legr00(1)=0; 
for k=2:n
    legr00(k)=prgr00(k-1);
end
s=0; 
for k=1:n
    if (prgr00(k)<e1) 
        p=k; s=s+raspr0(k); %‡ÁÏÂ ÔÓ ‰Ó 1 ÏÍÏ
    end
end
k=1; rprm(k)=s; q=1; 
for k=p:n 
    r=prgr00(k)-legr00(k);
    if ((prgr00(k)>e1) && (r<e1))
		rprm(q)=raspr0(k);
		q=q+1;
    end
    if (r>e1)
        m=raspr0(k)/r;
		t=r; j=0; 
        while (t>e) 
            j=j+1; 
            t=t-ht;
        end
            jk=j;
        for j=1:jk
		rprm(q)=m;
		q=q+1; 
        end
    end
end
prgr01=rpn; n=q-1; qg=n; 
for k=1:qg
    prgrm(k)=prgr01;
if (k>1) 
    legrm(k)=prgrm(k-1);
else
    legrm(k)=0;
end
    ras01(k)=rprm(k);
    prgr01=prgr01+de;
end
%prgrm=prgrm
ras01=ras01
if (vyb==0) 
    t=ras01;
elseif (vyb==1) 
    t=prgrm; 
elseif (vyb==2) 
    t=legrm;
end
end
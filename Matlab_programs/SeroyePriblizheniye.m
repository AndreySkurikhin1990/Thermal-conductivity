function r = SeroyePriblizheniye()
vybm=3; %выбор метода расчета
vybv=2; %выбор вещества
%для вермикулита
        vyfr=1; %выбор фракции
        vybu=1; %выбор укладки
        vybs=1; %выбор состояния
        vybmi=1; %выбор метода измерений 
%для шамота
        vybsha=1; %выбор вида шамота
        vystsha=1; %выбор стандартного шамота
        salosha=39e-2;
%для ИТОМ
        vybi=0; %выбор ИТОМ
switch (vybm)
    case 1 %метод Келлета, случай 1
        switch (vybv)
            case 1 %вермикулит
                t=Kellet_cas_1(vyfr,vybs,vybmi,vybu);
            case 2 %шамот
                t=Kellet_cas_1_Sha(vybsha,vystsha,salosha);
            case 3 %ИТОМ
                t=Kellet_cas_1_itom(vybi);
        end
    case 2
        switch (vybv)
            case 1 %вермикулит
                t=Kellet_cas_2(vyfr,vybs,vybmi,vybu);
            case 2
                t=Kellet_cas_2_Sha(vybsha,vystsha,salosha);
            case 3
                t=Kellet_cas_2_itom(vybi);
        end
    case 3
        switch (vybv)
            case 1 %для вермиулита
                t=tmp66(vyfv,vysv,vyuv,vmivmf);
            case 2 %для шамота
                t=tmp71(salosha,vybsha,vystsha);
            case 3 %для ИТОМ
                t=tmp72(vybi);
        end
end
r=t;
end
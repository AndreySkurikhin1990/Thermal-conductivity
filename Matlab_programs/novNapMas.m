function por = novNapMas(vyve, vymave, vyfv, vysv, vpmf, vpkf)
    switch (vyve) 
        case (0) %если выбран шамот
            ko=1e-2; porsha0=21.8*ko; porsha1=11.014*ko; porsha2=25.2*ko; porsha3=26.5*ko;
            porsha4=11.5*ko; porsha5=16.5*ko; %0 - Роучка, 1 - ШПД-41, 2 - ШБ-1 2-1, 3 - ШВ-1 1-1, 4 - ШПД, 5 - ШКУ-32 3-1
        switch (vymave) 
            case (0)
                por=porsha0; 
            case (1) 
                por=porsha1; 
            case (2) 
                por=porsha2;
            case (3) 
                por=porsha3; 
            case (4) 
                por=porsha4; 
            case (5) 
                por=porsha5;
        end
        case (1) %если выбран вермикулит
            if ((vyfv==0) || (vyfv==2)) %фракция 2-0,7 мм
                    if ((vysv==0) || (vysv==1)) %исходный или после повторных измерений
                            ko=1e-2; por207=66.35*ko; por16035=83.97*ko;
                            if (vpmf==0) %выбор пористости малой фракции
                                por=por207; 
                            elseif (vpmf==1) 
                                por=por16035; %выбор пористости малой фракции фракция 1,6-0,35 мм
                            end
                    elseif (vysv==2)
                        ko=1e-2; poro16035=84.36*ko;
                        por=poro16035; %фракция 1,6-0,35 мм
                    end
            elseif (vyfv==1) %фракция 8-4 мм
                if ((vysv==0)  || (vysv==1)) %исходный или после повторных измерений
                    ko=1e-2; poris84=55.75*ko; porin84=81.53*ko;
                                    if (vpkf==0)
                                        por=poris84; 
                                    elseif (vpkf==1) 
                                        por=porin84; %выбор пористости крупной фракции
                                    end
                elseif (vysv==2) %после обжига
                        ko=1e-2; poro84=86.61*ko;
                        por=poro84; 
                end
            end
        case (2) %если выбран ИТОМ
            ko=1e-2/2; pori440=(8e1+82e0)*ko; pori620=(75e0+78e0)*ko; 
            pori860=(65e0+68e0)*ko; pori1000=(62e0+65e0)*ko;
switch (vymave) 
        case (0)
            por=pori440; 
        case (1) 
            por=pori620; 
        case (2) 
            por=pori860; 
        case (3) 
            por=pori1000;
end
        case (3) %если выбран КВИ
    ko=1e-2; pork400=51.74*ko; pork500=52.07*ko; pork600=51.51*ko; pork700=39.75*ko;
	pork800=40.85*ko; pork900=39.37*ko; pork1000=36.07*ko;
                            switch (vymave)
                                case (4) 
                                    por=pork400; 
                                case (5) 
                                    por=pork500; 
                                case (6) 
                                    por=pork600; 
                                case (7) 
                                    por=pork700; 
                                case (8) 
                                    por=pork800; 
                                case (9) 
                                    por=pork900;
                                case (10) 
                                    por=pork1000;
                            end
    end
end
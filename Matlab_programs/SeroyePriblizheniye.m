function r = SeroyePriblizheniye()
vybm=3; %����� ������ �������
vybv=2; %����� ��������
%��� �����������
        vyfr=1; %����� �������
        vybu=1; %����� �������
        vybs=1; %����� ���������
        vybmi=1; %����� ������ ��������� 
%��� ������
        vybsha=1; %����� ���� ������
        vystsha=1; %����� ������������ ������
        salosha=39e-2;
%��� ����
        vybi=0; %����� ����
switch (vybm)
    case 1 %����� �������, ������ 1
        switch (vybv)
            case 1 %����������
                t=Kellet_cas_1(vyfr,vybs,vybmi,vybu);
            case 2 %�����
                t=Kellet_cas_1_Sha(vybsha,vystsha,salosha);
            case 3 %����
                t=Kellet_cas_1_itom(vybi);
        end
    case 2
        switch (vybv)
            case 1 %����������
                t=Kellet_cas_2(vyfr,vybs,vybmi,vybu);
            case 2
                t=Kellet_cas_2_Sha(vybsha,vystsha,salosha);
            case 3
                t=Kellet_cas_2_itom(vybi);
        end
    case 3
        switch (vybv)
            case 1 %��� ����������
                t=tmp66(vyfv,vysv,vyuv,vmivmf);
            case 2 %��� ������
                t=tmp71(salosha,vybsha,vystsha);
            case 3 %��� ����
                t=tmp72(vybi);
        end
end
r=t;
end
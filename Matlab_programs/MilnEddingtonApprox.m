function [ qrez ] = MilnEddingtonApprox()
pl=plot(te,ktep,'-b');
set(pl,'LineWidth',2);
hold on; grid on;
xlabel({'�����������, �'});
ylabel({'����������� ������������ ����������������, ��/(�*�)'}); 
title({'������ ����������� ������������� ��� �� �����������'});
end
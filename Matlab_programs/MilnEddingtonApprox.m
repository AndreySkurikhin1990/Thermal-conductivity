function [ qrez ] = MilnEddingtonApprox()
pl=plot(te,ktep,'-b');
set(pl,'LineWidth',2);
hold on; grid on;
xlabel({'Температура, К'});
ylabel({'Коэффициент кондуктивной теплопроводности, Вт/(м*К)'}); 
title({'График зависимости кондуктивного КТП от температуры'});
end
load('100Liver_x.mat')
load('100Liver_y.mat')
plot (x, y(:,2))
hold on
load('75Liver_y.mat')
plot (x, y(:,2))
load('50Liver_y.mat')
plot (x, y(:,2))
load('25Liver_y.mat')
plot (x, y(:,2))
load('10Liver_y.mat')
plot (x, y(:,2))
ylabel('[Ethanol] (mmol)','FontSize',18);
xlabel('time (min)','FontSize',18);
leg = legend ('100%','75%','50%','25%','10%','FontSize',14);
title(leg,'Liver Function','FontSize',14)
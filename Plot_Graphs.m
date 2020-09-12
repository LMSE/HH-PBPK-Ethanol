%% Plotting EtOH
figure

subplot(3,4,1)
plot (x, y(:,1));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
title('Adipose');

subplot(3,4,2)
plot (x, y(:,2));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
title('Blood');

subplot(3,4,3)
plot (x, y(:,3));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
title('Brain');

subplot(3,4,4)
plot (x, y(:,4));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
title('S.Int');

subplot(3,4,5)
plot (x, y(:,5));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
title('L.Int');

subplot(3,4,6)
plot (x, y(:,6));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
title('Heart');

subplot(3,4,7)
plot (x, y(:,7));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
title('Kidney');

subplot(3,4,8)
plot (x, y(:,8));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
title('Liver');

subplot(3,4,9)
plot (x, y(:,9));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
title('Lung');

subplot(3,4,10)
plot (x, y(:,10));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
title('Muscle');

%we skip pancreas, skin,spleen because they dont matter much

subplot(3,4,11)
plot (x, y(:,14));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
title('Stomach');

subplot(3,4,12)
plot (x, y(:,15));
hold on
plot (x, y(:,16));
plot (x, y(:,17));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
legend('Stomach','Sint','Lint')
title('Lumens');


figure

%% Plotting ACALD

subplot(3,4,1)
plot (x, y(:,18));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
title('Adipose');

subplot(3,4,2)
plot (x, y(:,19));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
title('Blood');

subplot(3,4,3)
plot (x, y(:,20));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
title('Brain');

subplot(3,4,4)
plot (x, y(:,21));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
title('S.Int');

subplot(3,4,5)
plot (x, y(:,22));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
title('L.Int');

subplot(3,4,6)
plot (x, y(:,23));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
title('Heart');

subplot(3,4,7)
plot (x, y(:,24));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
title('Kidney');

subplot(3,4,8)
plot (x, y(:,25));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
title('Liver');

subplot(3,4,9)
plot (x, y(:,26));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
title('Lung');

subplot(3,4,10)
plot (x, y(:,27));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
title('Muscle');

%we skip pancreas, skin,spleen because they dont matter much

subplot(3,4,11)
plot (x, y(:,31));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
title('Stomach');

subplot(3,4,12)
plot (x, y(:,32));
hold on
plot (x, y(:,33));
plot (x, y(:,34));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
legend('Stomach','Sint','Lint')
title('Lumens');











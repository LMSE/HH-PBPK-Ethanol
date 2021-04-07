%% Plotting EtOH
figure (2)

subplot(3,4,1)
plot (x, y(:,1));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
title('Adipose');
hold on

subplot(3,4,2)
plot (x, y(:,2));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
title('Blood');
hold on

subplot(3,4,3)
plot (x, y(:,3));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
title('Brain');
hold on

subplot(3,4,4)
plot (x, y(:,4));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
title('S.Int');
hold on

subplot(3,4,5)
plot (x, y(:,5));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
title('L.Int');
hold on

subplot(3,4,6)
plot (x, y(:,6));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
title('Heart');
hold on

subplot(3,4,7)
plot (x, y(:,7));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
title('Kidney');
hold on

subplot(3,4,8)
plot (x, y(:,8));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
title('Liver');
hold on

subplot(3,4,9)
plot (x, y(:,9));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
title('Lung');
hold on

subplot(3,4,10)
plot (x, y(:,10));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
title('Muscle');
hold on

%we skip pancreas, skin,spleen because they dont matter much

subplot(3,4,11)
plot (x, y(:,14));
xlabel('time(min)');
ylabel('concentration (mmol/L)');
title('Stomach');
hold on

% subplot(3,4,12)
% plot (x, y(:,15));
% hold on
% plot (x, y(:,16));
% plot (x, y(:,17));
% xlabel('time(min)');
% ylabel('concentration (mmol/L)');
% legend('Stomach','Sint','Lint')
% title('Lumens');
% hold on


%% Plotting ACALD
figure (3)

subplot(3,4,1)
plot (x, 1000*y(:,18));
xlabel('time(min)');
ylabel('concentration (umol/L)');
title('Adipose');
hold on

subplot(3,4,2)
plot (x, 1000*y(:,19));
xlabel('time(min)');
ylabel('concentration (umol/L)');
title('Blood');
hold on

subplot(3,4,3)
plot (x, 1000*y(:,20));
xlabel('time(min)');
ylabel('concentration (umol/L)');
title('Brain');
hold on

subplot(3,4,4)
plot (x, 1000*y(:,21));
xlabel('time(min)');
ylabel('concentration (umol/L)');
title('S.Int');
hold on

subplot(3,4,5)
plot (x, 1000*y(:,22));
xlabel('time(min)');
ylabel('concentration (umol/L)');
title('L.Int');
hold on

subplot(3,4,6)
plot (x, 1000*y(:,23));
xlabel('time(min)');
ylabel('concentration (umol/L)');
title('Heart');
hold on

subplot(3,4,7)
plot (x, 1000*y(:,24));
xlabel('time(min)');
ylabel('concentration (umol/L)');
title('Kidney');
hold on

subplot(3,4,8)
plot (x, 1000*y(:,25));
xlabel('time(min)');
ylabel('concentration (umol/L)');
title('Liver');
hold on

subplot(3,4,9)
plot (x, 1000*y(:,26));
xlabel('time(min)');
ylabel('concentration (umol/L)');
title('Lung');
hold on

subplot(3,4,10)
plot (x, 1000*y(:,27));
xlabel('time(min)');
ylabel('concentration (umol/L)');
title('Muscle');
hold on

%we skip pancreas, skin,spleen because they dont matter much

subplot(3,4,11)
plot (x, 1000*y(:,31));
xlabel('time(min)');
ylabel('concentration (umol/L)');
title('Stomach');
hold on

% subplot(3,4,12)
% plot (x, 1000*y(:,32));
% hold on
% plot (x, 1000*y(:,33));
% plot (x, 1000*y(:,34));
% xlabel('time(min)');
% ylabel('concentration (umol/L)');
% legend('Stomach','Sint','Lint')
% title('Lumens');










%Umulis Data
umulisx = [0 15 30 45 60 75 90 105 120 135 150 165];
umulisetoh = [0 5.3 6.2 5.7 5 4.1 3.2 2.5 1.9 1.3 1 0.7];
umulisacald = [0 2.2 2.7 2.9 3 2.9 2.75 2.5 2.2 1.5 1 0.75];


subplot (2,1,1)
plot (umulisx, umulisetoh, 'r-');
legend ('WMB-PBPK Model','Jones et al. 1998','Umulis Model')

subplot (2,1,2)
plot (umulisx, umulisacald,'r-');
legend ('WMB-PBPK Model','Jones et al. 1998','Umulis Model')
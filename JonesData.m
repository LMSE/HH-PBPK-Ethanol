jonesx = [0 15 30 45 60 75 90 105 120 135 150 165];
jonesetoh = [0 5.8 6.1 5.8 5 4.1 3.2 2.6 1.9 1.1 0.6 0.5];
jonesacald = [0 3.2 2.8 2.65 2.65 2.65 2.5 2.5 2.2 1.65 1.5 1.2];
jonesacaldmax = [0 4.6 4.3 3.75 3.8 3.8 3.6 3.3 2.65 2 1.95 1.4];
joneserror = jonesacaldmax - jonesacald;

subplot (2,1,1)
plot (jonesx, jonesetoh, 'x');


subplot (2,1,2)
errorbar (jonesx, jonesacald, joneserror, 'x');

% x = categorical({'100%','75%','50%','25%','10%'});
% x = reordercats(x,{'100%','75%','50%','25%','10%'});

x = [100, 75, 50, 25, 10];


%[100 75 50 25 10];

%Values: liver, c, u, br, sw, AUCet, AUCacald
vals = [195 0.0982 3.7195 0.4040 136.05 515.9670 371.9698;
        146.8492 0.0434 4.9922 0.5422 182.596 830.40 342.6625;
        72.4555 0.1691 5.5421 0.6019 202.7092 1170.7 226.0908;
        18.6777 0.0180 5.7146 .6206 209.0191 1399.5 86.1976;
        3.0066 0.0769 5.7494 .6244 210.2926 1465.5 13.0762];
%V
volume = [1.7723 1.8235 .3050 .6625 11.1565 1 1];

vals = vals./volume;

baseline = vals (1,:);

percentchange = vals./baseline;

% %b = bar(x,percentchange(:,1:5));
% b = bar(x,vals(:,1:5));
% legend ('Liver','Catalase','Urine','Breath','Sweat','FontSize',14);
% ylabel('Total Ethanol Elimination (mmol)','FontSize',18);
% %ylabel('% change in rate','FontSize',18);
% xlabel('% Liver Function','FontSize',18);


b = plot(x,vals(:,1:5));
    set ( gca, 'xdir', 'reverse' )
legend ('Liver','Catalase','Urine','Breath','Sweat','FontSize',14);
ylabel('Total Ethanol Elimination (mmol)','FontSize',18);
%ylabel('% change in rate','FontSize',18);
xlabel('% Liver Function','FontSize',18);


figure
%c = bar(x, percentchange(:,6:7));
c = bar(x, vals(:,6:7));
legend ('AUC(etoh)','AUC(acald)','FontSize',14);
ylabel('% change in AUC','FontSize',18);
xlabel('% Liver Function','FontSize',18);
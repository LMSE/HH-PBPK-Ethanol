
vals = [5445 7931; 5200 7759; 4835 7452; 4305 6831]; 


x = categorical({'40%','20%','12.5%','5.1%'});
x = reordercats(x,{'40%','20%','12.5%','5.1%'});

%b = bar(x,vals/1000);
bar(x,vals/1000)
ylabel('AUC (mol*min/L)','FontSize',18);
xlabel('Ethanol %','FontSize',18);

%leg = legend ('120 minutes','60 minutes')

title(legend ('120 minutes','60 minutes'),'Dosing Frequency')

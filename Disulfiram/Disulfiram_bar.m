
vals = [513; 895; 1815; 3942; 4525]; 


x = categorical({'0mg/L','2mg/L','4mg/L','6mg/L','8mg/L'});


b = bar(x,vals/1000);

ylabel('AUC (umol*min/L)','FontSize',18);
xlabel('[Disulfiram] (mg/L)','FontSize',18);


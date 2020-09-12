
vals = [510; 8796;848; 1542;1398; 3550;2128]; 


x = categorical({'ALDH2.1','ALDH2.2','ALDH2.3','ALDH2.4','ALDH2.5','ALDH2.6','ALDH2.7'});


b = bar(x,vals/1000);

ylabel('AUC (mmol*min/L)','FontSize',18);
xlabel('ALDH2 isoform','FontSize',18);


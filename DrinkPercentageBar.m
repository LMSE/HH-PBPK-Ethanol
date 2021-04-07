figure

x = categorical({'5.1%','12.5','20%','40%'});
x = reordercats(x,{'5.1%','12.5','20%','40%'});


%Values: EtOH, ACALD
AUC = [5427 1540; 5921 1562; 6194 1574; 6361 1576];

b = bar(x,AUC);
legend ('EtOH','Acetaldehyde', 'FontSize', 14);

yyaxis left
% set(gca,'ycolor','b')
ylabel('Ethanol AUC (mMol*min)','FontSize',18, 'Color','b');

yyaxis right
yticks([])
ylabel('Acetaldehyde AUC (uMol*min)','FontSize',18, 'Color','r');
xlabel('Drink Percentage','FontSize',18);
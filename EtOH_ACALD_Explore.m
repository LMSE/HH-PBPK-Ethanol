%% Load and Update Model, initialize diet intake
tic
clear
clc

global male

load('Harvey_1_03c.mat');
male = addReaction(male,'Skin_EX_etoh(c)_[bc]', 'reactionFormula', 'Skin_etoh[c]  <=> etoh[bc]', 'printLevel',0);
male = addReaction(male,'Lung_EX_etoh(br)_[bc]', 'reactionFormula', 'etoh[br]  <=> etoh[bc]', 'printLevel',0);
male = addReaction(male,'EX_etoh[br]', 'reactionFormula', 'etoh[br] -> ', 'printLevel',0);

male = addCOBRAConstraints(male,{'Skin_EX_etoh(c)_[bc]'},0,'c',[1], 'dsense','L','ConstraintID','slack_Skin_EX_etoh(c)_[bc]');
male = addCOBRAConstraints(male,{'Lung_EX_etoh(br)_[bc]'},0,'c',[1], 'dsense','L','ConstraintID','slack_Lung_EX_etoh(br)_[bc]');

diet=3000; %mmol/d

%% EtOH 

%all metabolites (non-slack) that involve EtOH, found by using Excel search for EtOH containing metabolites
EtOH_mets = [1569 1776 3607 3777 3989,...
4212 4826 5046 5213 5456 5956 7791 8020,...
8197 8400 10155 11236 11456 13904 14607 17799,...
23422 24217 46381 48133 53914 53915 53916 70762,...
70763 73018 80659 80660 80661 80662 80663 80664,...
81347 81348 81349 100133 100134 100369 100370 103788,...
103789 103790 103791 103792 103793 104707 104708 104709,...
104710 118111 118112 118113 118114 118115 125678 125679,...
125680 125681 126015 126016 126017]';

rxns_from_mets =... %reactions that are associated with EtOH_mets
[1482	1884	4467	4664	1607	1884	4467	4841	4664	5010	,...				
4841	5563	58551	5010	10053					,...
10617	14632	18452	24014	32802	61125	62318,...
14711	57891	58551	14344	57891	63555			,...	
5563	56999	5920	9099	9508	9666			,...
9099	10617	9508	9909	9666	10053				,...	
9909	14711	13723	14058	13723	14632	14058	14344	,...				
18452	18976	60357	18976	19500	24014	62318				,...	
31521	31522	32552	32552	32802	60357	61125	63555	63922	,...				
77223	77225	77223	77224	77225	30871	77223,...
30871	77223   30871	77225 3664	4664 3664	4664 3664	4467 3664	4467,...
1884	3664 1884	3664 3664	5010 3664	4841,...
3664	4841 32552	77386 32552	77386 32802	77386,...
32802	77386 9099	77381 9099	77381 9508	77381,...
9508	77381 9666	77381 9666	77381 10617	77381,...
10617	77381 9909	77381 10053	77381 17276	60357,...
17276	18452 17276	18452 17276	18976 17276	61125,...
13521	13723 13521	13723 13521	14058 13521	14058,...
13521	14632 13521	14632 13521	14344 81095 81096 81097]';

EtOH_rxns = unique(rxns_from_mets); %gives unique rxn IDs, lowest to high

%Relax EtOH constraints for further modification

%outputs reaction IDs related to EtOH, used 
% for i = 1:length(EtOH_mets)
%     male.mets(EtOH_mets(i));
%     findRxnIDs(male,findRxnsFromMets(male, male.mets(EtOH_mets(i)))) %outputs the rxn IDs correlated to each met
% end

% fix for transport rates
male.C(find(male.C == 20000)) = 1000000;
male.C(find(male.C == -20000)) = -1000000;

male = changeObjective(male, 'Liver_ALCD2if',1); %maximize ADH elimination
%male = changeObjective(male,'Colon_EX_etoh[luC]_[luLI]');



male.ub(strmatch('Diet_EX_etoh[d]',male.rxns))=-diet;
male.lb(strmatch('Diet_EX_etoh[d]',male.rxns))=-diet;

%exchange upper & lower bounds
% Liver can't pass 7200 mMol/day (actual MM limit is 6480mMol/day [this is 
% approx 20 standard drinks]), anything over that amount ends up in the 
% feces. This is a limitation of the model that can be fixed in the future. 

male.lb(strmatch('Pancreas_EX_etoh[luP]_[lu]',male.rxns))=0;

male.ub(strmatch('EX_etoh[br]',male.rxns))=.005*diet; %breath
male.lb(strmatch('EX_etoh[br]',male.rxns))=0.01;

male.lb(strmatch('Kidney_EX_etoh(e)_[bc]',male.rxns))=-.1*diet; %have to change from 14.94 otherwise limits urine output
male.ub(strmatch('EX_etoh[u]',male.rxns))=.1*diet;
male.lb(strmatch('EX_etoh[u]',male.rxns))=.03*diet;

male.ub(strmatch('EX_etoh[sw]',male.rxns))=.1*diet;
male.lb(strmatch('EX_etoh[sw]',male.rxns))=.03*diet;

%metabolism upper & lower bounds
male.ub(strmatch('Colon_CAT2p',male.rxns))=.02*diet;
male.lb(strmatch('Colon_CAT2p',male.rxns))=0;
% 
% male.ub(strmatch('Liver_ALCD2if',male.rxns))=.95*diet; %no need because FBA solves this
% male.lb(strmatch('Liver_ALCD2if',male.rxns))=.90*diet;
% 
male.ub(strmatch('Colon_ALCD2if',male.rxns))=0;
male.lb(strmatch('Colon_ALCD2if',male.rxns))=0;

male.ub(strmatch('Adipocytes_ALCD2if',male.rxns))=0;
male.lb(strmatch('Adipocytes_ALCD2if',male.rxns))=0;

male.ub(strmatch('Adipocytes_ALCD2yf',male.rxns))=0;
male.lb(strmatch('Adipocytes_ALCD2yf',male.rxns))=0;

%% ACALD

ACALD_mets = [1391 1392 1393  3540 3909 4827 4987 5150 5387 5791 5792 7792 8598];
% 1570 2579 5957 6918 7145 are not acetaldehyde

% To find reactions associated with each metabolite, listed in order
% for i = 1:length(ACALD_mets)
%     male.mets(ACALD_mets(i))
%     findRxnIDs(male,findRxnsFromMets(male, male.mets(ACALD_mets(i)))) %outputs the rxn IDs correlated to each met
% end

ACALD_rxns_all = [
    1393 1485 ...
    1393 1394 1482 1486 1854 3896 4416 80719 ...
    1394 1607 ...
    4416 4794 ...
    4794 5494 58421 ...
    10618 18453 24015 27019 32803 34199 39992 43159 48388 50162 52468 53545 54747 55819 56659 60934 69647 76835 ...
    5220 58421 ...
    5220 14287 72324 ...
    5494 56857 ...
    5836 5920 5922 6186 7035 9045 9100 80978 ...
    5836 6644 ...
    9100 10618 ...
    10754 10905 11411 11806];

ACALD_rxns = unique(ACALD_rxns_all);

% Getting rid of ACALD reactions not in PBPK model

male.ub(strmatch('Gall_ALDD2y',male.rxns))=0;
male.lb(strmatch('Gall_ALDD2y',male.rxns))=0;

male.ub(strmatch('Gall_DRPA',male.rxns))=0;
male.lb(strmatch('Gall_DRPA',male.rxns))=0;

male.ub(strmatch('Gall_r0186',male.rxns))=0;
male.lb(strmatch('Gall_r0186',male.rxns))=0;

male.ub(strmatch('Gall_ALDD2x',male.rxns))=0;
male.lb(strmatch('Gall_ALDD2x',male.rxns))=0;

male.ub(strmatch('Colon_ALCD2if',male.rxns))=0;
male.lb(strmatch('Colon_ALCD2if',male.rxns))=0;

male.ub(strmatch('Retina_EX_acald(e)_[bc]',male.rxns))=0;
male.lb(strmatch('Retina_EX_acald(e)_[bc]',male.rxns))=0;

male.ub(strmatch('Agland_EX_acald(e)_[bc]',male.rxns))=0;
male.lb(strmatch('Agland_EX_acald(e)_[bc]',male.rxns))=0;

male.ub(strmatch('Thyroidgland_EX_acald(e)_[bc]',male.rxns))=0;
male.lb(strmatch('Thyroidgland_EX_acald(e)_[bc]',male.rxns))=0;

male.ub(strmatch('Testis_EX_acald(e)_[bc]',male.rxns))=0;
male.lb(strmatch('Testis_EX_acald(e)_[bc]',male.rxns))=0;

male.ub(strmatch('Prostate_EX_acald(e)_[bc]',male.rxns))=0;
male.lb(strmatch('Prostate_EX_acald(e)_[bc]',male.rxns))=0;

male.ub(strmatch('CD4Tcells_EX_acald(e)_[bc]',male.rxns))=0;
male.lb(strmatch('CD4Tcells_EX_acald(e)_[bc]',male.rxns))=0;

male.ub(strmatch('Nkcells_EX_acald(e)_[bc]',male.rxns))=0;
male.lb(strmatch('Nkcells_EX_acald(e)_[bc]',male.rxns))=0;

male.ub(strmatch('Monocyte_EX_acald(e)_[bc]',male.rxns))=0;
male.lb(strmatch('Monocyte_EX_acald(e)_[bc]',male.rxns))=0;

male.ub(strmatch('Platelet_EX_acald(e)_[bc]',male.rxns))=0;
male.lb(strmatch('Platelet_EX_acald(e)_[bc]',male.rxns))=0;

male.ub(strmatch('RBC_EX_acald(e)_[bc]',male.rxns))=0;
male.lb(strmatch('RBC_EX_acald(e)_[bc]',male.rxns))=0;

male.ub(strmatch('Pancreas_EX_acald[luP]_[lu]',male.rxns))=0;
male.lb(strmatch('Pancreas_EX_acald[luP]_[lu]',male.rxns))=0;

male.ub(strmatch('GI_EX_acald[lu]_[d]',male.rxns))=0;
male.lb(strmatch('GI_EX_acald[lu]_[d]',male.rxns))=0;

    male.ub(strmatch('Liver_ALDD2x',male.rxns))=10;
    male.lb(strmatch('Liver_ALDD2x',male.rxns))=10;
        male.ub(strmatch('Liver_ALDD2y',male.rxns))=10;
    male.lb(strmatch('Liver_ALDD2y',male.rxns))=10;
%% Solving model and saving results

Solution=solveCobraLPCPLEX(male, 0, 0, 0, [], 0,'ILOGcomplex');

for i = 1:length(ACALD_rxns)
 ACALD.rxnID(i,1) = ACALD_rxns(i);
 ACALD.rxnNames(i,1)= male.rxnNames(ACALD_rxns(i))';
 ACALD.rxn(i,1)= male.rxns(ACALD_rxns(i))';
 ACALD.formulas(i,1)=printRxnFormula(male, male.rxns(ACALD_rxns(i)),false); %false turns off print
 ACALD.lb(i,1) = male.lb(ACALD_rxns(i))';
 ACALD.ub(i,1) = male.ub(ACALD_rxns(i))';
 ACALD.fluxes(i,1) = Solution.full((ACALD_rxns(i)))';
end

for ii = 1:length(EtOH_rxns)
 EtOH.rxnID(ii,1) = EtOH_rxns(ii);
 EtOH.rxnNames(ii,1)= male.rxnNames(EtOH_rxns(ii))';
 EtOH.rxn(ii,1)= male.rxns(EtOH_rxns(ii))';
 EtOH.formulas(ii,1)=printRxnFormula(male, male.rxns(EtOH_rxns(ii)),false); %false turns off print
 EtOH.lb(ii,1) = male.lb(EtOH_rxns(ii))';
 EtOH.ub(ii,1) = male.ub(EtOH_rxns(ii))';
 EtOH.fluxes(ii,1) = Solution.full((EtOH_rxns(ii)))';
end

toc
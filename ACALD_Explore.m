clear
clc

load ('harvey_1_03c.mat')

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

% Getting rid of reactions not in PBPK model

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


% SOlving model and saving results

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



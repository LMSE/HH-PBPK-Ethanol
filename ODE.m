function dC = ODE(c, age, sex, KEtOH,KACALD, v, q)

global kStom kSI kLI kStomSI kSILI kKid kPoop kLung aStom aLiver male bpmax iteration Solution rliv rlivALDH FBA aALDH 
global molliv molc molu molbr molsw

% ETOH Compartments
% 1.adipose 2.blood 3.brain 4.SI 5.LI 6.heart 7.kidney 8.liver 9.lung
% 10.muscle 11.pancreas 12.skin 13.spleen 14.stom 15.stom_lumen 
% 16.sint_lumen 17.lint_lumen
% ACALD Compartments 
% 18.adipose 19.blood 20.brain 21.SI 22.LI 23.heart 24.kidney 25.liver 26.lung
% 27.muscle 28.pancreas 29.skin 30.spleen 31.stom 32.stom_lumen 
% 33.sint_lumen 34.lint_lumen

r = maxrates(age, sex, c); %mM/min

% flow rate of blood into the liver: stom, s.int,l.int,pancreas,spleen,hep artery
qLiverIn=q(4)+q(5)+q(8)+q(11)+q(13)+q(14);

% concentration in vein going into liver: stom, s.int,l.int,pancreas,spleen,hep artery
cLiverInEtOH = (q(4)*c(4)/KEtOH(4) + q(5)*c(5)/KEtOH(5)+ q(8)*c(2) +q(11)*c(11)/KEtOH(11) + q(13)*c(13)/KEtOH(13) +q(14)*c(14)/KEtOH(14))/qLiverIn;
cLiverInACALD = (q(4)*c(21)/KACALD(4) + q(5)*c(22)/KACALD(5)+ q(8)*c(19) +q(11)*c(28)/KACALD(11) + q(13)*c(30)/KACALD(13) +q(14)*c(31)/KACALD(14))/qLiverIn;

% concentration in vein going back to lung with liver: adi,brain,heart,kid,liv,musc,skin   
        
cVenEtOH = (c(1)*q(1)/KEtOH(1) + c(3)*q(3)/KEtOH(3) + c(6)*q(6)/KEtOH(6) + c(7)*q(7)/KEtOH(7) + c(8)*qLiverIn/KEtOH(8) + c(10)*q(10)/KEtOH(10) + c(12)*q(12)/KEtOH(12))/q(9);
cVenACALD = (c(18)*q(1)/KACALD(1) + c(20)*q(3)/KACALD(3) + c(23)*q(6)/KACALD(6) + c(24)*q(7)/KACALD(7) + c(25)*qLiverIn/KACALD(8) + c(27)*q(10)/KACALD(10) + c(29)*q(12)/KACALD(12))/q(9);

dC=zeros(34,1); %initialize


%% Checks if we need to run FBA

global FBAcount ModelVer tol

MMMax = r(8)*60*24; %mmol/day, liver ethanol metabolism rate
MMMaxALDH = r(25)*60*24; %mmol/day, liver acald metabolism rate

if ModelVer==2 && ((iteration == 1 || abs(aLiver*MMMax - rliv)/(aLiver*MMMax) > tol) || abs(aLiver*MMMaxALDH - rlivALDH)/(aLiver*MMMaxALDH) > tol)
    %.1 --> 47 it @50min, up and down; .025 #114 @50 min, monotonic, no
    %woble; 0.05 #68@50min, no wobble   
    FBA = 1;
    FBAcount = FBAcount + 1;
    fprintf('etoh: %.2f (%.2f FBA) mM/day, acald: %.2f (%.2f FBA) mM/day\n',MMMax, rliv, MMMaxALDH, rlivALDH)
    fprintf('Running FBA again yields: \n')
else
    FBA = 0;
end

%% Running FBA

if (ModelVer ==2 && FBA == 1)
% Ethanol transport bounds
    male.ub(strmatch('Liver_ALCD2if',male.rxns))=aLiver*MMMax;
    male.lb(strmatch('Liver_ALCD2if',male.rxns))=aLiver*MMMax;
    %total metabolism, calculated based on liver metabolism rates assuming that the liver metabolizes 90-95% of all ethanol
        totalmetMax = MMMax/0.90;
        totalmetMin = MMMax/0.95;
    male.ub(strmatch('Diet_EX_etoh[d]',male.rxns))=0;
    male.lb(strmatch('Diet_EX_etoh[d]',male.rxns))=-100000;
    male.ub(strmatch('EX_etoh[br]',male.rxns))=.005*totalmetMax; %breath
    male.lb(strmatch('EX_etoh[br]',male.rxns))=.005*totalmetMax;
    male.lb(strmatch('Kidney_EX_etoh(e)_[bc]',male.rxns))=-.1*totalmetMax; %have to change from 14.94 otherwise limits urine output
    male.ub(strmatch('EX_etoh[u]',male.rxns))=.1*totalmetMax;
    male.lb(strmatch('EX_etoh[u]',male.rxns))=.03*totalmetMin;
    male.ub(strmatch('EX_etoh[sw]',male.rxns))=.1*totalmetMax;
    male.lb(strmatch('EX_etoh[sw]',male.rxns))=.03*totalmetMin;
% Ethanol metabolism upper & lower bounds
    male.ub(strmatch('Colon_CAT2p',male.rxns))=.02*totalmetMax;
    male.lb(strmatch('Colon_CAT2p',male.rxns))=0;
    male.ub(strmatch('Colon_ALCD2if',male.rxns))=0;
    male.lb(strmatch('Colon_ALCD2if',male.rxns))=0;
    male.ub(strmatch('Adipocytes_ALCD2if',male.rxns))=0;
    male.lb(strmatch('Adipocytes_ALCD2if',male.rxns))=0;
    male.ub(strmatch('Adipocytes_ALCD2yf',male.rxns))=0;
    male.lb(strmatch('Adipocytes_ALCD2yf',male.rxns))=0;
% ACALD transport bounds, lots of deletions because compartment not in PBPK model
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
    male.lb(strmatch('Pancreas_EX_acald[bc]_[luP]',male.rxns))=0;
    male.ub(strmatch('Pancreas_EX_acald[bc]_[luP]',male.rxns))=0;
    male.lb(strmatch('Liver_EX_acald[bpC]_[bpL]',male.rxns))=0;
    
    male.ub(strmatch('Liver_ALDD2x',male.rxns))=aLiver*MMMaxALDH;
    male.lb(strmatch('Liver_ALDD2x',male.rxns))=aLiver*MMMaxALDH;
    
% Solve FBA model
    Solution=solveCobraLPCPLEX(male, 0, 0, 0, [], 0,'ILOGcomplex');
    rliv=Solution.full(strmatch('Liver_ALCD2if',male.rxns)); %mmol/step
    rlivALDH=Solution.full(strmatch('Liver_ALDD2x',male.rxns));
    fprintf('etoh: %.2f (%.2f FBA) mM/day, acald: %.2f (%.2f FBA) mM/day\n',MMMax, rliv, MMMaxALDH, rlivALDH)
    fprintf('FBA #%.0f ran at iteration #%.0f (%.2f min) \n',FBAcount, iteration, iteration/10)
    fprintf('%.2f %.2f %.2f %.2f %.2f\n',molliv, molc, molu, molbr, molsw)
    fprintf('--------------------------------------------\n')
end

%% PBPK solver

if (ModelVer == 1) %(Base PK model)
%EtOH
    dC(1)=      (q(1)/v(1)) * (c(2) - c(1)/KEtOH(1));                          % 1.adipose
    dC(2)=      (q(9)/v(9)) * (c(9)/KEtOH(9) - c(2));                          % 2.artery
    dC(3)=      (q(3)/v(3)) * (c(2) - c(3)/KEtOH(3));                          % 3.brain   
    dC(4)=      (q(4)/v(4)) * (c(2) - c(4)/KEtOH(4))+ kSI*c(16);               % 4.s int
    dC(5)=      (q(5)/v(5)) * (c(2) - c(5)/KEtOH(5))+ kLI*c(17);               % 5.l int
    dC(6)=      (q(6)/v(6)) * (c(2) - c(6)/KEtOH(6));                          % 6.heart
    dC(7)=      (q(7)/v(7)) * (c(2) - c(7)/KEtOH(7));                          % kidney
    dC(8)=      (qLiverIn/v(8)) * (cLiverInEtOH-c(8)/KEtOH(8)) - aLiver*r(8);      % 8.liver, converted to mM/step
    dC(9)=      (q(9)/v(9)) * (cVenEtOH - c(9)/KEtOH(9));                          % 9.lung
    dC(10)=      (q(10)/v(10)) * (c(2) - c(10)/KEtOH(10));                     %10.muscle
    dC(11)=      (q(11)/v(11)) * (c(2) - c(11)/KEtOH(11));                     %11.pancreas
    dC(12)=      (q(12)/v(12)) * (c(2) - c(12)/KEtOH(12));                     %12.skin
    dC(13)=      (q(13)/v(13)) * (c(2) - c(13)/KEtOH(13));                     %13.spleen
    dC(14)=      (q(14)/v(14)) * (c(2) - c(14)/KEtOH(14)) +kStom*c(15) -aStom*r(14) ;%14.stomach
    dC(15)=      -kStom*c(15) - kStomSI*c(15);                             %15.stomach_lumen
    dC(16)=      -kSI*c(16) - kSILI*c(16)+kStomSI*c(15);                   %16.sint_lumen
    dC(17)=      -kLI*c(17) - kPoop*c(17)+kSILI*c(16);                     %17.lint_lumen
%ACALD
    dC(18)=      (q(1)/v(1)) * (c(19) - c(18)/KACALD(1));                         % 1.adipose
    dC(19)=      (q(9)/v(9)) * (c(26)/KACALD(9) - c(19));                         % 2.artery
    dC(20)=      (q(3)/v(3)) * (c(19) - c(20)/KACALD(3));                         % 3.brain   
    dC(21)=      (q(4)/v(4)) * (c(19) - c(21)/KACALD(4));              % 4.s int
    dC(22)=      (q(5)/v(5)) * (c(19) - c(22)/KACALD(5));              % 5.l int
    dC(23)=      (q(6)/v(6)) * (c(19) - c(23)/KACALD(6));                         % 6.heart
    dC(24)=      (q(7)/v(7)) * (c(19) - c(24)/KACALD(7));                         % kidney
    dC(25)=      (qLiverIn/v(8)) * (cLiverInACALD-c(25)/KACALD(8)) + aLiver*r(8)/360 - aALDH*r(25);     % 8.liver, converted to mM/step
    dC(26)=      (q(9)/v(9)) * (cVenACALD - c(26)/KACALD(9));                         % 9.lung
    dC(27)=      (q(10)/v(10)) * (c(19) - c(27)/KACALD(10));                     %10.muscle
    dC(28)=      (q(11)/v(11)) * (c(19) - c(28)/KACALD(11));                     %11.pancreas
    dC(29)=      (q(12)/v(12)) * (c(19) - c(29)/KACALD(12));                     %skin
    dC(30)=      (q(13)/v(13)) * (c(19) - c(30)/KACALD(13));                     %13.spleen
    dC(31)=      (q(14)/v(14)) * (c(19) - c(31)/KACALD(14));%14.stomach
    dC(32)=      0;                             %15.stomach_lumen
    dC(33)=      0;                   %16.sint_lumen
    dC(34)=      0;                     %17.lint_lumen
elseif (ModelVer == 2) %(GEM-PK model)
%EtOH
    dC(1)=      (q(1)/v(1)) * (c(2) - c(1)/KEtOH(1));                         % 1.adipose
    dC(2)=      (q(9)/v(9)) * (c(9)/KEtOH(9) - c(2));                         % 2.artery
    dC(3)=      (q(3)/v(3)) * (c(2) - c(3)/KEtOH(3));                         % 3.brain   
    dC(4)=      (q(4)/v(4)) * (c(2) - c(4)/KEtOH(4))+ kSI*c(16);              % 4.s int
    dC(5)=      (q(5)/v(5)) * (c(2) - c(5)/KEtOH(5))+ kLI*c(17) - Solution.full(strmatch('Colon_CAT2p',male.rxns))/(1440*20); 
    dC(6)=      (q(6)/v(6)) * (c(2) - c(6)/KEtOH(6));                         % 6.heart
    dC(7)=      (q(7)/v(7)) * (c(2) - c(7)/KEtOH(7)) - Solution.full(strmatch('EX_etoh[u]',male.rxns))/(1440*20);              % 7.kidney% 5.l int
    dC(8)=      (qLiverIn/v(8)) * (cLiverInEtOH-c(8)/KEtOH(8)) - aLiver*rliv/(24*60);  % 8.liver, converted to mM/step 
    dC(9)=      (q(9)/v(9)) * (cVenEtOH - c(9)/KEtOH(9)) - Solution.full(strmatch('EX_etoh[br]',male.rxns))/(1440*20);             % 9.lung
    dC(10)=      (q(10)/v(10)) * (c(2) - c(10)/KEtOH(10));                     %10.muscle
    dC(11)=      (q(11)/v(11)) * (c(2) - c(11)/KEtOH(11));                     %11.pancreas
    dC(12)=      (q(12)/v(12)) * (c(2) - c(12)/KEtOH(12)) - Solution.full(strmatch('EX_etoh[sw]',male.rxns))/(1440*20);    %12.skin
    dC(13)=      (q(13)/v(13)) * (c(2) - c(13)/KEtOH(13));                     %13.spleen
    dC(14)=      (q(14)/v(14)) * (c(2) - c(14)/KEtOH(14)) +kStom*c(15) -aStom*r(14);%14.stomach
    dC(15)=      -kStom*c(15) - kStomSI*c(15);                             %15.stomach_lumen
    dC(16)=      -kSI*c(16) - kSILI*c(16)+kStomSI*c(15);                   %16.sint_lumen
    dC(17)=      -kLI*c(17) - kPoop*c(17)+kSILI*c(16);                     %17.lint_lumen
%ACALD
    dC(18)=      (q(1)/v(1)) * (c(19) - c(18)/KACALD(1));                         % 1.adipose
    dC(19)=      (q(9)/v(9)) * (c(26)/KACALD(9) - c(19));                         % 2.artery
    dC(20)=      (q(3)/v(3)) * (c(19) - c(20)/KACALD(3));                         % 3.brain   
    dC(21)=      (q(4)/v(4)) * (c(19) - c(21)/KACALD(4));              % 4.s int
    dC(22)=      (q(5)/v(5)) * (c(19) - c(22)/KACALD(5)); 
    dC(23)=      (q(6)/v(6)) * (c(19) - c(23)/KACALD(6));                         % 6.heart
    dC(24)=      (q(7)/v(7)) * (c(19) - c(24)/KACALD(7))- Solution.full(strmatch('Colon_ALDD2xm',male.rxns))/(1440*20)-Solution.full(strmatch('Colon_ALDD2y',male.rxns))/(1440*20) - Solution.full(3896)/(1440*20);              % 7.kidney% 5.l int, 3896 because strmatch no work
    dC(25)=      (qLiverIn/v(8)) * (cLiverInACALD-c(25)/KACALD(8)) + aLiver*rliv/(1440*360)- aALDH*rlivALDH/(1440);  % 8.liver, converted to mM/step 
    %dC(25)=      (qLiverIn/v(8)) * (cLiverInACALD-c(25)/KACALD(8)) + %aLiver*rliv/(1440*360)- aALDH*r(25);  % 8.liver, converted to mM/step, works
    dC(26)=      (q(9)/v(9)) * (cVenACALD - c(26)/KACALD(9));             % 9.lung
    dC(27)=      (q(10)/v(10)) * (c(19) - c(27)/KACALD(10));                     %10.muscle
    dC(28)=      (q(11)/v(11)) * (c(19) - c(28)/KACALD(11));                     %11.pancreas
    dC(29)=      (q(12)/v(12)) * (c(19) - c(29)/KACALD(12));    %12.skin
    dC(30)=      (q(13)/v(13)) * (c(19) - c(30)/KACALD(13));                     %13.spleen
    dC(31)=      (q(14)/v(14)) * (c(19) - c(31)/KACALD(14));%14.stomach
    dC(32)=      0;                             %15.stomach_lumen
    dC(33)=      0;                   %16.sint_lumen
    dC(34)=      0;                     %17.lint_lumen
else
    error ('invalid model number')
end

% molliv = molliv + aLiver*rliv/(24*60*20)*v(8);
% molc = molc + Solution.full(strmatch('Colon_CAT2p',male.rxns))/(1440*20)*v(5);
% molu = molu + Solution.full(strmatch('EX_etoh[u]',male.rxns))/(1440*20)*v(7);
% molbr = molbr + Solution.full(strmatch('EX_etoh[br]',male.rxns))/(1440*20)*v(9);
% molsw = molsw + Solution.full(strmatch('EX_etoh[sw]',male.rxns))/(1440*20)*v(12);


end


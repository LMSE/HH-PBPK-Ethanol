w%% Initialize model
clear
clc
global male fat
global  kStom kSI kLI kStomSI kSILI kKid kPoop kLung aStom aLiver rliv male FBA ModelVer aALDH tol

% ETOH Compartments
% 1.adipose 2.blood 3.brain 4.SI 5.LI 6.heart 7.kidney 8.liver 9.lung
% 10.muscle 11.pancreas 12.skin 13.spleen 14.stom 15.stom_lumen 
% 16.sint_lumen 17.lint_lumen
% ACALD Compartments 
% 18.adipose 19.blood 20.brain 21.SI 22.LI 23.heart 24.kidney 25.liver 26.lung
% 27.muscle 28.pancreas 29.skin 30.spleen 31.stom 32.stom_lumen 
% 33.sint_lumen 34.lint_lumen

%% Inputs
ModelVer = 1; %1 = PBPK, 2 = GEM-PBPK
tol = .01;

tMax = 500; %minutes
drinkPercent = .4;
gkg=0.25; % grams/kg of EtOH drank, stdrinks = 14g/standard drink, gkg = sdrinks*14/kg

%multidose example
%multi = [20 14; 90 14; 95 14;100 14];
multi = [];
n = size(multi,1);

bodyMass = 74.5; age = 25.6;  sex = 0; %0=male, 1=female (not yet supported)
height = 180; %backcalculated from BMI of 26.35

fat = 20;       %fat percent by weight
kStom = 0.7135*drinkPercent^2-0.0985*drinkPercent + 0.0112;  %absorption rate by Stom (.013 or .065)
kStomSI = 1.953*drinkPercent^2-0.168*drinkPercent+0.0255;    %rate of movement from Stom_lumen to SI_lumen
kSI = -0.006*drinkPercent^2+0.0686*drinkPercent+0.0615;        %absorption rate by SI
kLI = 0.00;        %absorption rate by LI
kSILI = 0;      %rate of movement from SI_lumen to LI_lumen
kKid = 0.0;       %rate of elimination to urine
kPoop = 0;      %rate of elimination to feces, 0 because no ethanol in feces
kLung = 0.0;      %rate of elimination via lungs
aStom = 1;      %amount of reaction in Stomach (1 is 100%)
aLiver = 1.25;     %amount of reaction in Liver (1 is 100%)

ALDHtype = 1;   %1.WT 2.E504K (East Asian) 3.I41V (African) 4.P92T (Latino) 5. T244M(S Asian) 6.V304M (Latino) 7.R338W (Finnish) 
mglDisulfiram = 0; %mg/L

Disulfiram = 1000*mglDisulfiram/296.539; %uM of Disulfiram, ALDH inhibitor
aALDH = aALDHtype(ALDHtype)*(0.0015*Disulfiram^2 - 0.0734*Disulfiram + 1); %correlation based off Kitson et al, 1978 on Sheep liver
%% Initialize Physiological parameters
BSA=0.007184.*bodyMass.^.425.*height.^.725; %(m^2)
v = organVolM(age, sex, height, bodyMass, BSA); % volume and mass distribution of organs
q = organFlow(age, sex, height, BSA); % flow of blood into the organs


%% Initialize parameters

%EtOH
lipophilicityEtOH = -0.31; %logP EtOH
fracUnboundEtOH = .99; % fraction unbound EtOH
molarMassEtOH = 46.0684;
pEthanol = 0.789;
volDrink = gkg*bodyMass/(pEthanol*drinkPercent*1000); % volume of drink in litres
volStomLumen = 1.10; % volume of gastric fluid in L
KEtOH = organPartition(lipophilicityEtOH, fracUnboundEtOH); % partitions of the organs

%ACALD
lipophilicityACALD = -0.34; %logP ACALD
fracUnboundACALD = .99; % fraction unbound ACALD
molarMassACALD = 44.05; 
KACALD = organPartition(lipophilicityACALD, fracUnboundACALD); % partitions of the organs

c0 = zeros (34,1);
c0(15) = 0.8*gkg*bodyMass/molarMassEtOH*1000/(volStomLumen); %mMol, initial EtOH into stomach, 80% bioavailability
c0(19) = 0; %no initial blood ACALD concentration

%% Modifying Rxn Bounds

if (ModelVer == 2)
    %Load model and add missing reactions
    load('Harvey_1_03c.mat');
    male = addReaction(male,'Skin_EX_etoh(c)_[bc]', 'reactionFormula', 'Skin_etoh[c]  <=> etoh[bc]', 'printLevel',0);
    male = addReaction(male,'Lung_EX_etoh(br)_[bc]', 'reactionFormula', 'etoh[br]  <=> etoh[bc]', 'printLevel',0);
    male = addReaction(male,'EX_etoh[br]', 'reactionFormula', 'etoh[br] -> ', 'printLevel',0);
   male = addReaction(male,'Pancreas_EX_acald[bc]_[luP]','reactionFormula','acald[bc] <=> Pancreas_acald[luP]','printLevel',0);
   male = addReaction(male,'Liver_EX_acald[bpC]_[bpL]','reactionFormula','Colon_acald[bpC] <=> Liver_acald[bpL]','printLevel',0);
    %adding on slack_ variables
    male = addCOBRAConstraints(male,{'Skin_EX_etoh(c)_[bc]'},0,'c',[1], 'dsense','L','ConstraintID','slack_Skin_EX_etoh(c)_[bc]');
    male = addCOBRAConstraints(male,{'Lung_EX_etoh(br)_[bc]'},0,'c',[1], 'dsense','L','ConstraintID','slack_Lung_EX_etoh(br)_[bc]');
   male = addCOBRAConstraints(male,{'Pancreas_EX_acald[bc]_[luP]'},0,'c',[1], 'dsense','L','ConstraintID','Pancreas_EX_acald[luP]');
   male = addCOBRAConstraints(male,{'Liver_EX_acald[bpC]_[bpL]'},0,'c',[1], 'dsense','L','ConstraintID','slack_Liver_EX_acald[bpL]');


%List of EtOH & ACALD related reactions
EtOH_rxns = [1482;1607;1884;3664;4467;4664;4841;5010;5563;5920;9099;9508;9666;9909;10053;10617;13521;13723;14058;14344;14632;14711;17276;18452;18976;19500;24014;30871;31521;31522;32552;32802;56999;57891;58551;60357;61125;62318;63555;63922;77223;77224;77225;77381;77386];
ACALD_rxns = [1393,1394,1482,1485,1486,1607,1854,3896,4416,4794,5220,5494,5836,5920,5922,6186,6644,7035,9045,9100,10618,10754,10905,11411,11806,14287,18453,24015,27019,32803,34199,39992,43159,48388,50162,52468,53545,54747,55819,56659,56857,58421,60934,69647,72324,76835,80719,80978]';

%fix for transport rate limit
male.C(find(male.C == 20000)) = 70000;
male.C(find(male.C == -20000)) = -70000;

male = changeObjective(male, 'Liver_ALDD2x',1); %maximize ADH elimination
else
    FBA=0;
end

%for tracking how much is excreted in each process
global molliv molc molu molbr molsw
molliv = 0;
molc = 0;
molu = 0;
molbr =0;
molsw =0;


%% Solver


itMax = tMax*10;
t = 0:itMax;
C = zeros(tMax+1,34);
dC = zeros(tMax,34);
step = tMax/itMax;
C(1,:) = c0; %initial conditions

global iteration FBAcount
FBAcount = 0;

for iteration = 1:itMax
    dC(iteration,:)=ODE(C(iteration, :)', age, sex, KEtOH,KACALD, v, q);
    if n~=0 
        for j = 1:n
            t_new = multi(j,1);
            if iteration == t_new*10
                fprintf('dose')%should v this 1000 be here?? check units                                
                C((iteration),15) = multi(j,2)/molarMassEtOH*1000/(volStomLumen) + C((iteration),15);
            end
        end
    end
    C(iteration+1,:) = C(iteration,:) + dC(iteration,:)*step;
end

%% Plotting
figure (1)

x=t/(itMax/tMax);y=C;
subplot (2,1,1)
plot (x, y(:,2), '-');
xlabel('time(min)','FontSize',14);
ylabel('concentration (mmol/L)', 'FontSize',14);
title('EtOH');
hold on

subplot (2,1,2)
plot (x, y(:,19)*1000, '-');
ylabel('concentration (umol/L)','FontSize',14);
xlabel('time(min)','FontSize',14);
title('ACALD');
hold on

%JonesData
%AverageError (jonesx, jonesacald, y(:,19), 10)
Plot_Graphs;
% AUC(x, y(:,2)) %uM*min/L, Blood
% AUC(x, y(:,19)*1000) %uM*min/L, Blood
AUC(x, y(:,8)) %uM*min/L, Liver
AUC(x, y(:,25)*1000) %uM*min/L, Liver


sound(sin(1:3000))

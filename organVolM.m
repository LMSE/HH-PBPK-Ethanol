function v = organVolM(age, sex, height, bodyMass, BSA)
 global fat
 % mass (kg), correlations from Felix Stader 2019, lumens have no mass

 m =[0.01*fat*bodyMass; % 1.adipose
     exp(.067*BSA-.0025*age-.38*sex+1.7);  % 2.blood
     exp(-0.0075*age+0.0078*height-0.97);  % 3.brain
     .45*3E-6*height^2.49;                  % 4.sint, 
     .55*3E-6*height^2.49;                  % 5.lint
     0.34*BSA+0.0018*age-0.36;             % 6.heart
     -0.00038*age-0.056*sex+0.33;          % 7.kidney
     exp(0.87*BSA - 0.0014*age - 1);       % 8.liver
     exp(0.028*height+0.0077*age-5.6);     % 9.lung
     17.9*BSA-0.0667*age-5.68*sex-1.22;    %10.muscle
     0.103;                                %11.pancreas
     exp(-0.0058*age-0.37*sex+1.13);       %12.skin, dermis
     exp(1.13*BSA-3.93);                   %13.spleen
     1.05;                                 %14.stomach, not from Stader
     0;                                    %15.stomach_lumen
     0;                                    %16.sint_lumen
     0];                                   %17.lint_lumen

 %redistributing leftover weight to intestines and muscles
fix = bodyMass - sum(m);
m(4) = m(4)+0.1*fix;
m(5) = m(5)+0.1*fix;
m(12) = m(12)+0.8*fix;
 
 % densities (kg/L)
 
p = [0.916; % 1.adipose
     1.060; % 2.blood
     1.035; % 3.brain
     1.044; % 4.gut
     1.044; % 5.gut
     1.030; % 6.heart     
     1.050; % 7.kidney     
     1.080; % 8.liver, [Heinemann et al]     
     1.050; % 9.lung     
     1.041; %10.muscle     
     1.045; %11.pancreas     
     1.116; %12.skin, dermis (epidermis & hypodermis average out to ~same)     
     1.054; %13.spleen
     1.050; %14.stomach     
     1;     %15.stomach_lumen
     1;     %16.sint_lumen
     1];    %17.lint_lumen     
v = m./p;

end


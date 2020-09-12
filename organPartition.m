function K = organPartition(logP, fracUnbound)
%% Partition Coefficients [from Peters Textbook]
    
 %   V_nlp  V_php  V_water
   
CK =[.79    .002   .180; % 1.Adipose
     .0035  .00225 .945; % 2.blood/plasma
     .051   .0565  .770; % 3.brain
     .0487  .0163  .718; % 4.SInt
     .0487  .0163  .718; % 5.LInt
     .0115  .0166  .758; % 6.heart
     .0207  .0162  .783; % 7.kidney     
     .0348  .0252  .751; % 8.liver     
     .003   .009   .811; % 9.lung     
     .0238  .0072  .760; %10.muscle
     .0723  .0188  .660; %11.pancreas     
     .0284  .0111  .718; %12.skin     
     .0201  .0198  .788; %13.spleen     
     .0338  .0182  .784; %14.stomach
     0 0 0;              %15.stomach_lumen
     0 0 0;              %16.sint_lumen
     0 0 0];             %17.lint_lumen
     
    POW=10^logP;
    fut=1/(1+((1-fracUnbound)*0.5/(fracUnbound)));
    
    % non adipose
    K(:,1)=((POW*(CK(:,1)+.3*CK(:,2))+CK(:,3)+.7*CK(:,2))/(POW*(CK(2,1)+0.3*CK(2,2))+CK(2,3)+.7*CK(2,2))*(fracUnbound/fut));
    % adipose, assume logD ~ logP and fut = 1 for adipose due to minimal binding proteins
    K(1,1)=((POW*(CK(1,1)+.3*CK(1,2))+CK(1,3)+.7*CK(1,2))/(POW*(CK(2,1)+0.3*CK(2,2))+CK(2,3)+.7*CK(2,2))*(fracUnbound/1));

end


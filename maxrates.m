function r = maxrates(age, sex, cIn)

    %Vm(mM/min) Km (mM)      
VmKm = [
% Alcohol Dehydrogenase    
        0       1;      % 1.adipose
        0       1;      % 2.blood
        0       1;      % 3.brain   
        0       1;      % 4.s int
        0       1;      % 5.l int
        0       1;      % 6.heart
        0       1;      % 7.kidney
       2.2   1  ;      % 8.liver, Umulis et al
        0       1;      % 9.lung
        0       1;      %10.muscle
        0       1;      %11.pancreas
        0       1;      %12.skin
        0       1;      %13.spleen
        .68     41;     %14.stomach, Toroghi et al
        0       1;      %15.stomach_lumen
        0       1;      %16.sint_lumen
        0       1;     %17.lint_lumen
% Aldehyde Dehydrogenase
        0       1;      %18.adipose
        0       1;      %19.blood
        0       1;      %20.brain   
        0       1;      %21.s int
        0       1;      %22.l int
        0       1;      %23.heart
        0       1;      %24.kidney
        2.7    1.6;     %25.liver, Umulis et al
        0       1;      %26.lung
        0       1;      %27.muscle
        0       1;      %28.pancreas
        0       1;      %29.skin
        0       1;      %30.spleen
        0       1;      %31.stomach
        0       1;      %32.stomach_lumen
        0       1;      %33.sint_lumen
        0       1];     %34.lint_lumen

    
    
    
    
r(:,1) = (VmKm(:,1).*(cIn))./(VmKm(:,2)+cIn); %[mM/min]
    
if sex == 0
    if age <= 25
        adh = 1;
    elseif age > 55
        adh = 0.5;
    else
        adh = -0.00102*age^2 + 0.067*age - 0.0663;
    end
else
    if age <= 25
        adh = 1;
    elseif age > 55
        adh = .38;
    else
        adh = -0.0003*age^2 + 0.0168*age + .397;
    end
end

r = adh*r; %age and sex effects on ADH enzyme

end


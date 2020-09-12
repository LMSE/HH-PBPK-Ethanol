function [A]=AUC(x,y)
    n = length(x);
    A=0;
    for i=1:(n-2)
       A = A+0.5*(x(i+1)-x(i))*(y(i+1)+y(i));
    end
end

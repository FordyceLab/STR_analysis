function [KDmicros] = analyticalKDs(params)
%Given experiments of a = KD,both,tot , b= KD,core,tot c =KD,flank,tot
%This fhucntion returns the analyical solution for the microscopic KDs of
%testingstate, core and flank state [1 2 3]

a=params(1);
b=params(2);
c=params(3);
KDmicros=zeros(3,1);

KDmicros(1)=a*b*c/(a*b+a*c-b*c);
KDmicros(2)=(-a*b-a*c+b*c)/(b*(a-c));
KDmicros(3)=(-a*b-a*c+b*c)/(c*(a-b));

end


function [c,ceq] = paramsConstr(params,paramIndex,paramIndexEnd,toleratedfacs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
ceq=0;
c=zeros(1,4);

for i=1:4
    paramsPart=params(paramIndex(i):paramIndexEnd(i));
    minVal=min(paramsPart);
    diffFac=paramsPart/minVal;
    c(i)=max(diffFac)-toleratedfacs(i);
end



end


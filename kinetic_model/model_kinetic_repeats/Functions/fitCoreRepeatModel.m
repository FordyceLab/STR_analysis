function [bestParams,bestResidual] = fitCoreRepeatModel(koffData,KDData,weights,seqIdentifier,regionsAdded,paramIndex,numGuesses)
%FITCOREREPEATMODEL(KOFFDATA,KDDATA,SEQIDENTIFIER,PARAMINDEX); Summary of this function goes here
%   Detailed explanation goes here
numParams=1+1+2*regionsAdded(1)+2*regionsAdded(2);
paramIndexEnd=paramIndex(2:end)-1;
paramIndexEnd(4)=numParams;

lb=0*ones(1,numParams);
ub=inf*ones(1,numParams);
guessMin=10^-6*ones(1,numParams);
guessMax=10^6*ones(1,numParams);

lb(1)=0;
lb(2)=10^6;

lb(paramIndex(3):paramIndexEnd(3))=0; %p_core_rel 
lb(paramIndex(4):paramIndexEnd(4))=0; %p_flank_rel

ub(2)=10^6;

%ub(1)=0.1; %konmax
%ub(2)=10^6; %koff,M
%ub(paramIndex(1):paramIndexEnd(1))=10^2; %koffmicro_cores
%ub(paramIndex(2):paramIndexEnd(2))=10^2; %koffmicro_flank 
%ub(paramIndex(3):paramIndexEnd(3))=100; %p_core_rel 
%ub(paramIndex(4):paramIndexEnd(4))=100; %p_flank_rel




guessMin(1)=10^-6; %konmax
guessMin(2)=10^6; %koff,M
guessMin(paramIndex(1):paramIndexEnd(1))=10^-2; %koffmicro_cores
guessMin(paramIndex(2):paramIndexEnd(2))=10^-2; %koffmicro_flank 
guessMin(paramIndex(3):paramIndexEnd(3))=10^-1; %p_core_rel 
guessMin(paramIndex(4):paramIndexEnd(4))=10^-1; %p_flank_rel


guessMax(1)=0.01; %konmax
guessMax(2)=10^6; %koff,M
guessMax(paramIndex(1):paramIndexEnd(1))=10^2; %koffmicro_cores
guessMax(paramIndex(2):paramIndexEnd(2))=10^3; %koffmicro_flank 
guessMax(paramIndex(3):paramIndexEnd(3))=10; %p_core_rel 
guessMax(paramIndex(4):paramIndexEnd(4))=10; %p_flank_rel


funDiff = @(x) diffSquaredCoreRepeatModel(x,koffData,KDData,weights,seqIdentifier,paramIndex);

options=optimoptions('fmincon','Display','final','TolX',1e-7,'TolFun',1e-7,'MaxFunctionEvaluations',10^6,'MaxIter',10^5);
paramsAll=zeros(numGuesses,numParams);
resAll=zeros(numGuesses,1);
parfor i= 1:numGuesses
    
    if i==1
        params0=lb;
    elseif i==2
        params0=guessMax*2;
    else
        %params0=guessMin+(guessMax-guessMin).*rand(1,numParams);
        params0=guessMin.*10.^((log10(guessMax)-log10(guessMin)).*rand(1,numParams));
    end
    
    [paramsCurr,res] = fmincon(funDiff,params0,[],[],[],[],lb,ub,[],options);
   
     resAll(i)=res;
   
        
    paramsAll(i,:)=paramsCurr;
    
    
end

[bestResidual,bestind]=min(resAll);
bestParams=paramsAll(bestind,:);

end


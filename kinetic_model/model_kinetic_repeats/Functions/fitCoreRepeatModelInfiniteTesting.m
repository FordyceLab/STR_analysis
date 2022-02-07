function [bestParams,bestResidual,paramsAll,resAll] = fitCoreRepeatModelInfiniteTesting(koffData,KDData,weights,seqIdentifier,regionsAdded,paramIndex,numGuesses,lb,ub)
%FITCOREREPEATMODEL(KOFFDATA,KDDATA,SEQIDENTIFIER,PARAMINDEX); Summary of this function goes here
%   Detailed explanation goes here
numParams=1+2*regionsAdded(1)+2*regionsAdded(2);
paramIndexEnd=paramIndex(2:end)-1;
paramIndexEnd(4)=numParams;





guessMin(1)=10^-6; %konmax
guessMin(paramIndex(1):paramIndexEnd(1))=10^-3; %koffmicro_core
guessMin(paramIndex(2):paramIndexEnd(2))=10^-3; %koffmico_flank 
guessMin(paramIndex(3):paramIndexEnd(3))=10^-2; %p_core_rel 
guessMin(paramIndex(4):paramIndexEnd(4))=10^-2; %p_flank_rel


guessMax(1)=0.01; %konmax
guessMax(paramIndex(1):paramIndexEnd(1))=10^3; %koffmicro_core
guessMax(paramIndex(2):paramIndexEnd(2))=10^3; %koff_micro_flank
guessMax(paramIndex(3):paramIndexEnd(3))=100; %p_core_rel 
guessMax(paramIndex(4):paramIndexEnd(4))=100; %p_flank_rel


funDiff = @(x) diffSquaredCoreRepeatModelInfiniteTesting(x,koffData,KDData,weights,seqIdentifier,paramIndex);

options=optimoptions('fmincon','Display','off','TolX',1e-7,'TolFun',1e-7,'MaxFunctionEvaluations',10^5,'MaxIter',10^4);
paramsAll=zeros(numGuesses,numParams);
resAll=inf*ones(1,numGuesses);
parfor i= 1:numGuesses
    warning('off');
    if i==1
        params0=guessMin/2;
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


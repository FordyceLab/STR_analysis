function [rates,boundFracOn,boundFracOff] = modelRatesVecTimesODEs(kaMax,kOffMicroDeep,pvals,invTimes,dissocConc,assocConc,dissocMeasTimes,assocMeasTimes,retslope,analytic)
%This function return rates=[kd,ka] for a multistep sequential binding
%process parameterised by (kaMax,kOffMicroDeep,pvals), where pvals is a
%vector with probabilities, pvals(i)= probability to transiation from state
%i+1 to i, where state 1 is the completeley specifically bound state (with microscopic offrate kOffMicroDeep)
%Times:
%%Times Extended model, so that every state has a time spent in this state
%(effectively rates of transitions, not transition probabilities), solved
%by solving the system of ODEs instead of mean first passage time (DNA
%bound is all state 1 to n-1)

nsteps=numel(pvals);

dim=nsteps+2-1;



rateMatrAoff=getRateMatrODEs(pvals,[kOffMicroDeep invTimes]);
rateMatrAoffSteady=rateMatrAoff;
rateMatrAoffSteady(dim,1:dim)=rateMatrAoffSteady(dim,1:dim)-kaMax*dissocConc;%Correct to get the concentration of the dissocation experiment


rateMatrAon=rateMatrAoff;
rateMatrAon(dim,1:dim)=rateMatrAon(dim,1:dim)-kaMax*assocConc;


bon=zeros(dim,1);
bon(dim)=kaMax*assocConc;
%boff=zeros(dim,1);
boffSteady=zeros(dim,1);
boffSteady(dim)=kaMax*dissocConc;
xzeroOn=zeros(dim,1);
xstarOffSteady=-rateMatrAoffSteady\boffSteady;
xzeroOff=xstarOffSteady;
lastwarn('');
if analytic==1
    slopeon=sum(rateMatrAon*xzeroOn+bon);
    slopeoff=sum(rateMatrAoff*xzeroOff);
else

    
    xstarOn=-rateMatrAon\bon;
    
    %xstarOff=-rateMatrAoff\boff;
    xstarOff=zeros(dim,1);
    
    [Von,Don]=eig(rateMatrAon);
    [Voff,Doff]=eig(rateMatrAoff);
    
    numAssoc=numel(assocMeasTimes);
    numDissoc=numel(dissocMeasTimes);
    
    funcValsOn=zeros(dim,numAssoc);
    funcValsOff=zeros(dim,numDissoc);
    for i=1:numAssoc
        expMatr=Von*diag(exp(diag(Don*assocMeasTimes(i))))/Von;%
        %expMatr=expm(rateMatrAon*assocMeasTimes(i));
        funcValsOn(:,i)=xstarOn+expMatr*(xzeroOn-xstarOn);
    end
    for i=1:numDissoc
        expMatr=Voff*diag(exp(diag(Doff*dissocMeasTimes(i))))/Voff;%
        %expMatr=expm(rateMatrAoff*dissocMeasTimes(i));
        funcValsOff(:,i)=xstarOff+expMatr*(xzeroOff-xstarOff);
    end
    
    boundFracOn=sum(funcValsOn(1:dim,:))';
    boundFracOff=sum(funcValsOff(1:dim,:))';
    
    %y intercept is 0
    slopeon=regress(boundFracOn,assocMeasTimes);
    
    polyoff=regress(boundFracOff, [dissocMeasTimes ones(numDissoc,1)]);
    slopeoff=polyoff(1);
end
ka=slopeon/assocConc;
if retslope==1
    rates=[-slopeoff,ka];
    
else
   %kd=-slopeoff/(1+slopeoff/(ka*dissocConc));
   kd=-slopeoff/sum(xzeroOff);
   rates=[kd,ka]; 
    
end



if ~isempty(lastwarn())
    rates=[inf,inf];
else


end        

end


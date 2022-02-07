function [KD,koff,kon] = getRatesCoreRepeatComparing(params,times,solutionConc)
%getRatesCoreRepeat gets equllibrium KD and effeictice rates koff and kon,
%for the  4 state model described by params
%params = [konmax invTimeTestingstate koffmicro_core koffmicro_flank p_core p_flank]
% index 0: Free protein, Index 1: testing state, index 2:core index 3: flank
%The free / not bound protein state is not handled explicitly, but it is
%instead used that x(0)=1-x(1)-x(2)-x(3). There by the non-zero b vector for the assocation experiment,
%Using this implementation and notations give an invertible A matrix which
%is conventient
slopePoints=5;

bon = [params(1)*solutionConc;0;0];
boff = [0;0;0];

%bOff = [-params(1)*solutionConc;0;0];

%A=[-params(1)*solutionConc params(2)*(1-params(5)-params(6)) 0 0;params(1)*solutionConc -params(2) params(3) params(4);0 params(2)*params(5) -params(3) 0;0 params(2)*params(6) 0 -params(4)];
A=[-params(2) params(3) params(4);params(2)*params(5) -params(3) 0;params(2)*params(6) 0 -params(4)];
Aoff=A;
%dx(1)/dt should also have a term params(1)*solutionConc*x(0), that not handled in the experssiomn above. Correct by adding and using x(0)=1-x(1)-(x2)-x(3) on the first row
A(1,:)= A(1,:)-params(1)*solutionConc;
Aon=A;

x0on=[0;0;0];
%Add an extra line saying all concentrations should sum to 1

x0off=Aon\(-bon);

funcValsOn=solveODEs(Aon,bon,x0on,times);
funcValsOff=solveODEs(Aoff,boff,x0off,times);

onCurve=sum(funcValsOn(2:3,:));

theSum=sum(x0off);

KD = solutionConc*(1-theSum)/theSum;


%survivalFuncOff=sum(funcValsOff)/theSum;
survivalFuncOff=sum(funcValsOff(2:3,:));
survivalFuncOff=survivalFuncOff/survivalFuncOff(1);
%Formula for expecation value using the cdf (1-survival function)
offTime=trapz(times,survivalFuncOff);
koffIntegral=1/offTime;

p=polyfit(times(1:slopePoints),survivalFuncOff(1:slopePoints),1);
koffSlope=-p(1);
p=polyfit(times(1:slopePoints),onCurve(1:slopePoints),1);
konSlope=p(1)/solutionConc;

%Try with mean first passagetimes

%tiem from dissocatied to core OR flank

onMeanFirstPassage=((1/(params(1)*solutionConc))+1/params(2))/(params(5)+params(6));

konMeanFirstPassage=1/(onMeanFirstPassage*solutionConc);

offMeanFirstPassageCore=((1-params(6))*(1/params(3))+params(6)*(1/params(4))+(1/params(2)))/(1-params(5)-params(6));

offMeanFirstPassageFlank=((1-params(5))*(1/params(4))+params(5)*(1/params(3))+(1/params(2)))/(1-params(5)-params(6));

offMeanFirstPassage=x0off(2)/(x0off(2)+x0off(3))*offMeanFirstPassageCore+x0off(3)/(x0off(2)+x0off(3))*offMeanFirstPassageFlank;
offMeanFirstPassageInitial=params(5)/(params(5)+params(6))*offMeanFirstPassageCore+params(6)/(params(5)+params(6))*offMeanFirstPassageFlank;

koffMeanFirstPassage=1/offMeanFirstPassage;

koffMeanFirstPassageInitial=1/offMeanFirstPassageInitial;


end


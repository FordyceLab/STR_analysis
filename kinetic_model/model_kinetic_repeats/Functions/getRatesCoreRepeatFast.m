function [KD,koff] = getRatesCoreRepeatFast(params,solutionConc)
%getRatesCoreRepeat gets equllibrium KD and effeictice rates koff and kon,
%for the  4 state model described by params
%params = [konmax invTimeTestingstate koffmicro_core koffmicro_flank p_core p_flank]
% index 0: Free protein, Index 1: testing state, index 2:core index 3: flank
%The free / not bound protein state is not handled explicitly, but it is
%instead used that x(0)=1-x(1)-x(2)-x(3). There by the non-zero b vector for the assocation experiment,
%Using this implementation and notations give an invertible A matrix which
%is conventient. Note that koff and kon are inverses of mean first passage
%times at equllibrium, so that koff != koffIntitial, and KD!=koff/kon
slopePoints=5;

bon = [params(1)*solutionConc;0;0];



%Add an extra line saying all concentrations should sum to 1

x0off=Aon\(-bon);

theSum=sum(x0off);

KD = solutionConc*(1-theSum)/theSum;





%Try with mean first passagetimes

%time from dissocatied to core OR flank

%onMeanFirstPassage=((1/(params(1)*solutionConc))+1/params(2))/(params(5)+params(6));

%konMeanFirstPassage=1/(onMeanFirstPassage*solutionConc);

offMeanFirstPassageCore=((1-params(6))*(1/params(3))+params(6)*(1/params(4))+(1/params(2)))/(1-params(5)-params(6));

offMeanFirstPassageFlank=((1-params(5))*(1/params(4))+params(5)*(1/params(3))+(1/params(2)))/(1-params(5)-params(6));

offMeanFirstPassage=x0off(2)/(x0off(2)+x0off(3))*offMeanFirstPassageCore+x0off(3)/(x0off(2)+x0off(3))*offMeanFirstPassageFlank;
%offMeanFirstPassageInitial=params(5)/(params(5)+params(6))*offMeanFirstPassageCore+params(6)/(params(5)+params(6))*offMeanFirstPassageFlank;

koffMeanFirstPassage=1/offMeanFirstPassage;

%koffMeanFirstPassageInitial=1/offMeanFirstPassageInitial;

%kon=konMeanFirstPassage;
koff=koffMeanFirstPassage;

end


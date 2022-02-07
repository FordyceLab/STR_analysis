function [KD,koff] = getRatesCoreRepeatFastConvertParamsKDmicro(params,solutionConc)
%getRatesCoreRepeat gets equllibrium KD and effeictice rates koff and kon,
%for the  4 state model described by params
%input params=[%konmax %KDmicro,M,testing KDmicro_core KDmicro_flank p_core_rel p_flank_rel]
%But is converted in beginning of function to
%params = [konmax invTimeTestingstate koffmicro_core koffmicro_flank p_core p_flank]
% index 0: Free protein, Index 1: testing state, index 2:core index 3: flank
%The free / not bound protein state is not handled explicitly, but it is
%instead used that x(0)=1-x(1)-x(2)-x(3). There by the non-zero b vector for the assocation experiment,
%Using this implementation and notations give an invertible A matrix which
%is conventient. Note that koff and kon are inverses of mean first passage
%times at equllibrium, so that koff != koffIntitial, and KD!=koff/kon
paramsOrig=params;
koffM=params(2)*params(1);
invTimeTestingstate=koffM+koffM*params(5)+koffM*params(6);
p_core=koffM*params(5)/invTimeTestingstate;
p_flank=koffM*params(6)/invTimeTestingstate;
koffmicro_core=params(3)*koffM*params(5);
koffmicro_flank=params(4)*koffM*params(6);
if params(5)==0
    koffmicro_core=1;%doesn't matter just not 0
end
if params(6)==0
    koffmicro_flank=1;%doesn't matter just not 0
end


params(2)=invTimeTestingstate;
params(5)=p_core;
params(6)=p_flank;
params(3)=koffmicro_core;
params(4)=koffmicro_flank;

% % % % A=[-params(2) params(3) params(4);params(2)*params(5) -params(3) 0;params(2)*params(6) 0 -params(4)];
% % % % %dx(1)/dt should also have a term params(1)*solutionConc*x(0), that not handled in the experssiomn above. Correct by adding and using x(0)=1-x(1)-(x2)-x(3) on the first row
% % % % A(1,:)= A(1,:)-params(1)*solutionConc;
% % % % Aon=A;
% % % % 
% % % % bon = [params(1)*solutionConc;0;0];
% % % % 
% % % % 
% % % % 
% % % % %Add an extra line saying all concentrations should sum to 1
% % % % lastwarn('', '');
% % % % x0off=Aon\(-bon);
% % % %  [warnMsg, warnId] = lastwarn();
% % % %  
% % % %  if ~isempty(warnId)
% % % %      aa=1 +1  ;
% % % %      KD=NaN;
% % % %      koff=NaN;
% % % %  else
% % % %      theSum=sum(x0off);
% % % %      
% % % %      KD = solutionConc*(1-theSum)/theSum;
% % % %      

if p_core==0 %only
   KD =paramsOrig(2)*paramsOrig(4)/(paramsOrig(4)+1);
elseif p_flank==0
    
    KD =paramsOrig(2)*paramsOrig(3)/(paramsOrig(3)+1);
else
    KD = paramsOrig(2)*paramsOrig(3)*paramsOrig(4)/(paramsOrig(3)*paramsOrig(4)+paramsOrig(4)+paramsOrig(3));
    
end
     
     
     
     
     %Try with mean first passagetimes
     
     %time from dissocatied to core OR flank
     
     %onMeanFirstPassage=((1/(params(1)*solutionConc))+1/params(2))/(params(5)+params(6));
     
     %konMeanFirstPassage=1/(onMeanFirstPassage*solutionConc);
     
     offMeanFirstPassageCore=((1-params(6))*(1/params(3))+params(6)*(1/params(4))+(1/params(2)))/(1-params(5)-params(6));
     
     offMeanFirstPassageFlank=((1-params(5))*(1/params(4))+params(5)*(1/params(3))+(1/params(2)))/(1-params(5)-params(6));
     
     %offMeanFirstPassage=x0off(2)/(x0off(2)+x0off(3))*offMeanFirstPassageCore+x0off(3)/(x0off(2)+x0off(3))*offMeanFirstPassageFlank;
     offMeanFirstPassage=paramsOrig(4)/(paramsOrig(3)+paramsOrig(4))*offMeanFirstPassageCore+paramsOrig(3)/(paramsOrig(3)+paramsOrig(4))*offMeanFirstPassageFlank;
     %offMeanFirstPassageInitial=params(5)/(params(5)+params(6))*offMeanFirstPassageCore+params(6)/(params(5)+params(6))*offMeanFirstPassageFlank;
     
     koffMeanFirstPassage=1/offMeanFirstPassage;
     
     %koffMeanFirstPassageInitial=1/offMeanFirstPassageInitial;
     
     %kon=konMeanFirstPassage;
     koff=koffMeanFirstPassage;
 end



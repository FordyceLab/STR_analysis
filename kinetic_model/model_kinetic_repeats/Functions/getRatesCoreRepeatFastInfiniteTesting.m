function [KD,koff] = getRatesCoreRepeatFastInfiniteTesting(params)
%getRatesCoreRepeat gets equllibrium KD and effeictice rates koff and kon,
%for the  4 state model described by params
%params = [konmax koffmicro_core koffmicro_flank p_core_rel p_flank_rel]
%But is converted in beginning of function to
%params = [konmax invTimeTestingstate koffmicro_core koffmicro_flank p_core p_flank]
% index 0: Free protein, Index 1: testing state, index 2:core index 3: flank
%The free / not bound protein state is not handled explicitly, but it is
%instead used that x(0)=1-x(1)-x(2)-x(3). There by the non-zero b vector for the assocation experiment,
%Using this implementation and notations give an invertible A matrix which
%is conventient. Note that koff and kon are inverses of mean first passage
%times at equllibrium, so that koff != koffIntitial, and KD!=koff/kon
%paramsOrig=params;

p_core = params(4)/(1+params(4)+params(5));
p_flank= params(5)/(1+params(4)+params(5));

%params(4)=p_core;
%params(5)=p_flank;

%if params(4)==0
%    koffmicro_core=1;%doesn't matter just not 0
%end
%if params(5)==0
%    koffmicro_flank=1;%doesn't matter just not 0
%end


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
% 
% if p_core==0 %only
%    KD =paramsOrig(2)*paramsOrig(4)/(paramsOrig(4)+1);
% elseif p_flank==0
%     
%     KD =paramsOrig(2)*paramsOrig(3)/(paramsOrig(3)+1);
% else
%     KD = paramsOrig(2)*paramsOrig(3)*paramsOrig(4)/(paramsOrig(3)*paramsOrig(4)+paramsOrig(4)+paramsOrig(3));
%     
% end
     
     
     
     
     %Try with mean first passagetimes
     
     %time from dissocatied to core OR flank
     
     onMeanFirstPassage=((1/(params(1))))/(p_core+p_flank);
     
     konMeanFirstPassage=1/(onMeanFirstPassage);
     
     offMeanFirstPassageCore=((1-p_flank)*(1/params(2))+p_flank*(1/params(3)))/(1-p_core-p_flank);
     
     offMeanFirstPassageFlank=((1-p_core)*(1/params(3))+p_core*(1/params(2)))/(1-p_core-p_flank);
     
     %offMeanFirstPassage=x0off(2)/(x0off(2)+x0off(3))*offMeanFirstPassageCore+x0off(3)/(x0off(2)+x0off(3))*offMeanFirstPassageFlank;
     
     KDmicroCoreProp=params(2)/p_core; %Something that is proportional to KDMicroCore
     KDmicroFlankProp=params(3)/p_flank;%Also proportional (same prop constant as above)
     
     fracCoreEq=KDmicroFlankProp/(KDmicroCoreProp+KDmicroFlankProp);
     fracFlankEq=KDmicroCoreProp/(KDmicroCoreProp+KDmicroFlankProp);
     
     if p_core==0
         fracCoreEq=0;
         fracFlankEq=1;
         
     elseif p_flank==0
         fracCoreEq=1;
         fracFlankEq=0;
     end
     
     offMeanFirstPassage=fracCoreEq*offMeanFirstPassageCore+fracFlankEq*offMeanFirstPassageFlank;
     
     offMeanFirstPassageInitial=p_core/(p_core+p_flank)*offMeanFirstPassageCore+p_flank/(p_core+p_flank)*offMeanFirstPassageFlank;
     
     koffMeanFirstPassageMeasured=1/offMeanFirstPassage;
     
     koffMeanFirstPassageInitial=1/offMeanFirstPassageInitial;
     
     %kon=konMeanFirstPassage;
     koff=koffMeanFirstPassageMeasured;
     
     KD=koffMeanFirstPassageInitial/konMeanFirstPassage;
 end



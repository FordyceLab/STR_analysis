function [onCoreFirstPassage] = getPassagetimeCoreInfiniteTesting(params,solutionConc)
%returns first passage time from dissocatiated to the core for the 4 state
%model defined by
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

t1E = 1/(params(1)*solutionConc);

onCoreFirstPassage = ((1-p_flank)*(t1E)+p_flank*(1/params(3)))/p_core;


 end



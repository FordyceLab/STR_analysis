function [funcVals] = solveODEs(A,b,x0,times)
%Solves a system of ordinary differnatial eqations on the form dx/dt =
%Ax(t)+b , where A is a matrix rates (generator matrix), see https://en.wikipedia.org/wiki/Matrix_differential_equation
    xstar=-A\b;
    dim=numel(b);
   
    
    %[Von,Don]=eig(A);
    
    numTimes=numel(times);
    
    funcVals=zeros(dim,numTimes);
    for i=1:numTimes
        %expMatr=Von*diag(exp(diag(Don*assocMeasTimes(i))))/Von;%
        %expMatr=expm2(rateMatrAon*assocMeasTimes(i));
        expMatr=expm(A*times(i));
        funcVals(:,i)=xstar+expMatr*(x0-xstar);
    end

end


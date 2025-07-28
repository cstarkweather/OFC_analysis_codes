% AUROC: [auc, pval, shuff_auc] = auROC(x, y, N, NUM_ITER)
%  auROC computes the area under the ROC given two distributions x and y, 
%   and optionally does significance tests.  Significance tests (by 
%   permutation tests) are only performed if one of the optional 2nd or 3rd
%   outputs are saved in the call.
%
%  Inputs (required):
%   - x and y: two distributions. Need not be the same size or dimensions.
%  Inputs (optional):
%   - N: number of threshold values used to construct ROC. [def=100]
%   - NUM_ITER: # of shufflings to perform for any signif tests [def=5000]
%  Outputs:
%   - auc: area under ROC curve, such that 1 indicates x > y, reliably.
%   - pval: fraction of shufflings which are further from 0.5 than true auc
%   - shuff_auc: the vector of shuffled auc's.  This is useful for
%                constructing confidence intervals.


function [auc, pval, shuff_auc] = auROC(x, y, N, NUM_ITER)
%% Check variables
if nargin<3 || isempty(N)
    N = 100; %default number of threshold values used to construct ROC.
end
if nargin<4
    NUM_ITER = 5000;
end
if nargout>=2
    testflag = 1; %flag to do the simulations
else
    testflag = 0; 
end

%% Compute ROC and calculate the area under the curve
x = x(:);  %make sure that that 
y = y(:);
lenx = length(x);
leny = length(y);
threshlo = min([x; y]);
threshhi = max([x; y]);
thresh = linspace(threshlo,threshhi,N);
fa = zeros(1,N);	% allocate the false alarm vector
hit = zeros(1,N);
for i = 1:N
  fa(N-i+1) = sum(y > thresh(i));
  hit(N-i+1) = sum(x > thresh(i));
end
fa = fa/leny;
hit = hit/lenx;
fa(1) = 0;
hit(1) = 0;
fa(N) = 1;
hit(N) = 1;
auc = trapz(fa,hit); 

%% if necessary, bootstrap the ROC for significance testing
% The idea is to randomly assign the same values into two groups and
% recompute ROCs.
if testflag
    z = [x; y];
    shuff_auc = zeros(NUM_ITER,1);  %preallocate the shuffled auc vector
    for i = 1:NUM_ITER
        inds = randperm(lenx+leny);
        shuff_auc(i) = auROC( z(inds(1:lenx)), z(inds(1+lenx:end)), N ); %recurse!
    end
    shuff_auc = sort(shuff_auc);
    pval = sum( abs(shuff_auc-0.5) > abs(auc-0.5) ) / NUM_ITER; 
     %i.e., fraction of shuffled AUC's more distant from 0.5 than the true auc is distant from 0.5.
end
    
% if nargout < 2
%     varargout = [];
% elseif nargout == 2
%     varargout{1} = pval;
% elseif nargout == 3
%     varargout{1} = pval;
%     varargout{2} = shuff_auc;
% end

end
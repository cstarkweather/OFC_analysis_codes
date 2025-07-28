% KMEANSINIT: seeds = kmeansinit(data,k)
% This function takes the data to be clustered and the number of clusters
% and computes initial seeds for the k clusters, according to the k-means++
% algorithm: http://en.wikipedia.org/wiki/K-means%2B%2B 
% The idea is to sequentially choose clusters form the data points to be
% the initial seeds, with the probability of a datum being chosen
% increasing in proportion to the squared-distance to the nearest seed
% already chosen.
% 
% Inputs (required)
%  - data: m x n matrix of m data points and n parameters.
%  - k: number of seeds to be produced.
% Output
%  - seeds: k x n matrix of k seeds and n parameters.

% Vinod Rao

function seeds = kmeansinit(data,k)

[m,n] = size(data);

chosen = zeros(k,1);
seeds = zeros(k,n);

% 1) Choose one center uniformly at random from among the data points.
chosen(1) = ceil(rand*m);

fulldist = pdist(data,'euclidean');
fulldist = squareform(fulldist);
fulldist(fulldist==0)=inf;
for i = 2:k
    % 2) For each data point x, compute D(x), the distance between x and
    % the nearest center that has already been chosen.
    D = fulldist(chosen(chosen>0),:);
    if length(D)~=numel(D)
        D = min(D);
    end
    % 3) Choose one new data point at random as a new center, using a weighted
    % probability distribution where a point x is chosen with probability proportional to D(x)2.
    notchosen = setdiff(1:m,chosen);  %ensure new
    D = D(notchosen).^2;              %and square
    prob = D./sum(D);                 %prob
    newseedind = randsample(length(prob),1,true,prob); %weight prob
    chosen(i) = notchosen(newseedind);
end

seeds = data(chosen,:);
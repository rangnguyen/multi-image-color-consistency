function C = extract_brochure_palette(I, k, sigma)
% increase one more to remove dark palette color
k = k + 1;
[m, n, b] = size(I);
I = reshape(I, [], 3);
bin = 15;
colorTransform = makecform('srgb2lab');
I = applycform(I, colorTransform);
[ws, X] = im3Dhist(I,bin);


% Initialize the cluster centers
cinits = zeros(k, size(X,2));
cw = ws;
N = size(X, 1);
sigma2 = sigma^2;
for i = 1:k
    [~,id] = max(cw);
    cinits(i,:) = X(id,:);
    d2 = repmat(cinits(i,:), N, 1) - X; 
    d2 = sum(d2 .* d2,2);
    cw = cw .* (1 - exp(-d2/sigma2));
end


opt.weight = ws;
[~, C, ~] = fkmeans(X, cinits, opt);


% sort by lightness
[~,id] = sort(C(:,1), 'descend');

C = C(id,:);
C = C(1:k-1,:);

colorTransform = makecform('lab2srgb');
C = applycform(C, colorTransform);

function [W, C, labels, ids] = im3Dhist(I, bin)

[W, C, ~, labels] = histcn_modified(I, bin, bin, bin, 'AccumData', I, 'Fun', @mean);
ids = find(W ~= 0);
W = W(ids);
C = C(ids,:);

function [count, meanval, edges, id] = histcn_modified(X, varargin)
% function [count edges mid loc] = histcn(X, edge1, edge2, ..., edgeN)
%
% Purpose: compute n-dimensional histogram
%
% INPUT
%   - X: is (M x N) array, represents M data points in R^N
%   - edgek: are the bin vectors on dimension k, k=1...N.
%     If it is a scalar (Nk), the bins will be the linear subdivision of
%     the data on the range [min(X(:,k)), max(X(:,k))] into Nk
%     sub-intervals
%     If it's empty, a default of 32 subdivions will be used
%
% OUTPUT
%   - count: n-dimensional array count of X on the bins, i.e.,
%         count(i1,i2,...,iN) = cardinal of X such that
%                  edge1(i1) <= X(:,i1) < edge1(i1)+1 and
%                       ...
%                  edgeN(iN) <= X(:,iN) < edgeN(iN)+1
%   - edges: (1 x N) cell, each provides the effective edges used in the
%     respective dimension
%   - mid: (1 x N) cell, provides the mid points of the cellpatch used in
%     the respective dimension
%   - loc: (M x N) array, index location of X in the bins. Points have out
%     of range coordinates will have zero at the corresponding dimension.
%
% DATA ACCUMULATE SYNTAX:
%   [ ... ] = histcn(..., 'AccumData', VAL);
%   where VAL is M x 1 array. Each VAL(k) corresponds to position X(k,:)
%   will be accumulated in the cell containing X. The accumulate result
%   is returned in COUNT.
%   NOTE: Calling without 'AccumData' is similar to having VAL = ones(M,1)
%
%   [ ... ] = histcn(..., 'AccumData', VAL, 'FUN', FUN);
%     applies the function FUN to each subset of elements of VAL.  FUN is
%     a function that accepts a column vector and returns
%     a numeric, logical, or char scalar, or a scalar cell.  A has the same class
%     as the values returned by FUN.  FUN is @SUM by default.  Specify FUN as []
%     for the default behavior.
%
% Usage examples:
%   M = 1e5;
%   N = 3;
%   X = randn(M,N);
%   [N edges mid loc] = histcn(X);
%   imagesc(mid{1:2},N(:,:,ceil(end/2)))
%
% % Compute the mean on rectangular patch from scattered data
%   DataSize = 1e5;
%   Lat = rand(1,DataSize)*180;
%   Lon = rand(1,DataSize)*360;
%   Data = randn(1,DataSize);
%   lat_edge = 0:1:180;
%   lon_edge = 0:1:360;
%   meanData = histcn([Lat(:) Lon(:)], lat_edge, lon_edge, 'AccumData', Data, 'Fun', @mean);
%
% See also: HIST, ACCUMARRAY
% 
% Bruno Luong: <brunoluong@yahoo.com>
% Last update: 25/August/2011

if ndims(X)>2
    error('histcn: X requires to be an (M x N) array of M points in R^N');
end
DEFAULT_NBINS = 32;

AccumData = [];
Fun = {};

% Looks for where optional parameters start
% For now only 'AccumData' is valid
split = find(cellfun('isclass', varargin, 'char'), 1, 'first');
if ~isempty(split)
    for k = split:2:length(varargin)
        if strcmpi(varargin{k},'AccumData')
            AccumData = varargin{k+1};
        elseif strcmpi(varargin{k},'Fun')
            Fun = varargin(k+1); % 1x1 cell
        end
    end
    varargin = varargin(1:split-1);
end

% Get the dimension
nd = size(X,2);
edges = varargin;
if nd<length(edges)
    nd = length(edges); % wasting CPU time warranty
else
    edges(end+1:nd) = {DEFAULT_NBINS};
end

% Allocation of array loc: index location of X in the bins
loc = zeros(size(X));
% Loop in the dimension
for d=1:nd
    ed = edges{d};
    Xd = X(:,d);
    if isempty(ed)
        ed = DEFAULT_NBINS;
    end
    if isscalar(ed) % automatic linear subdivision
        ed = linspace(min(Xd),max(Xd),ed+1);
    end
    edges{d} = ed;
    % Call histc on this dimension
    [~, loc(:,d)] = histc(Xd, ed, 1);
end % for-loop

bin = length(edges{1});
if nd == 2
    id = (loc(:,1)-1)*bin + loc(:,2);
else
    id = (loc(:,1)-1)*bin*bin + (loc(:,2)-1)*bin + loc(:,3);
end

count = accumarray(id, 1);
meanval = zeros(size(count,1), size(AccumData,2));
for i = 1:size(AccumData,2)
    meanval(:,i) = accumarray(id, AccumData(:, i), [], Fun{:});
end





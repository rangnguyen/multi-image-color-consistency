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

[W, C, ~, labels] = histcn(I, bin, bin, bin, 'AccumData', I, 'Fun', @mean);
ids = find(W ~= 0);
W = W(ids);
C = C(ids,:);

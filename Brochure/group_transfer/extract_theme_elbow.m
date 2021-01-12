function [C, H, labels, Wo, Co] = extract_theme_elbow(I, mode, LUT, bin)
% increase one more to remove dark palette color
% k = k + 1;
[m, n, b] = size(I);
I = reshape(I, [], 3);
I = RGB2Lab(I, LUT);

[Wo, Co] = im3Dhist(I,bin, mode);

ids = find(Wo ~= 0);
ws = Wo(ids);
ws = ws/sum(ws);
X = Co(ids,:);
% init distortion 
distortion=zeros(7,1);

% Only 1 cluster
Cold = mean(X, 1);
distortion(1) = sum(ws.*sum((X - repmat(Cold, [size(X,1), 1])).^2,2));
for k = 2:7
    % Initialize the cluster centers
    cinits = zeros(k, size(X,2));
    cw = ws;
    N = size(X, 1);
    for i = 1:k
        [~,id] = max(cw);
        cinits(i,:) = X(id,:);
        d2 = repmat(cinits(i,:), N, 1) - X; 
        d2 = sum(d2 .* d2,2);
%         cw = cw .* (1 - exp(-d2/(80^2)));
        d2 = d2/max(d2);
        cw = cw .* (d2.^2);    
    end
    opt.weight = ws;
    [~, ~, temp] = fkmeans(X, cinits, opt);
    distortion(k) = sum(temp);
end

variance=distortion(1:end-1)-distortion(2:end);
distortion_percent=cumsum(variance)/(distortion(1)-distortion(end));

[r,~]=find(distortion_percent>0.93);
K=r(1,1)+1;
% Find the right one with K cluster
 % Initialize the cluster centers
cinits = zeros(K, size(X,2));
cw = ws;
N = size(X, 1);
for i = 1:K
    [~,id] = max(cw);
    cinits(i,:) = X(id,:);
    d2 = repmat(cinits(i,:), N, 1) - X; 
    d2 = sum(d2 .* d2,2);
%     cw = cw .* (1 - exp(-d2/(80^2)));
    d2 = d2/max(d2);
    cw = cw .* (d2.^2);    
end
opt.weight = ws;
[L, C, ~] = fkmeans(X, cinits, opt);
% C = cinits;
H = zeros(K,1);
for i = 1:K
    H(i) = sum(ws(L == i));
end

H = H / sum(H(:));

% compute the labels
D = zeros(size(I,1), K);
if mode == 2
    for i = 1:K
        D(:,i) = sum((I(:,2:3) - repmat(C(i,:), [size(I,1), 1])).^2, 2);
    end
else
    for i = 1:K
        D(:,i) = sum((I - repmat(C(i,:), [size(I,1), 1])).^2, 2);
    end
end
[~, labels] = min(D, [], 2);
labels = reshape(labels, [m, n]);

% covert to rgb
if mode == 2
    C = [70*ones(K, 1) C];
end
colorTransform = makecform('lab2srgb');
C = applycform(C, colorTransform);

function [W, C] = im3Dhist(I, bin, mode)
if mode == 2
    [W, C, ~, labels] = histcn(I(:,2:3), bin, bin, 'AccumData', I(:,2:3), 'Fun', @mean);
else
    [W, C, ~, labels] = histcn(I, bin, bin, bin, 'AccumData', I, 'Fun', @mean);
end



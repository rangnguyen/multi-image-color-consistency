function [C, H, X, ws] = extract_common_theme(X, ws, k, sigma, mode, LUT, flag)
% increase one more to remove dark palette color
% k = k + 1;
if flag == 0  
    I = []; 
    for i = 1:length(X)    
        I = [I; reshape(X{i}, [], 3)];
    end
    bin = 15;
    I = RGB2Lab(I, LUT);
    clear X;
    [ws, X] = im3Dhist(I,bin, mode);
end


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
[L, C, ~] = fkmeans(X, cinits, opt);
% C = cinits;
H = zeros(k,1);
for i = 1:k
    H(i) = sum(ws(L == i));
end

% sort by dominant
[~,id] = sort(H, 'descend');

C = C(id,:);
H = H(id);
H = H / sum(H(:));

% covert to rgb
if mode == 2
    C = [70*ones(k, 1) C];
end
colorTransform = makecform('lab2srgb');
C = applycform(C, colorTransform);

function [W, C, labels, ids] = im3Dhist(I, bin, mode)
if mode == 2
    [W, C, ~, labels] = histcn(I(:,2:3), bin, bin, 'AccumData', I(:,2:3), 'Fun', @mean);
else
    [W, C, ~, labels] = histcn(I, bin, bin, bin, 'AccumData', I, 'Fun', @mean);
end
ids = find(W ~= 0);
W = W(ids);
C = C(ids,:);
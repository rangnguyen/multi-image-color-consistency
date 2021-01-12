function [Pout, Lout, M] = solve_optimal_all_palette(Pin, reference, W, Hall, Wall, m, gamma, eta)
tic;
n = length(Pin);
Pout = Pin;
P = cell2mat(Pin(:));
colorTransform = makecform('srgb2lab');
for i = 1:n
    Pin{i} = applycform(Pin{i}, colorTransform);
end
if n > 1    
    gamma = gamma /((n-1)/2);
end

if sum(reference) == 0
    reference = reference + 1;
end

no_of_color = 0;
for i = 1:n
    if reference(i)
        no_of_color = no_of_color + size(Pin{i},1);
    end
end
m = min(m, no_of_color);
[Lout, M, ~] = solve_optimal_all_palette_rec(Pin, reference, W, Hall, Wall, m, n, gamma, eta);

colorTransform = makecform('lab2srgb');
M = applycform(M, colorTransform);
% median
% D = pdist2(M,P);
% [~,idx] = min(D, [], 2);
% M = P(idx,:);

% vivid 
% M = choose_vivid(P, Lout, m);

for i = 1:n   
    idx = Lout{i} > 0;
    Pout{i}(idx,:) = M(Lout{i}(idx),:);  
end
toc;


function M = choose_vivid(P, L, m)
M = zeros(m,3);
L = cell2mat(L(:));
[~,S,~] = rgb2hsv(P);
for i = 1:m
    id = L==i;
    si = S(id);
    Pi = P(id,:);
    [~,idx] = max(si);
    M(i,:) = Pi(idx,:);
end

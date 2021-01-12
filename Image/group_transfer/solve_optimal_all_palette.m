function [Pout, Lout, M] = solve_optimal_all_palette(Pin, reference, W, Hall, Wall, m, gamma, eta)
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

for i = 1:n   
    idx = Lout{i} > 0;
    Pout{i}(idx,:) = M(Lout{i}(idx),:);  
end


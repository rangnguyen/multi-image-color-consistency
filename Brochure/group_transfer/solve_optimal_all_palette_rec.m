function [Lout, M, sum_val] = solve_optimal_all_palette_rec(Pin, reference, W, Hall, Wall, m, n , gamma, eta)
% To fast cut down unimpossible path
init_min = Inf();
% load(['penalty_' num2str(m) '.mat']);
% Initialize the cluster centers
% M = init_mean(cell2mat(Pin(:)), cell2mat(W(:)), m);
M = init_mean(Hall, Wall, m);
% k-median
M = choose_median(Pin,M);

% loop 10 times
old_val = 0;
lambda = gamma/50;
for t = 1:10
    sum_val = 0;
    % solve for the assignment
    for i = 1:n
        [Lout{i}, val] = solve_optimal_individual_palette(Pin{i}, W{i}, M, lambda, gamma, eta, init_min);
        sum_val = sum_val + val;
    end
    % re-compute the group color theme (mean colors)
    M = mean_solve(Pin, reference, W, Lout, m, n, lambda);
    disp(['Iteration ' num2str(t) ', val: ' num2str(sum_val)]);
    if abs(old_val - sum_val) < 10
        break;
    else
        old_val = sum_val;
    end
end



function M = init_mean(P, W, m)
if m == 1
    M = mean(P, 1);
else
    cinits = zeros(m, size(P,2));
    cw = W;
    N = size(P, 1);
    sigma = 80;
    sigma2 = sigma^2;
    for i = 1:m
        [~,id] = max(cw);
        cinits(i,:) = P(id,:);
        d2 = repmat(cinits(i,:), N, 1) - P; 
        d2 = sum(d2 .* d2,2);
%         cw = cw .* (1 - exp(-d2/sigma2));
        d2 = d2/max(d2);
        cw = cw .* (d2.^2); 
    end

    % initialize for the group theme
    opt.weight = W;
    [~,M,~] = fkmeans(P, cinits, opt);
end
if size(P,2) == 2
    M = [70*ones(m, 1) M];
end

function M = mean_solve2(P, reference, W, L, m, n , lambda)
A = cell2mat(P(:));
B = cell2mat(W(:));
C = cell2mat(L(:));
M = zeros(m,3);
for i = 1:m
    a = A(C==i,:);
    b = B(C==i,:);
    M(i,:) = sum(a.*repmat(b,1, 3), 1)./repmat(sum(b),1,3);
end
M = choose_median(P,M);

function M = mean_solve(P, reference, W, L, m, n, lambda)
A = zeros(m,m);
B = zeros(m,3);
M = zeros(m,3);

for i = 1:n
    if reference(i) == 0
        continue;
    end
    Pi = P{i}; Wi = W{i}; Li = L{i};
    ni = size(Pi,1);
    % first term
    for j = 1:ni
        if(Li(j) > 0)
            A(Li(j),Li(j)) = A(Li(j),Li(j)) + Wi(j);
            B(Li(j),:) = B(Li(j),:) + Wi(j)*Pi(j,:);
        end
    end
    
    % second term
    for j1 = 1:ni-1
        for j2 = j1+1:ni
            if(Li(j1) > 0 && Li(j2) > 0)
                A(Li(j1),Li(j1)) = A(Li(j1),Li(j1)) + lambda;
                A(Li(j1),Li(j2)) = A(Li(j1),Li(j2)) - lambda;
                A(Li(j2),Li(j2)) = A(Li(j2),Li(j2)) + lambda;
                A(Li(j2),Li(j1)) = A(Li(j2),Li(j1)) - lambda;
                B(Li(j1),:) = B(Li(j1),:) + lambda*(Pi(j1,:)-Pi(j2,:));
                B(Li(j2),:) = B(Li(j2),:) + lambda*(Pi(j2,:)-Pi(j1,:));
            end
        end
    end
end


% solve least squares
for i = 1:3
    M(:,i) = A\B(:,i);
end

M(isnan(M)) = 0;
% k-median
M = choose_median(P,M);

function M = choose_median(Pin, M)
P = cell2mat(Pin(:));
D = pdist2(M,P);
[~,idx] = min(D, [], 2);
M = P(idx,:);



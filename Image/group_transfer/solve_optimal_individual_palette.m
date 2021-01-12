function [min_label_com, min_obj_func] = solve_optimal_individual_palette(P, W, M, lambda, gamma, eta, init_min)
m = size(M,1);
% brute-force all possible cases
min_obj_func = init_min;
D1 = pdist2(P, P);
D2 = psub2(P,M);
D3 = [zeros(size(P,1), 1), repmat(W, [1,m]).*(pdist2(P, M).^2)];

n = size(P,1);
num_of_all_combination = (m+1)^n;
label = zeros(n,1);
min_label_com = label;
label(n) = -1;
for idx = 1:num_of_all_combination    
    label(n) = label(n) + 1;
    curId = n;
    while(label(curId) > m)
        label(curId) = 0;
        curId = curId - 1;
        label(curId) = label(curId) + 1;
    end
    term4 = sum((label==0).*W);
    val = eta*term4;
    if val >= min_obj_func
        continue;
    end
    for i = 1:n
        val = val + D3(i,label(i)+1); % the first term
    end
    % cut down the unsolution
    if val >= min_obj_func
        continue;
    end
    term2 = 0;
    term3 = 0;
    for ii = 1:n-1
        for jj = ii+1:n
                term2 = term2 + D2(ii,jj,label(ii)+1,label(jj)+1);
            if label(ii) == label(jj) && label(ii) > 0
                term3 = term3 + D1(ii,jj);
            end
        end
    end

    val = val + lambda*term2 + gamma*term3;
    if val < min_obj_func
        min_obj_func = val;
        min_label_com = label;
    end

end

function D = psub2(P, M)
p = size(P,1);
m = size(M,1);
D = zeros(p, p, m + 1, m + 1);
for i1 = 1:p-1
    for i2 = i1+1:p
        for i3 = 1:m
            for i4 = 1:m
                D(i1,i2,i3+1,i4+1) = sum((P(i1,:) - P(i2,:) - M(i3,:) + M(i4,:)).^2);
            end
        end
    end
end


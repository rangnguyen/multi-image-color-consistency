function [Iout, Pout] = Lab_transfer(Iin, C1, C2, LUT1, LUT2)
% Convert from RGB to Lab

[m, n, b] = size(Iin);
Iin = reshape(Iin, [m*n, b]);
Pout = C2;

Iin = RGB2Lab(Iin, LUT1);
C1 = RGB2Lab(C1, LUT1);
C2 = RGB2Lab(C2, LUT1);


Iout = Iin;
% Iout(:,1) = L_transfer(Iin(:,1), C1(:,1), C2(:,1));
Iout(:,2:3) = ab_transfer(Iin(:,2:3), C1(:,2:3), C2(:,2:3));


% Convert from Lab to RGB
Iout = Lab2RGB(Iout, LUT2);
Iout = reshape(Iout, [m,n,b]);



function Lout = L_transfer(Lin, L1, L2)
L1 = sort(L1, 'descend');
L2 = sort(L2, 'descend');

% Insert the white and black level
L1 = [100; L1; 0];
L2 = [100; L2; 0];

n = size(Lin, 1);
Lout = Lin;
for i = 1:n
    j = 2;
    while(L1(j) > Lin(i))
        j = j + 1;
    end
    d1 = L1(j-1) - Lin(i);
    d2 = Lin(i) - L1(j);
    d = d1 + d2;
    Lout(i) = d2/d * L2(j-1) + d1/d * L2(j);
end

function [C1, C2, ids] = Lab_matching_distance(C1, C2)

K = size(C1, 1);
P = perms(1:K);
M = size(P, 1);
D = zeros(M, 1);

for i = 1:M
    for j = 1:K
        D(i) = D(i) + sum((C1(j,:)-C2(P(i,j),:)).^2);
    end
end
[~, id] = min(D);
C2 = C2(P(id, :),:);
ids = P(id,:);


function [C1, C2, ids] = Lab_matching_combine(C1, C2, H1, H2)

K = size(C1, 1);
P = perms(1:K);
M = size(P, 1);
D = zeros(M, 1);

for i = 1:M
    for j = 1:K
        D(i) = D(i) + max(H1(j), H2(P(i,j)))*sum((C1(j,:)-C2(P(i,j),:)).^2);
    end
end
[~, id] = min(D);
C2 = C2(P(id, :),:);
ids = P(id,:);
        
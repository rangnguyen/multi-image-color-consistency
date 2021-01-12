function Io = ab_transfer(Ii, C1, C2)

[m, b] = size(Ii);
% remove close color
k = size(C1,1);
eps = 0.0001;

W = zeros(m, k);
for i = 1:k
    D = zeros(m, 1);
    for j = 1:b
        D = D + (Ii(:,j) - C1(i,j)).^2;
    end
    W(:,i) = 1./(D + eps);
end

Io = zeros(size(Ii));
% Normalize W
sumW = sum(W, 2);
for j = 1:k
    W(:,j) = W(:,j) ./ sumW;
end

for i = 1:k
    for j = 1:b
        Io(:,j) = Io(:,j) + W(:,i) .* (Ii(:,j) + C2(i,j) - C1(i,j));
    end
end


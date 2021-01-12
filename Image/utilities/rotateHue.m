function P = rotateHue(P, val)
[m,n, b] = size(P);
if(b > 1)
    P = reshape(P, [m*n, b]);
    P(P<0) = 0;
    P(P>1) = 1;
end
P = rgb2hsv(P);
H = P(:,1);
H = H + val;
H(H<0) = 1+H(H<0);
H(H>1) = H(H>1)-1;
P(:,1) = H;
P = hsv2rgb(P);
if (b > 1)
    P = reshape(P, [m, n, b]);
end
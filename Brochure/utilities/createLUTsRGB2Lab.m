function createLUTsRGB2Lab()

colorTransform = makecform('srgb2lab');
N2 = 256*256;
N3 = N2 * 256;
LUT1 = zeros(N3, 3);
[X, Y] = meshgrid(0:255, 0:255);
Z = [X(:), Y(:)];
for i = 1:256
    K = [(i-1)*ones(N2,1), Z]/255;
    T = int16(applycform(K, colorTransform));
    LUT1((i-1)*N2+1:i*N2,:) = T;
    disp(i);
end

colorTransform = makecform('lab2srgb');
N = 101*256*256;
LUT2 = zeros(N, 3);
[X, Y] = meshgrid(-127:128, -127:128);
Z = [X(:), Y(:)];
for i = 1:101
    K = [(i-1)*ones(N2,1), Z];
    T = uint8(applycform(K, colorTransform)*255);
    LUT2((i-1)*N2+1:i*N2,:) = T;
    disp(i);
end


save('data', 'LUT1', 'LUT2');

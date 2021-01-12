function I = insertImagesIntoBrochure(I, J, name)

mask = getMaskFromBrochure(name);
m = size(mask,1);
for i = 1:m
    K = J{i};
    hm = mask(i,2) - mask(i,1)+1;
    wm = mask(i,4) - mask(i,3)+1;
    [hi, wi, ~] = size(K);
    if((hm ~= hi) || (wm ~= wi))
        K = imresize(K, [hm, wm]);
    end
    I(mask(i,1):mask(i,2), mask(i,3):mask(i,4), :) = K;
end
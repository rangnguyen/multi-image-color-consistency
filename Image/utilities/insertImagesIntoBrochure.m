function I = insertImagesIntoBrochure(I, J, name)

mask = getMaskFromBrochure(name);
m = size(mask,1);
for i = 1:m
    I(mask(i,1):mask(i,2), mask(i,3):mask(i,4), :) = J{i};
end
function I = excludeImagesFromBrochure(I, name)

mask = getMaskFromBrochure(name);
m = size(mask,1);
expand = 5;
for i = 1:m
    is = max(1, (mask(i,1)-expand));
    ie = min(size(I,1), mask(i,2)+expand);
    js = max(1, (mask(i,3)-expand));
    je = min(size(I,2), mask(i,4)+expand);
    I(is:ie,js:je, :) = 0;
end
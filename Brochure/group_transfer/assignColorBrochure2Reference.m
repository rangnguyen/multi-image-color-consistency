function [Pout, Match] = assignColorBrochure2Reference(Pbro, Pref, flag)

Pout = Pref;
Match = zeros(size(Pref, 1), 1);

if flag == 0
    return;
end

Hsvbro = rgb2hsv(Pbro);
Hsvref = rgb2hsv(Pref);


[vals, ids] = sort(Hsvbro(:,2), 'descend');
pivot = [];
for i = 1:length(vals)
    if(vals(i) > 0.6)
        j = 1;
        while(j <= length(pivot))
            if HueDistance(Hsvbro(ids(i),1), Hsvbro(pivot(j),1)) < 0.03
                break;
            end
            j = j + 1;
        end
        if (j > length(pivot))
            pivot = [pivot; ids(i)];
        end
    end
end

if isempty(pivot)
    return;
end
    
d = zeros(size(pivot));
for i = 1:size(Pref, 1)
    if Hsvref(i,2) > 0.25
        for j = 1:length(pivot)
            d(j) = HueDistance(Hsvbro(pivot(j),1), Hsvref(i,1));
        end
        [val, id] = min(d);
        if val < 0.1
            Pout(i,:) = Pbro(pivot(id),:);
            Match(i) = pivot(id);
        end
    end
end

function d = HueDistance(a, b)
d = abs(a-b);
if d > 0.5
    d = 1-d;
end
function Pout = assignColorReference2InputManual(Pin, Match, Prefbro)
Pout = Pin;
% Minimum the 2d distance
for t = 1:length(Match)   
    mask = Match{t} > 0;
    Pout{t}(mask,:) = Prefbro(Match{t}(mask),:);    
end

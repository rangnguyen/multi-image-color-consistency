function P = updatePalette(P, Match, NP)
t = Match > 0;
P(t,:) = NP(Match(t),:);
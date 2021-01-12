function unmatch = checkUnMatchBrochurePalette(MatchBro, MatchPre)

unmatch = ones(length(MatchPre), 1);
for i = 1:length(MatchPre)
    flag = 0;
    for j = 1:length(MatchPre{i})
        if MatchPre{i}(j) > 0 && MatchBro(MatchPre{i}(j)) > 0
            flag = 1;
            break;
        end
    end
    if flag
        unmatch(i) = 0;
    end
end

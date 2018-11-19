function idx = dispatch(nEval, nCore)
nEvals = zeros(nCore,1);
tmp = floor(nEval / nCore);
if tmp > 0
    nEvals = nEvals + tmp;
end

r  = rem(nEval, nCore);
if r > 0
    nEvals(1:r) = nEvals(1:r) + 1;
end
idx = cell(length(nEvals),1);
offset = 0;
for i = 1:length(nEvals)
    idx{i} = offset + (1:nEvals(i));
    offset = idx{i}(end);
end
end
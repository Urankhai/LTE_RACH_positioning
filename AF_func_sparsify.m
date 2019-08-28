function p_new = AF_func_sparsify(p, th)
p_new = p;
for k = 1:length(p)
    if abs(p(k)) < th
        p_new(k) = 0;
    else
        p_new(k) = p(k)*(abs(p(k)) - th)/abs(p(k));
    end
end

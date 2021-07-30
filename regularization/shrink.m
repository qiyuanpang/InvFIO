function ans1 = shrink(x, r)
    absx = abs(x);
    ans1 = (x./absx) .* max(absx - r, 0);
end
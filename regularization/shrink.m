function ans = shrink(x, r)
    absx = abs(x);
    ans = (x./absx) .* max(absx - r, 0);
end
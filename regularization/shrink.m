function ans1 = shrink(x, r)
    absx = abs(x);
    ans1 = (x./absx) .* max(absx - r, 0);
    ans1(isnan(ans1)) = 0;
end
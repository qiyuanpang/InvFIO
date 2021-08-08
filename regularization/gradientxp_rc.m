function v = gradientxp_rc(f, h)
    v = (4*gradientxp(f,h)-gradientxp_2h(f,h))/3;
end

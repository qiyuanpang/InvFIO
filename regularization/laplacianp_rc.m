function v = laplacianp_rc(f,h)
    v = (4*laplacianp(f,h)-laplacianp_2h(f,h))/3;
end

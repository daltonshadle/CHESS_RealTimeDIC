function eps = FitStrain(coord,disp)
    p=polyfit(coord,disp,1);
    eps=p(1);
end
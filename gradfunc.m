function ypr = gradfunc(x)
    syms x1 x2
    for i = 1:length(x1)
        for j = 1:length(x2)
            F(i, j) = (x1(i) - 1)^2 + (x2(j) - 1)^2 - x1(i) * x2(j);
        end
    end
    g=gradient(F);
    ypr=double(subs(g,[x1 x2],[x(1) x(2)]));
end
function H = hessianfunc(x)
    syms x1 x2
    for i = 1:length(x1)
        for j = 1:length(x2)
            F(i, j) = (x1(i) - 1)^2 + (x2(j) - 1)^2 - x1(i) * x2(j);
        end
    end
    H = hessian(F, [x1, x2]);
    H = double(subs(H, [x1, x2], [x(1), x(2)]));
end

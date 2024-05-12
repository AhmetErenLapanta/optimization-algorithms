function y = func(x)
    x1 = x(1);
    x2 = x(2);
    for i = 1:length(x1)
        for j = 1:length(x2)
            y = (x1(i) - 1)^2 + (x2(j) - 1)^2 - x1(i) * x2(j);
        end
    end
end
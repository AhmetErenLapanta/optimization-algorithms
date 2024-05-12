clear all
close all
clc

n = 3;
x1= -n^2:0.01:n^2;
x2= -n^2:0.01:n^2;
F = zeros(length(x1), length(x2));
for i = 1:length(x1)
    for j = 1:length(x2)
        F(i, j) = (x1(i) - 1)^2 + (x2(j) - 1)^2 - x1(i) * x2(j);
    end
end
realFMin = min(min(F));
mesh(x1,x2,F)

figure
contourf(x1,x2,F)
hold on

num_points = 3;
initial_points = randn(2, num_points); 
colors = hsv(num_points); 
epsilon = 10^(-4);

%% Newton-Raphson
fprintf('Newton-Raphson Algorithm\n');

for i = 1:num_points
    x = initial_points(:, i); % Assign initial point
    color = colors(i, :); % Assign color
    
    tic
    fprintf('k=1, x1=%f, x2=%f, f(x)=%f\n',x(1),x(2),func(x))
    plot(x(1),x(2),'r.')
    x_next=x-inv(hessianfunc(x))*gradfunc(x);
    fprintf('k=2, x1=%f, x2=%f, f(x)=%f, error=%f\n',x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
    plot(x_next(1),x_next(2),'r*')
    k = 3;
    while(abs(func(x_next)-func(x)) > epsilon)
        x = x_next;
        x_next = x - inv(hessianfunc(x)) * gradfunc(x);
        fprintf('k=%d, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n', k, x_next(1), x_next(2), func(x_next), abs(func(x_next)-func(x)));
        plot(x_next(1), x_next(2), 'r*');
        k = k + 1;
    end
    toc
    title('Newton-Raphson Algorithm')
    set(gca,'fontsize',35)
    set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);
end


%% Hestenes-Stiefel Algorithm
figure
contourf(x1,x2,F)
hold on

fprintf('Hestenes-Stiefel Algorithm\n');

for i = 1:num_points
    x = initial_points(:, i);
    color = colors(i, :);
    
    tic
    fprintf('k=1, x1=%f, x2=%f, f(x)=%f\n',x(1),x(2),func(x))
    plot(x(1),x(2),'r.')
    % TO BE WRITTEN...
    g=gradfunc(x);
    d=-g;
    
    % alpha-argmin procedure
    alpha=0:0.01:1;
    funcalpha=zeros(length(alpha),1);
    
    for i=1:length(alpha)
        funcalpha(i)=func(x+alpha(i)*d);
    end
    [val,ind]=min(funcalpha);
    alpha=alpha(ind);
    % end of alpha-argmin procedure
    x_next=x+alpha*d;
    g_next=gradfunc(x_next);
    beta=(g_next'*(g_next-g))/(d'*(g_next-g));
    d_next=-g_next+beta*d;
    
    fprintf('k=2, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
    plot(x_next(1),x_next(2),'r*')
    k=3;
    
    while(abs(func(x_next)-func(x))>epsilon)
        x=x_next;
        g=g_next;
        d=d_next;
    
        % alpha-argmin procedure
        alpha=0:0.01:1;
        funcalpha=zeros(length(alpha),1);
        for i=1:length(alpha)
            funcalpha(i)=func(x+alpha(i)*d);
        end
        [val,ind]=min(funcalpha);
        alpha=alpha(ind);
        % end of alpha-argmin procedure
    
        x_next=x+alpha*d;
        g_next=gradfunc(x_next);
        beta=(g_next'*(g_next-g))/(d'*(g_next-g));
        d_next=-g_next+beta*d;
    
        fprintf('k=%d, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',k,x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
        
        plot(x_next(1),x_next(2),'r*')
        k=k+1;
    end
    toc
    title('Hestenes-Stiefel Algorithm')
    set(gca,'fontsize',35)
    set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);
end


%% Polak-Ribiere Algorithm
figure
contourf(x1,x2,F)
hold on

fprintf('Polak-Ribiere Algorithm\n');
for i = 1:num_points
    x = initial_points(:, i); % Assign initial point
    color = colors(i, :); % Assign color

    tic
    fprintf('k=1, x1=%f, x2=%f, f(x)=%f\n',x(1),x(2),func(x))
    plot(x(1),x(2),'r.')
    g=gradfunc(x);
    d=-g;
    alpha=0:0.01:1;
    funcalpha=zeros(length(alpha),1);
    for i=1:length(alpha)
        funcalpha(i)=func(x+alpha(i)*d);
    end
    [val,ind]=min(funcalpha);
    alpha=alpha(ind);
    x_next=x+alpha*d;
    fprintf('k=2, x1=%f, x2=%f, f(x)=%f, error=%f\n',x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
    plot(x_next(1),x_next(2),'r*')
    k=3;
    
    while(abs(func(x_next)-func(x))>epsilon)
        x=x_next;
        g=gradfunc(x);
        beta=(g'*(g-gradfunc(x_next)))/(d'*(g-gradfunc(x_next)));
        d_next=-g+beta*d;
        alpha=0:0.01:1;
        funcalpha=zeros(length(alpha),1);
        for i=1:length(alpha)
            funcalpha(i)=func(x+alpha(i)*d_next);
        end
        [val,ind]=min(funcalpha);
        alpha=alpha(ind);
        x_next=x+alpha*d_next;
        fprintf('k=%d, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',k,x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
        plot(x_next(1),x_next(2),'r*')
        k=k+1;
    end
    toc
    title('Polak-Ribiere Algorithm')
    set(gca,'fontsize',35)
    set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);
end

%% Fletcher-Reeves Algorithm
figure
contourf(x1,x2,F)
hold on
fprintf('Fletcher-Reeves Algorithm\n');

for i = 1:num_points
    x = initial_points(:, i); 
    color = colors(i, :);
    tic
    fprintf('k=1, x1=%f, x2=%f, f(x)=%f\n',x(1),x(2),func(x))
    plot(x(1),x(2),'r.')
    g=gradfunc(x);
    d=-g;
    alpha=0:0.01:1;
    funcalpha=zeros(length(alpha),1);
    for i=1:length(alpha)
        funcalpha(i)=func(x+alpha(i)*d);
    end
    [val,ind]=min(funcalpha);
    alpha=alpha(ind);
    x_next=x+alpha*d;
    fprintf('k=2, x1=%f, x2=%f, f(x)=%f, error=%f\n',x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
    plot(x_next(1),x_next(2),'r*')
    k=3;
    
    while(abs(func(x_next)-func(x))>epsilon)
        x=x_next;
        g=gradfunc(x);
        beta=(norm(g)^2)/(norm(gradfunc(x))^2);
        d_next=-g+beta*d;
        alpha=0:0.01:1;
        funcalpha=zeros(length(alpha),1);
        for i=1:length(alpha)
            funcalpha(i)=func(x+alpha(i)*d_next);
        end
        [val,ind]=min(funcalpha);
        alpha=alpha(ind);
        x_next=x+alpha*d_next;
        fprintf('k=%d, x1=%f, x2=%f, f(x)=%f, abs. error=%f\n',k,x_next(1),x_next(2),func(x_next),abs(func(x_next)-func(x)))
        plot(x_next(1),x_next(2),'r*')
        k=k+1;
    end
    toc
    title('Fletcher-Reeves Algorithm')
    set(gca,'fontsize',35)
    set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 2);
end
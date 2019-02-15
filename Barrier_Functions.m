% Начало работы с SeDuMi
checkDependency('spotless');
prog = spotsosprog;

solver = @spot_sedumi;
solver_name = 'sedumi';

% x = [x1;x2]
x1 = msspoly('x');
x2 = msspoly('y');
x = [x1;x2];
prog = prog.withIndeterminate(x);

% Динамика
f = [x2; -x1 + x1^3 - x2];

% Создание барьерной функции
d = 4;
[prog,B] = prog.newFreePoly(monomials(x,0:d));

% Нахождение производной от В
Bdot = diff(B,x)*f;

% Данные полученные при решение предыдущих частей задания
g_X0 = 0.5^2 - (x1 - 1.5)^2 - x2^2;
g_Xu = 0.4^2 - (x1+1)^2 - (x2+1)^2;

% Нахождение постоянных
[prog,L1] = prog.newSOSPoly(monomials(x,0:d-2));
prog = prog.withSOS(-B - L1*g_X0 - 1e-4);
                                       
[prog,L2] = prog.newSOSPoly(monomials(x,0:d-2));
prog = prog.withSOS(B - L2*g_Xu);

prog = prog.withSOS(-Bdot);

% Решение задачи суммы квадратов
options = spot_sdp_default_options();
options.verbose = 1;
sol = prog.minimize(0,solver,options);

% Проверка на наличие решения
if ~sol.isPrimalFeasible
    error('Нет решения');
end

if strcmp(solver_name,'sedumi')
    if sol.info.solverInfo.feasratio < 0
        error('Нет решения');
    end
end

% Проверка для В
xprint = msspoly('x',2);
disp(' ');
disp('Barrier function:');
B_sol = clean(subs(sol.eval(B),x,xprint),1e-5)

% Построение графиков
figure
[X,Y] = meshgrid(linspace(-3,3,100),linspace(-3,3,100));

gs_X0 = msubs(g_X0,x,[X(:),Y(:)]');
contour(X,Y,reshape(double(gs_X0),100,100),[0 0],'LineWidth',3);
hold on

gs_Xu = msubs(g_Xu,x,[X(:),Y(:)]');
contour(X,Y,reshape(double(gs_Xu),100,100),[0 0],'r','LineWidth',3)

gs_B = msubs(sol.eval(B),x,[X(:),Y(:)]');
contour(X,Y,reshape(double(gs_B),100,100),[0 0],'b','LineWidth',3)

% Построение векторного поля
[X,Y] = meshgrid(linspace(-3,3,50),linspace(-3,3,50));
x1dots = reshape(double(msubs(f(1),x,[X(:),Y(:)]')),50,50);
x2dots = reshape(double(msubs(f(2),x,[X(:),Y(:)]')),50,50);
x1dots = 0.1*x1dots./(sqrt(x1dots.^2 + x2dots.^2));
x2dots = 0.1*x2dots./(sqrt(x1dots.^2 + x2dots.^2));
quiver(X,Y,x1dots,x2dots,'AutoScale','off','Color','k');

% Легенда
title('Barrier functions')
xlabel('x_1');
ylabel('x_2');

legend('X_0','X_u','0 level-set of B(x)');
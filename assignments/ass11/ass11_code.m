format long;
close all;
clear all;
clc

tol=1e-13;
s=1e-6;
max_counter=1000;
xmu = 0.1;
ysign=-1;
delta=1e-3;
h=1e-3;

struct = load('data.mat');
data = struct.data;
[L3,C3,eigL3,veigL3]=RTBP_eq_points(xmu,tol);
x0=L3(1);

T_i=data(100,3):0.001:round(data(end,3),3);
size_T=size(T_i);
x0=1.123959614214410;
y0_prime=-0.175736293599871;
options = odeset('RelTol',1e-13,'AbsTol',1e-13);  % ODE solver options for accuracy
orbits_data=zeros(size_T(2),5);
counter=1;

for T=T_i
    [x,y_prime]=newton_variational(x0,y0_prime,xmu,T,tol);
    [~,X]=ode45(@(t,y) f_variational_eq(t,y,xmu,1), [0,T], [x 0 0 y_prime,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1], options);
    r1 = sqrt((x-xmu)^2);
    r2 = sqrt((x-xmu+1)^2);
    C=2*Omega([x,0],xmu,r1,r2) - y_prime^2;
    k=X(end,5)+X(end,10)+X(end,15)+X(end,20)-2;
    orbits_data(counter,1)=C;
    orbits_data(counter,2)=x;
    orbits_data(counter,3)=y_prime;
    orbits_data(counter,4)=T;
    orbits_data(counter,5)=k;
    counter=counter+1;
    % Imprime los valores actuales de T y x
    fprintf('T: %.6f, x: %.6f\n', T, x);

    x0=x;
    y0_prime=y_prime;
end

% Extrae los datos necesarios desde la matriz orbits_data
x_values = orbits_data(:, 2); % x está en la segunda columna de orbits_data
T_values = orbits_data(:, 4); % T está en la cuarta columna de orbits_data
C_values = orbits_data(:, 1); % C está en la primera columna de orbits_data
k_values = orbits_data(:, 5); % k está en la quinta columna de orbits_data

% Grafica T(x)
figure;
plot(x_values, T_values, 'g-o', 'LineWidth', 1.5); % Línea verde con marcadores
xlabel('x', 'FontSize', 14); % Etiqueta del eje x más grande
ylabel('T(x)', 'FontSize', 14); % Etiqueta del eje y más grande
title('Gráfica de T(x)', 'FontSize', 16); % Título con fuente más grande
grid on;
saveas(gcf, 'Grafica_Tx.png'); % Guarda como PNG

% Grafica C(x)
figure;
plot(x_values, C_values, 'g-o', 'LineWidth', 1.5); % Línea verde con marcadores
xlabel('x', 'FontSize', 14); % Etiqueta del eje x más grande
ylabel('C(x)', 'FontSize', 14); % Etiqueta del eje y más grande
title('Gráfica de C(x)', 'FontSize', 16); % Título con fuente más grande
grid on;
saveas(gcf, 'Grafica_Cx.png'); % Guarda como PNG

% Grafica k(x)
figure;
plot(x_values, k_values, 'g-o', 'LineWidth', 1.5); % Línea verde con marcadores
xlabel('x', 'FontSize', 14); % Etiqueta del eje x más grande
ylabel('k(x)', 'FontSize', 14); % Etiqueta del eje y más grande
title('Gráfica de k(x)', 'FontSize', 16); % Título con fuente más grande
grid on;
saveas(gcf, 'Grafica_kx.png'); % Guarda como PNG

%%
xmu = 0.01;
T_i=3.114802556760205*2:0.0001:3.114802556760205*2+0.0001*1000;
size_T=size(T_i);
x0=1.033366313746765;
y0_prime=-.05849376854515592;
options = odeset('RelTol',1e-13,'AbsTol',1e-13);  % ODE solver options for accuracy
orbits_data=zeros(size_T(2),5);
counter=1;

for T=T_i
    [x,y_prime]=newton_variational(x0,y0_prime,xmu,T,tol);
    [~,X]=ode45(@(t,y) f_variational_eq(t,y,xmu,1), [0,T], [x 0 0 y_prime,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1], options);
    r1 = sqrt((x-xmu)^2);
    r2 = sqrt((x-xmu+1)^2);
    C=2*Omega([x,0],xmu,r1,r2) - y_prime^2;
    k=X(end,5)+X(end,10)+X(end,15)+X(end,20)-2;
    orbits_data(counter,1)=C;
    orbits_data(counter,2)=x;
    orbits_data(counter,3)=y_prime;
    orbits_data(counter,4)=T;
    orbits_data(counter,5)=k;
    counter=counter+1;
    % Imprime los valores actuales de T y x
    fprintf('T: %.6f, x: %.6f\n', T, x);

    x0=x;
    y0_prime=y_prime;
end
%%
% Extrae los datos necesarios desde la matriz orbits_data
x_values = orbits_data(1:407, 2); % x está en la segunda columna de orbits_data
T_values = orbits_data(1:407, 4); % T está en la cuarta columna de orbits_data
C_values = orbits_data(1:407, 1); % C está en la primera columna de orbits_data
k_values = orbits_data(1:407, 5); % k está en la quinta columna de orbits_data

% Grafica T(x)
figure;
plot(x_values, T_values, 'g-o', 'LineWidth', 1.5); % Línea verde con marcadores
xlabel('x', 'FontSize', 14); % Etiqueta del eje x más grande
ylabel('T(x)', 'FontSize', 14); % Etiqueta del eje y más grande
title('Gráfica de T(x)', 'FontSize', 16); % Título con fuente más grande
grid on;
saveas(gcf, 'Grafica_Tx.png'); % Guarda como PNG

% Grafica C(x)
figure;
plot(x_values, C_values, 'g-o', 'LineWidth', 1.5); % Línea verde con marcadores
xlabel('x', 'FontSize', 14); % Etiqueta del eje x más grande
ylabel('C(x)', 'FontSize', 14); % Etiqueta del eje y más grande
title('Gráfica de C(x)', 'FontSize', 16); % Título con fuente más grande
grid on;
saveas(gcf, 'Grafica_Cx.png'); % Guarda como PNG

% Grafica k(x)
figure;
plot(x_values, k_values, 'g-o', 'LineWidth', 1.5); % Línea verde con marcadores
xlabel('x', 'FontSize', 14); % Etiqueta del eje x más grande
ylabel('k(x)', 'FontSize', 14); % Etiqueta del eje y más grande
title('Gráfica de k(x)', 'FontSize', 16); % Título con fuente más grande
grid on;
saveas(gcf, 'Grafica_kx.png'); % Guarda como PNG

%% STATE 2
tol=1e-13;
s=1e-6;
max_counter=10000;
xmu = 0.01;
ysign=-1;
delta=1e-3;
h=1e-3;

T_i=(2*3.114802556760205):0.0001:2*3.125;
size_T=size(T_i);
x0=1.033366313746765;
y0_prime=-0.05849376854515592;
options = odeset('RelTol',1e-13,'AbsTol',1e-13);  % ODE solver options for accuracy
orbits_data=zeros(size_T(2),5);
counter=1;




el = 0;
for T = T_i
    % Imprimir valores antes de llamar a newton_variational
    fprintf('Antes de Newton:\n x0 = %.8f\n y0_prime = %.8f\n\n, T = %.6f\n\n', x0, y0_prime, T);
    
    [x_new, y_prime] = newton_variational(x0, y0_prime, xmu, T, tol);
    
    [~, X] = ode45(@(t, y) f_variational_eq(t, y, xmu, 1), [0, T], [x_new, 0, 0, y_prime, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1], options);
    
    r1 = sqrt((x_new - xmu)^2);
    r2 = sqrt((x_new - xmu + 1)^2);
    C = 2 * Omega([x_new, 0], xmu, r1, r2) - y_prime^2;
    k = X(end, 5) + X(end, 10) + X(end, 15) + X(end, 20) - 2;
    
    orbits_data(counter, 1) = C;
    orbits_data(counter, 2) = x_new;
    orbits_data(counter, 3) = y_prime;
    orbits_data(counter, 4) = T;
    orbits_data(counter, 5) = k;
    counter = counter + 1;

    % Incrementar el contador de iteraciones
    el = el + 1;

    % Actualizar valores
    x0 = x_new;
    y0_prime = y_prime;

    % Imprimir valores después de la actualización
    fprintf('Después de actualizar:\n x0 = %.8f\n y0_prime = %.8f\n\n', x0, y0_prime);
    
    % Salir del bucle si se cumplen las condiciones
    
end


%%
% Extrae los datos necesarios desde la matriz orbits_data
x_values = orbits_data(:, 2); % x está en la segunda columna de orbits_data
T_values = orbits_data(:, 4); % T está en la cuarta columna de orbits_data
C_values = orbits_data(:, 1); % C está en la primera columna de orbits_data
k_values = orbits_data(:, 5); % k está en la quinta columna de orbits_data

% Grafica T(x)
figure;
plot(x_values, T_values, 'm-o', 'LineWidth', 1.5); % Línea magenta con marcadores
xlabel('x', 'FontSize', 14); % Etiqueta del eje x más grande
ylabel('T(x)', 'FontSize', 14); % Etiqueta del eje y más grande
title('Gráfica de T(x)', 'FontSize', 16); % Título con fuente más grande
grid on;


% Grafica C(x)
figure;
plot(x_values, C_values, 'm-o', 'LineWidth', 1.5); % Línea magenta con marcadores
xlabel('x', 'FontSize', 14); % Etiqueta del eje x más grande
ylabel('C(x)', 'FontSize', 14); % Etiqueta del eje y más grande
title('Gráfica de C(x)', 'FontSize', 16); % Título con fuente más grande
grid on;


% Grafica k(x)
figure;
plot(x_values, k_values, 'm-o', 'LineWidth', 1.5); % Línea magenta con marcadores
xlabel('x', 'FontSize', 14); % Etiqueta del eje x más grande
ylabel('k(x)', 'FontSize', 14); % Etiqueta del eje y más grande
title('Gráfica de k(x)', 'FontSize', 16); % Título con fuente más grande
grid on;


%% PART 2


T_i = 2 * 3.114802556760205;
xmu = 0.01;
h = 1e-3;
ysign = -1;
delta=1e-3;
tol=1e-13;

% INITIAL CONDITIONS
T = 6.22960511352041;
x0 = 1.0333662809371982;
y0_prime = -0.05849370327223067;
options = odeset('RelTol',1e-13,'AbsTol',1e-13);  
counter=1;
delta=0.0001:0.0001:0.001;
size_delta=size(delta);
matrix=zeros(size_delta(2),8);

% Imprimir valores
fprintf('T = %.15f\n', T);
fprintf('x0 = %.15f\n', x0);
fprintf('y0_prime = %.15f\n', y0_prime);

for deltaS = delta
    [x,y_prime]=newton_variational(x0,y0_prime,xmu,T,tol);
    [~,X]=ode45(@(t,y) f_variational_eq(t,y,xmu,1), [0,T], [x 0 0 y_prime,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1], options);
    DG = [X(end,5)      X(end,6)  X(end,7)  X(end,8);
            X(end,9)  X(end,10) X(end,11) X(end,12);
            X(end,13) X(end,14) X(end,15) X(end,16);
            X(end,17) X(end,18) X(end,19) X(end,20)];
    
    
    [veps,vaps]=eig(DG);
    tol2=1e-3;
    unstable_vap=0;
    unstable_vep=0;
    stable_vap=0;
    stable_vep=0;
    for i=1:4
        if real(vaps(i,i))>1+tol2
            unstable_vap=real(vaps(i,i));
            unstable_vep=sign(veps(2,i))*real(veps(:,i));
        end
        if real(vaps(i,i))<1-tol2
            stable_vap=real(vaps(i,i));
            stable_vep=sign(veps(2,i))*real(veps(:,i));
        end
    end
    s=1e-6+deltaS;
    x0_unstable_plus=[x0 0 0 y0_prime]'+s*unstable_vep;
    x0_unstable_minus=[x0 0 0 y0_prime]'-s*unstable_vep;
    x0_stable_plus=[x0 0 0 y0_prime]'+s*stable_vep;
    x0_stable_minus=[x0 0 0 y0_prime]'-s*stable_vep;
    
    [cross_points1,total_time1,cross_times1,orbit1]=poincare(1,x0_unstable_plus',xmu,1,h,tol);
    [cross_points2,total_time2,cross_times2,orbit2]=poincare(1,x0_unstable_minus',xmu,1,h,tol);
    [cross_points3,total_time3,cross_times3,orbit3]=poincare(1,x0_stable_plus',xmu,-1,h,tol);
    [cross_points4,total_time4,cross_times4,orbit4]=poincare(1,x0_stable_minus',xmu,-1,h,tol);

    figure(1)
    hold on;
    plot(orbit1(1,:),orbit1(2,:),"red")
    legend('Wu+')
    hold off;
    figure(2)
    hold on;
    plot(orbit2(1,:),orbit2(2,:),"blue")
    legend('Wu-')
    hold off;
    figure(3)
    hold on;
    plot(orbit3(1,:),orbit3(2,:),"green")
    legend('Ws+')
    hold off;
    figure(4)
    hold on;
    plot(orbit4(1,:),orbit4(2,:),"magenta")
    legend('Ws-')
    hold off;

    figure(5)
    hold on;
    plot(orbit1(1,:),orbit1(2,:),"red")
    hold on;
    plot(orbit2(1,:),orbit2(2,:),"blue")
    hold on;
    plot(orbit3(1,:),orbit3(2,:),"green")
    hold on;
    plot(orbit4(1,:),orbit4(2,:),"magenta")
    legend('Wu+','Wu-','Ws+','Ws-')
    hold off;
    
    matrix(counter,:)=[orbit1(2,end),orbit1(3,end),orbit2(2,end),orbit2(3,end),orbit3(2,end),orbit3(3,end),orbit4(2,end),orbit4(3,end)];
    counter=counter+1;
end
    figure(6)
    hold on
    plot(matrix(:,1),matrix(:,2),"red")
    hold on
    plot(matrix(:,3),matrix(:,4),"blue")
    hold on;
    plot(matrix(:,5),matrix(:,6),"green")
    hold on;
    plot(matrix(:,7),matrix(:,8),"magenta")
    legend('Wu+','Wu-','Ws+','Ws-')
    xlabel('y')
    ylabel('x''')
    %%
    {
    figure(7)
    plot(matrix(:,1),matrix(:,2),"red")
    legend('Wu+')
    xlabel('y')
    ylabel('x''')
    }

%%
T_i = 2 * 3.114802556760205;
xmu = 0.01;
h = 1e-3;
ysign = -1;
delta=1e-3;
tol=1e-13;

% INITIAL CONDITIONS
T = 6.22960511352041;
N = 1e-3; % Step size for t
time_points = linspace(0, T, floor(T / N) + 1); % Creates the interval

x0 = 1.0333662809371982;
y0_prime = -0.05849370327223067;
options = odeset('RelTol',1e-13,'AbsTol',1e-13);  
counter=1;
% delta=0.0001:0.0001:0.001;
% size_delta=size(delta);
% Parameters (Ensure x0, y0_prime, xmu, T, and tol are defined)

[x,y_prime]=newton_variational(x0,y0_prime,xmu,T,tol);
[~,X]=ode45(@(t,y) f_variational_eq(t,y,xmu,1), [0,T], [x 0 0 y_prime,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1], options);
DG = [X(end,5)      X(end,6)  X(end,7)  X(end,8);
        X(end,9)  X(end,10) X(end,11) X(end,12);
        X(end,13) X(end,14) X(end,15) X(end,16);
        X(end,17) X(end,18) X(end,19) X(end,20)];


[veps,vaps]=eig(DG);
tol2=1e-3;
unstable_vap=0;
unstable_vep=0;
stable_vap=0;
stable_vep=0;
for i=1:4
    if real(vaps(i,i))>1+tol2
        unstable_vap=real(vaps(i,i));
        unstable_vep=-sign(veps(2,i))*real(veps(:,i));
    end
    if real(vaps(i,i))<1-tol2
        stable_vap=real(vaps(i,i));
        stable_vep=-sign(veps(2,i))*real(veps(:,i));
    end
end
% s=1e-6+deltaS;
% Print results
fprintf('Unstable Eigenvalue: %.15f\n', unstable_vap);
fprintf('Unstable Eigenvector:\n');
disp(unstable_vep);

fprintf('Stable Eigenvalue: %.15f\n', stable_vap);
fprintf('Stable Eigenvector:\n');
disp(stable_vep);







% Parameters
N_points = 30;  % Number of points along each direction
s = 1e-6;  % Scaling factor
x_initial = [x, 0, 0, y_prime]';  % Starting point of the periodic orbit
T = 6.22960511352041;  % Period of the orbit
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);  % ODE solver options

% Storage for endpoints
endpoints = zeros(N_points, 8);  % [Wu+ y, Wu+ x', Wu- y, Wu- x', Ws+ y, Ws+ x', Ws- y, Ws- x']

% Plot 1: All endpoints together
figure(1);
hold on;
title('Endpoints for All Manifolds');
xlabel('y'); ylabel('x''');

% Individual Figures for Orbits
figure(2); hold on; title('Orbits of Wu+ Branch'); xlabel('x'); ylabel('y');
figure(3); hold on; title('Orbits of Wu- Branch'); xlabel('x'); ylabel('y');
figure(4); hold on; title('Orbits of Ws+ Branch'); xlabel('x'); ylabel('y');
figure(5); hold on; title('Orbits of Ws- Branch'); xlabel('x'); ylabel('y');

% Combined Plot for All Manifolds
figure(6);
hold on;
legend_handles = [];  % To store plot handles for the legend

% Loop for Wu+ (positive unstable branch)
figure(2);  % For individual plot of Wu+
for j = 1:N_points
    factor = (j - 1) / N_points;  % j/N (from 0 to 1)
    p_j = x_initial + s * unstable_vap^factor * unstable_vep;  % Initial point along unstable eigenvector
    
    % Compute orbit for Wu+
    [cross_points, total_time, cross_times, orbit] = poincare(1, p_j', xmu, 1, h, tol);

    % Store endpoint
    endpoints(j, 1:2) = [orbit(2, end), orbit(3, end)];

    % Plot endpoints for Wu+
    figure(1);
    plot(orbit(2, end), orbit(3, end), 'ro');  % Red dots for Wu+ endpoints

    % Plot individual Wu+
    figure(2);
    plot(orbit(1, :), orbit(2, :), 'red');

    % Combined plot Wu+
    figure(6);
    if j == 1
        handle = plot(orbit(1, :), orbit(2, :), 'red', 'DisplayName', 'Wu+');
        legend_handles(1) = handle;  % Store handle for Wu+
    else
        plot(orbit(1, :), orbit(2, :), 'red');
    end
end

% Loop for Wu- (negative unstable branch)
figure(3);  % For individual plot of Wu-
for j = 1:N_points
    factor = (j - 1) / N_points;  % j/N (from 0 to 1)
    p_j = x_initial - s * unstable_vap^factor * unstable_vep;  % Initial point along unstable eigenvector

    % Compute orbit for Wu-
    [cross_points, total_time, cross_times, orbit] = poincare(1, p_j', xmu, 1, h, tol);

    % Store endpoint
    endpoints(j, 3:4) = [orbit(2, end), orbit(3, end)];

    % Plot endpoints for Wu-
    figure(1);
    plot(orbit(2, end), orbit(3, end), 'go');  % Blue dots for Wu- endpoints

    % Plot individual Wu-
    figure(3);
    plot(orbit(1, :), orbit(2, :), 'green');

    % Combined plot Wu-
    figure(6);
    if j == 1
        handle = plot(orbit(1, :), orbit(2, :), 'green', 'DisplayName', 'Wu-');
        legend_handles(2) = handle;  % Store handle for Wu-
    else
        plot(orbit(1, :), orbit(2, :), 'green');
    end
end

% Loop for Ws+ (positive stable branch)
figure(4);  % For individual plot of Ws+
for j = 1:N_points
    factor = (j - 1) / N_points;  % j/N (from 0 to 1)
    p_j = x_initial + s * stable_vap^(-factor) * stable_vep;  % Initial point along stable eigenvector

    % Compute orbit for Ws+ (backward in time)
    [cross_points, total_time, cross_times, orbit] = poincare(1, p_j', xmu, -1, h, tol);

    % Store endpoint
    endpoints(j, 5:6) = [orbit(2, end), orbit(3, end)];

    % Plot endpoints for Ws+
    figure(1);
    plot(orbit(2, end), orbit(3, end), 'bo');  % Green dots for Ws+ endpoints

    % Plot individual Ws+
    figure(4);
    plot(orbit(1, :), orbit(2, :), 'blue');

    % Combined plot Ws+
    figure(6);
    if j == 1
        handle = plot(orbit(1, :), orbit(2, :), 'blue', 'DisplayName', 'Ws+');
        legend_handles(3) = handle;  % Store handle for Ws+
    else
        plot(orbit(1, :), orbit(2, :), 'blue');
    end
end

% Loop for Ws- (negative stable branch)
figure(5);  % For individual plot of Ws-
for j = 1:N_points
    factor = (j - 1) / N_points;  % j/N (from 0 to 1)
    p_j = x_initial - s * stable_vap^(-factor) * stable_vep;  % Initial point along stable eigenvector

    % Compute orbit for Ws- (backward in time)
    [cross_points, total_time, cross_times, orbit] = poincare(1, p_j', xmu, -1, h, tol);

    % Store endpoint
    endpoints(j, 7:8) = [orbit(2, end), orbit(3, end)];

    % Plot endpoints for Ws-
    figure(1);
    plot(orbit(2, end), orbit(3, end), 'mo');  % Magenta dots for Ws- endpoints

    % Plot individual Ws-
    figure(5);
    plot(orbit(1, :), orbit(2, :), 'magenta');

    % Combined plot Ws-
    figure(6);
    if j == 1
        handle = plot(orbit(1, :), orbit(2, :), 'magenta', 'DisplayName', 'Ws-');
        legend_handles(4) = handle;  % Store handle for Ws-
    else
        plot(orbit(1, :), orbit(2, :), 'magenta');
    end
end

% Combined plot legend
figure(6);
legend(legend_handles, {'Wu+', 'Wu-', 'Ws+', 'Ws-'}, 'Location', 'best');
title('Combined Orbits for All Manifolds');
xlabel('x');
ylabel('y');
hold off;

% Plot 7: Continuous Wu+ in y-x' space
figure(7);
hold on;
title('Wu+ in y-x'' Space');
xlabel('y');
ylabel('x''');
% Directly append the first point for continuity
plot([endpoints(:, 1); endpoints(1, 1)], [endpoints(:, 2); endpoints(1, 2)], 'r-', 'LineWidth', 1.5);  % Red line
hold off;

% Plot 8: Continuous Wu- in y-x' space
figure(8);
hold on;
title('Wu- in y-x'' Space');
xlabel('y');
ylabel('x''');
plot([endpoints(:, 3); endpoints(1, 3)], [endpoints(:, 4); endpoints(1, 4)], 'b-', 'LineWidth', 1.5);  % Blue line
hold off;

% Plot 9: Continuous Ws+ in y-x' space
figure(9);
hold on;
title('Ws+ in y-x'' Space');
xlabel('y');
ylabel('x''');
plot([endpoints(:, 5); endpoints(1, 5)], [endpoints(:, 6); endpoints(1, 6)], 'g-', 'LineWidth', 1.5);  % Green line
hold off;

% Plot 10: Continuous Ws- in y-x' space
figure(10);
hold on;
title('Ws- in y-x'' Space');
xlabel('y');
ylabel('x''');
plot([endpoints(:, 7); endpoints(1, 7)], [endpoints(:, 8); endpoints(1, 8)], 'm-', 'LineWidth', 1.5);  % Magenta line
hold off;

% Plot 11: Combined Continuous Endpoints for All Manifolds
figure(11);
hold on;
title('Combined Endpoints for All Manifolds');
xlabel('y');
ylabel('x''');
plot([endpoints(:, 1); endpoints(1, 1)], [endpoints(:, 2); endpoints(1, 2)], 'r-', 'LineWidth', 1.5, 'DisplayName', 'Wu+');  % Red line
plot([endpoints(:, 3); endpoints(1, 3)], [endpoints(:, 4); endpoints(1, 4)], 'g-', 'LineWidth', 1.5, 'DisplayName', 'Wu-');  % Blue line
plot([endpoints(:, 5); endpoints(1, 5)], [endpoints(:, 6); endpoints(1, 6)], 'b-', 'LineWidth', 1.5, 'DisplayName', 'Ws+');  % Green line
plot([endpoints(:, 7); endpoints(1, 7)], [endpoints(:, 8); endpoints(1, 8)], 'm-', 'LineWidth', 1.5, 'DisplayName', 'Ws-');  % Magenta line
legend('Location', 'best');
hold off;



%%
% Loop for Wu+ (positive unstable branch)
Jacobi_constants = zeros(N_points, 1); % Pre-allocate storage for Jacobi constants

for j = 1:N_points
    factor = (j - 1) / N_points;  % j/N (from 0 to 1)
    p_j = x_initial + s * unstable_vap^factor * unstable_vep;  % Initial point along unstable eigenvector
    
    % Compute orbit for Wu+
    [~, ~, ~, orbit] = poincare(1, p_j', xmu, 1, h, tol);

    % Calculate Jacobi constant
    x = orbit(1, end);  % Position x
    y = orbit(2, end);  % Position y
    x_dot = orbit(3, end);  % Velocity in x
    y_dot = orbit(4, end);  % Velocity in y
    r1 = sqrt((x - xmu)^2 + y^2);
    r2 = sqrt((x - xmu + 1)^2 + y^2);
    Omega_value = 0.5 * (x^2 + y^2) + (1 - xmu) / r1 + xmu / r2;  % Potential energy
    Jacobi_constants(j) = 2 * Omega_value - (x_dot^2 + y_dot^2);
end

% Compute the derivative of the Jacobi constant
Jacobi_derivative = diff(Jacobi_constants);

% Plot the derivative
figure;
plot(1:N_points-1, Jacobi_derivative, 'b-', 'LineWidth', 1.5);
title('Derivative of the Jacobi Constant for Wu+');
xlabel('Orbit Index');
ylabel('dC_J / dIndex');
grid on;
% Set y-axis limits for a range around 10^-6
ylim([-1e-9, 1e-9]); % Adjust the y-axis to range [-5e-13, 5e-13]


%%
function [x0,y0_prime]=newton_variational(x0,y0_prime,xmu,T,tol)
    deltaT=0.0001;
    options = odeset('RelTol',1e-13,'AbsTol',1e-13);  % ODE solver options for accuracy
    [~,X]=ode45(@(t,y) f_variational_eq(t,y,xmu,1), [0,T/2], [x0 0 0 y0_prime,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1], options);
    x_prime=X(end,3);
    y=X(end,2);
    counter=1;
    while abs(x_prime)>tol || abs(y)>tol
        A=[X(end,9),X(end,12);X(end,13),X(end,16)]; 
        b=[-y;-x_prime];
        delta=A\b;
        x0=x0+delta(1);
        y0_prime=y0_prime+delta(2);
        [~,X]=ode45(@(t,y) f_variational_eq(t,y,xmu,1), [0,T/2], [x0 0 0 y0_prime,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1], options);
        x_prime=X(end,3);
        y=X(end,2);
    end
end


% Function defining the system of differential equations
function df = f_variational_eq(~,x,mu,dir)
  df = zeros(20,1);
  r1 = sqrt((x(1)-mu)^2 + x(2)^2);
  r2 = sqrt((x(1)-mu+1)^2 + x(2)^2);

  df(1) = x(3);
  df(2) = x(4);
  df(3) = 2*x(4) + x(1) - ((1-mu)*(x(1)-mu)/(r1^3)) - mu*(x(1)-mu+1)/(r2^3);
  df(4) = -2*x(3) + x(2)*(1 - (1-mu)/(r1^3) - mu/(r2^3));

  df(5)=x(13);
  df(6)=x(14);
  df(7)=x(15);
  df(8)=x(16);
  df(9)=x(17);
  df(10)=x(18);
  df(11)=x(19);
  df(12)=x(20);
  df(13)=Omegaxx([x(1) x(2)],mu)*x(5)+Omegaxy([x(1) x(2)],mu)*x(9)+2*x(17);
  df(14)=Omegaxx([x(1) x(2)],mu)*x(6)+Omegaxy([x(1) x(2)],mu)*x(10)+2*x(18);
  df(15)=Omegaxx([x(1) x(2)],mu)*x(7)+Omegaxy([x(1) x(2)],mu)*x(11)+2*x(19);
  df(16)=Omegaxx([x(1) x(2)],mu)*x(8)+Omegaxy([x(1) x(2)],mu)*x(12)+2*x(20);
  df(17)=Omegaxy([x(1) x(2)],mu)*x(5)+Omegayy([x(1) x(2)],mu)*x(9)-2*x(13);
  df(18)=Omegaxy([x(1) x(2)],mu)*x(6)+Omegayy([x(1) x(2)],mu)*x(10)-2*x(14);
  df(19)=Omegaxy([x(1) x(2)],mu)*x(7)+Omegayy([x(1) x(2)],mu)*x(11)-2*x(15);
  df(20)=Omegaxy([x(1) x(2)],mu)*x(8)+Omegayy([x(1) x(2)],mu)*x(12)-2*x(16);
  if dir == -1
      df=-df;
  end
end

function df = f(~,x,mu,dir)
  df = zeros(4,1);
  r1 = sqrt((x(1)-mu)^2 + x(2)^2);
  r2 = sqrt((x(1)-mu+1)^2 + x(2)^2);

  df(1) = x(3);
  df(2) = x(4);
  df(3) = 2*x(4) + x(1) - ((1-mu)*(x(1)-mu)/(r1^3)) - mu*(x(1)-mu+1)/(r2^3);
  df(4) = -2*x(3) + x(2)*(1 - (1-mu)/(r1^3) - mu/(r2^3));

  if dir == -1
      df=-df;
  end
end

function [L3,C3,eigL3,veigL3]=RTBP_eq_points(xmu,tol)
    %%%%%L3%%%%%
    xi=1-(7*xmu)/(12);
    aux=0;
    while(abs(xi-aux)>tol)
        aux=xi;
        xi=F3(xi,xmu);
    end
    L3=[xmu+xi,0,0,0];
    xL3=L3(1);
    C3 = 2*Omega([xL3 0],xmu,abs(xL3-xmu),abs(xL3-xmu+1));
    DG3 = zeros(4,4); DG3(1,3)=1; DG3(2,4)=1; DG3(3,4)=2; DG3(4,3)=-2;
    DG3(3,1) = Omegaxx([xL3 0],xmu);
    DG3(3,2) = Omegaxy([xL3 0],xmu);
    DG3(4,1) = Omegaxy([xL3 0],xmu);
    DG3(4,2) = Omegayy([xL3 0],xmu);
    [veigL3,eigL3] = eig(DG3);
    eigp3=zeros(4,1);
    for i = 1:4
        eigp3(i)=eigL3(i,i);
    end
    eigL3=eigp3;
end


function res = F(x,xmu,dy)
  [cross_points,~,~,~]=poincare(1,[x 0 0 dy],xmu,1,1e-2,1e-13);
  res = cross_points(3);
end


function dy=dy_i(x,xmu,C,ysign)
    r1 = sqrt((x-xmu)^2);
    r2 = sqrt((x-xmu+1)^2);
    dy= ysign*sqrt(2*Omega([x 0], xmu, r1, r2)-C);
end

function new=F1(old,xmu)
    new = ((xmu*(1-old)^2)/(3-2*xmu-old*(3-xmu-old)))^(1/3);
end
function new=F2(old,xmu)
    new = ((xmu*(1+old)^2)/(3-2*xmu+old*(3-xmu+old)))^(1/3);
end
function new=F3(old,xmu)
    new = (((1-xmu)*(1+old)^2)/(1+2*xmu+old*(2+xmu+old)))^(1/3);
end

function res = Omega(x,mu,r1,r2)
  res = (1/2)*(x(1)^2 + x(2)^2) + (1-mu)/r1 + mu/r2 + (1/2)*(mu*(1-mu));
end

function res = Omegaxx(x,mu)
  res =  (mu-1)/((mu-x(1))^2+x(2)^2)^(3/2)-mu/((x(1)-mu+1)^2+x(2)^2)^(3/2)-(3*(2*mu-2*x(1))^2*(mu-1))/(4*((mu-x(1))^2 ... 
      +x(2)^2)^(5/2))+(3*mu*(2*x(1)-2*mu+2)^2)/(4*((x(1)-mu+1)^2+x(2)^2)^(5/2))+1;
end

function res = Omegayy(x,mu)
  res = (mu-1)/((mu-x(1))^2+x(2)^2)^(3/2)-mu/((x(1)-mu+1)^2+x(2)^2)^(3/2)+(3*mu*x(2)^2)/((x(1)-mu+1)^2 ... 
      +x(2)^2)^(5/2)-(3*x(2)^2*(mu-1))/((mu-x(1))^2+x(2)^2)^(5/2)+1;
end

function res = Omegaxy(x,mu)
  res = (3*x(2)*(2*mu-2*x(1))*(mu-1))/(2*((mu-x(1))^2+x(2)^2)^(5/2)) ... 
      +(3*mu*x(2)*(2*x(1)-2*mu+2))/(2*((x(1)-mu+1)^2+x(2)^2)^(5/2));
end

function val = g(x)
%NEW POINCARÉ SECTION, THIS IS X = 0.2, g(x) must be zero
    val = x(1)-0.2;  
end


function [x1,x2]=varying(x,xmu,C,ysign,delta)
    x1=x;
    x2=x+delta;
    dy1=dy_i(x1,xmu,C,ysign);
    dy2=dy_i(x2,xmu,C,ysign);
    while F(x1,xmu,dy1)*F(x2,xmu,dy2)>=0
        x1=x2;
        x2=x2+delta;
        dy1=dy2;
        dy2=dy_i(x2,xmu,C,ysign);
        %F(x1,xmu,dy1)
    end
end

function [c,t]=bisection_method(x1,x2,xmu,C,ysign,tol,max_counter)
    dy1=dy_i(x1,xmu,C,ysign);
    dy2=dy_i(x2,xmu,C,ysign);
    bool=0;
    if F(x1,xmu,dy1)==0
        c=x1;
    bool=1;
    end
    if F(x2,xmu,dy2)==0 && bool==0
        c=x2;
    bool=1;
    end
    counter=0;
    while abs(x2-x1)>tol && counter<max_counter && bool==0
        c=(x2+x1)/2;
        dy1=dy_i(x1,xmu,C,ysign);
        dyc=dy_i(c,xmu,C,ysign);
        if F(c,xmu,dyc)==0 && bool==0
            bool=1;
        end
        if bool==0
        if F(x1,xmu,dy1)*F(c,xmu,dyc) <0
            x2=c;
        else
            x1=c;
        end
        end
        counter=counter+1;
    end
    c=(x2+x1)/2;
    [~,t,~,~]=poincare(2,[c 0 0 dy_i(c,xmu,C,ysign)],xmu,1,1e-4,1e-13);
end

function [cross_points,total_time,cross_times,orbit]=poincare(n_crossing,x0,xmu,idir,h,tol)

options = odeset('RelTol',1e-10,'AbsTol',1e-10);  % ODE solver options for accuracy

% Initialize arrays to store crossing points and times
cross_points = zeros(n_crossing, 4);  % Pre-allocate for speed: stores [position, velocity]
total_time = 0;                       % Initialize total time counter
cross_times = zeros(n_crossing, 1);   % Pre-allocate for speed: stores crossing times
orbit=[x0'];
for i = 1:n_crossing
    % For subsequent crossings, the initial condition is the previous crossing point
    if i ~= 1
        x0 = cross_points(i-1, :);  % Set x0 to the previous crossing point
    end
    val=6;
    found=0;
    while found==0
    % Solve the ODE using ode45 from t=0 to t=h with initial condition x0
    [t, x1] = ode45(@(t,y) f(t,y,xmu,idir), [total_time,total_time+val], x0, options);
    index=2;
    size_sol=size(x1);
    while index<=size_sol(1)-1 && found==0
        if g(x1(index,:))*g(x1(index+1,:))<0  
            found=1;
        else
        index=index+1;
        end
    end
    orbit=[orbit,x1(1:index,:)'];
    total_time=t(index);
    
    x0=x1(index,:);
    end
    % Approximate the crossing point
    approx = x1(index, :);
    % Refine the crossing point using Newton's method until the error is within tolerance
    counter=0;
    while abs(g(approx)) > tol
        counter=counter+1;
        % Compute correction delta using Newton's method
        delta = -g(approx) / approx(3);
        total_time = total_time + idir * delta;  % Update the total time
        if counter==300 || abs(approx(4))<1e-8
            total_time=-1;
            return
        end
        % Choose the direction for ODE integration based on the sign of delta
        if delta < 0
            [t, new_approx] = ode45(@(t,y) f(t,y,xmu,-1), [0, abs(delta)], approx, options);
        end
        if delta > 0
            [t, new_approx] = ode45(@(t,y) f(t,y,xmu,1), [0, abs(delta)], approx, options);
        end
        approx = new_approx(end, :);  % Update the approximation
        
    end
    % Store the refined crossing point and time
    cross_points(i, :) = approx;
    cross_times(i) = total_time;
    orbit=[orbit,approx'];
    end
end

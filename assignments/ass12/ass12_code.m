format long;
close all;
clear all;
clc

% Parameters
tol = 1e-13;
s = 1e-6;
max_counter = 1000;
xmu = 0.1;
ysign = -1;
delta = 1e-3;
h = 1e-3;
% forward
idir = 1;

% Initial conditions
x0 = 3.0;
C = 4.033150586970220;

% Compute initial y_prime
y_prime = dy_i(x0, xmu, C, ysign);
disp(y_prime);

% Poincaré options and computation
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
[cross_points, ~, ~, orbit1] = poincare(100, [x0, 0.0, 0.0, y_prime], xmu, idir, h, tol);

% Filter crossing points where y' < 0
% % even cross points
% cross_points_filtered = cross_points(2:2:end,:);
% Filtrar puntos de cruce donde y' < 0
cross_points_filtered = cross_points(cross_points(:,4) < 0, :);
% Figure 1: Plotting the orbit and filtered crossing points
figure;
plot(orbit1(1,:), orbit1(2,:), 'b-', 'LineWidth', 1.5); % Orbit in blue
hold on;
plot(cross_points_filtered(:,1), cross_points_filtered(:,2), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5); % Filtered points in red
xlabel('x');
ylabel('y');
title('Orbit and Filtered Poincaré Section Crossings (y'' < 0)');
legend('Orbit', 'Poincaré Crossings (y'' < 0)');
grid on;
hold off;

% Figure 2: Plotting only the filtered Poincaré crossing points
figure;
plot(cross_points_filtered(1:end,1), cross_points_filtered(1:end,2), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
xlabel('x');
ylabel('y');
title('Filtered Poincaré Section Crossings (y'' < 0)');
grid on;


%%
% Cargar datos
orbit1 = readmatrix('zvc_orbit1.txt');
orbit2 = readmatrix('zvc_orbit2.txt');
orbit3 = readmatrix('zvc_orbit3.txt');

orbits = {orbit1, orbit2, orbit3};

% Figura principal
figure;
hold on;
colors = lines(length(orbits)); % Generar colores distintos

% Dibujar órbitas principales
for i = 1:length(orbits)
    plot(orbits{i}(:, 1), orbits{i}(:, 2), 'DisplayName', ['Orbit ' num2str(i)], 'Color', colors(i, :));
end

% Dibujar órbita y puntos Poincaré
x0 = 2.48448275;
C = 4;

% Calcular y_prime
y_prime = dy_i(x0, xmu, C, ysign);
disp(y_prime);

% Opciones de Poincaré y cálculo
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
[cross_points, ~, ~, orbit1] = poincare(500, [x0, 0.0, 0.0, y_prime], xmu, idir, h, tol);

plot(orbit1(1, :), orbit1(2, :), 'magenta', 'LineWidth', 1.5);

xlabel('x');
ylabel('y');
title('PSP');
% xlim([-0.6, 0.8])
% ylim([-0.6, 0.6])

% % Filtrar puntos cruzados
% cross_points_filtered = cross_points(cross_points(:,4) < 0, :);
% plot(cross_points_filtered(:,1), cross_points_filtered(:,2), 'r-', 'MarkerSize', 8, 'LineWidth', 1.5); 
% 
% % Sub-ejes para el zoom con superposición
% zoom_axes = axes('Position', [0.6, 0.6, 0.3, 0.3]); % Ajusta la posición y el tamaño
% set(zoom_axes, 'Color', 'white'); % Fondo opaco blanco
% box on;
% hold on;
% 
% % Dibujar el zoom
% for i = 1:length(orbits)
%     plot(zoom_axes, orbits{i}(:, 1), orbits{i}(:, 2), 'Color', colors(i, :));
% end
% 
% plot(zoom_axes, orbit1(1, :), orbit1(2, :), 'magenta', 'LineWidth', 1.5);
% plot(zoom_axes, cross_points_filtered(:,1), cross_points_filtered(:,2), 'r-', 'MarkerSize', 8, 'LineWidth', 1.5);
% 
% % Configuración de los límites del zoom
% xlim(zoom_axes, [-8, 8]);
% ylim(zoom_axes, [-8, 8]);
% 
% hold off;

%%


% Parameters
tol = 1e-13;
s = 1e-6;
max_counter = 1000;
xmu = 0.1;
ysign = -1;
delta = 1e-3;
h = 1e-3;
idir = 1; % Forward integration
C = 4;

ranges = {
    [-1.05, -0.74],   % Range 1
    [-0.452, 0.63],   % Range 2
    [1.63, 10]        % Range 3
};

% Generate approximately 10 points for each range
x_vals = [];
for i = 1:length(ranges)
    x_vals = [x_vals, linspace(ranges{i}(1), ranges{i}(2), 20)];
end
orbits_poincare = {};
cross_points_all = [];

% Load ZVCs
orbit1 = readmatrix('zvc_orbit1.txt');
orbit2 = readmatrix('zvc_orbit2.txt');
orbit3 = readmatrix('zvc_orbit3.txt');
orbits = {orbit1, orbit2, orbit3};

% Loop through initial conditions and compute Poincaré crossings
for i = 1:length(x_vals)
    x0 = x_vals(i);

    % Compute initial y_prime
    y_prime = dy_i(x0, xmu, C, ysign);

    % Compute Poincaré Section
    options = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
    [cross_points, ~, ~, orbit] = poincare(300, [x0, 0.0, 0.0, y_prime], xmu, idir, h, tol);

    % Filter crossing points where y' < 0
    cross_points_filtered = cross_points(cross_points(:, 4) < 0, :);
    cross_points_all = [cross_points_all; cross_points_filtered(:, 1:2)]; % Collect all crossings
    orbits_poincare{i} = orbit; % Save orbit for this initial condition
end

%%
figure;
hold on;

% Plot ZVCs
for i = 1:length(orbits)
    plot(orbits{i}(:, 1), orbits{i}(:, 2), 'Color', 'blue', 'LineWidth', 1.2, 'DisplayName', ['ZVC ' num2str(i)]);
end

% Plot PSP Points
plot(cross_points_all(:, 1), cross_points_all(:, 2), 'r.', 'LineWidth', 1.5, 'DisplayName', 'Poincaré Points');
xlabel('x', 'FontSize', 14); % Make x-label larger
ylabel('y', 'FontSize', 14); % Make y-label larger
title('Zero Velocity Curves (ZVC) and Poincaré Section Points', 'FontSize', 16); % Make title larger
legend show;

% Set axis limits
xlim([-0.6, 0.8]); % Set x-axis range
ylim([-0.6, 0.6]);  % Set y-axis range

% Adjust tick sizes
set(gca, 'FontSize', 12); % Make tick labels larger

grid on;
hold off;


%%
% Parameters
tol = 1e-13; % Tolerance for Newton's method and integration
xmu = 0.1;   % µ value
C = 4;       % Jacobi constant
h = 1e-3;    % Step size for integration
ysign = -1;  % y' sign
idir = +1;   % Forward integration
x0 = 3.2171; % Initial x value

% Step 1: Compute initial y_prime
y0_prime = dy_i(x0, xmu, C, ysign);

% Step 2: Compute time to Poincaré section and full period
[cross_points, times, ~, ~] = poincare(1, [x0, 0.0, 0.0, y0_prime], xmu, idir, h, tol);
T = times(1) * 2; % Full period

% Step 3: Refine x0 and y0_prime using Newton's method
[x, y_prime] = newton_variational(x0, y0_prime, xmu, T, tol);

% Step 4: Integrate the periodic orbit with variational equations
options = odeset('RelTol', tol, 'AbsTol', tol);
[~, X] = ode45(@(t, y) f_variational_eq(t, y, xmu, 1), [0, T], ...
    [x, 0, 0, y_prime, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1], options);

% Step 5: Extract Monodromy Matrix (M)
M = reshape(X(end, 5:end), [4, 4]);

% Step 6: Compute Eigenvalues of M
eigenvalues = eig(M);

% Step 7: Display Results
fprintf('Refined Initial Conditions:\n');
fprintf('x0 = %.6f, y_prime = %.6f\n', x, y_prime);
fprintf('Full Period (T): %.6f\n', T);
fprintf('Monodromy Matrix (M):\n');
disp(M);
fprintf('Eigenvalues of M:\n');
disp(eigenvalues);
fprintf('Trace of M: %.6f\n', trace(M));
fprintf('Determinant of M: %.6f\n', det(M));

% Step 8: Plot Periodic Orbit
[t_orbit, orbit] = ode45(@(t, y) f(t, y, xmu, 1), [0, T], [x, 0, 0, y_prime], options);
figure;
plot(orbit(:, 1), orbit(:, 2), 'b-', 'LineWidth', 1.5);
xlabel('x');
ylabel('y');
title('Refined Periodic Orbit');
grid on;
%%
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
% in this case we want poincaré section to be x' = 0 and y'<0 (checked
% later)
    val = x(3);
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
    disp(['Dimensiones iniciales de orbit: ', num2str(size(orbit))]);
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
         r1 = sqrt((x0(1)-xmu)^2 + x0(2)^2);
         r2 = sqrt((x0(1)-xmu+1)^2 + x0(2)^2);
         
        derivative = 2*x0(4) + x0(1) - ((1-xmu)*(x0(1)-xmu)/(r1^3)) - xmu*(x0(1)-xmu+1)/(r2^3);
        delta = -g(approx) / derivative;
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


  





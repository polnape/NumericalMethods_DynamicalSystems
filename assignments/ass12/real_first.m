clc
clear all
close all
format long

xmu = 0.5;
ysign=-1;
delta=1e-3;
h=1e-3;
tol=1e-13;
max_iter=100;
n=200;
delta_s=1e-2;

[L1,C1,eigL1,L2,C2,eigL2,L3,C3,eigL3]=RTBP_eq_points(xmu,tol);
L4 =[ -0.5 + xmu,sqrt(3)/2]; 
C4= 2*Omega(L4,xmu,sqrt((L4(1)-xmu)^2 + L4(2)^2),sqrt((L4(1)-xmu+1)^2 + L4(2)^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C>C1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C=C1+0.1;
C = C1+0.1;
number=10000;
y=linspace(0,5,number);
i=1;
aux1=0;
aux2=0;
list_changes=[];
for y0=y
    aux2=2*Omega([xmu,y0],xmu,sqrt((xmu-xmu)^2 + y0^2),sqrt((xmu-xmu+1)^2 + y0^2))-C;
    if aux1*aux2<0
        list_changes=[list_changes, y0];
    end
    aux1=aux2;
end
aux1=0;
aux2=0;
for y0=y
    aux2=2*Omega([xmu-1,y0],xmu,sqrt((xmu-1-xmu)^2 + y0^2),sqrt((xmu-1-xmu+1)^2 + y0^2))-C;
    if aux1*aux2<0
        list_changes=[list_changes, y0];
    end
    aux1=aux2;
end
ic1=bisection_method(xmu,list_changes(1),list_changes(1)-(y(2)-y(1)),xmu,C,tol,100);
ic2=bisection_method(xmu,list_changes(2),list_changes(2)-(y(2)-y(1)),xmu,C,tol,100);
ic3=bisection_method(xmu-1,list_changes(3),list_changes(3)-(y(2)-y(1)),xmu,C,tol,100);
orbit1=pseudoarc(xmu,ic1,xmu,C,tol, max_iter,delta_s,2000);
orbit2=pseudoarc(xmu,ic2,xmu,C,tol, max_iter,delta_s,2000);
orbit3=pseudoarc(xmu-1,ic3,xmu,C,tol, max_iter,delta_s,2000);

figure
plot(orbit1(:,1),orbit1(:,2),'blue')
hold on;
plot(orbit2(:,1),orbit2(:,2),'blue')
hold on;
plot(orbit3(:,1),orbit3(:,2),'blue')
scatter([xmu,xmu-1],0,'filled','red')
xlabel('x')
ylabel('y')
title('C>C1')


writematrix(orbit1, 'zvc_orbit1.txt', 'Delimiter', 'tab');
writematrix(orbit2, 'zvc_orbit2.txt', 'Delimiter', 'tab');
writematrix(orbit3, 'zvc_orbit3.txt', 'Delimiter', 'tab');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
number=10000;
y=linspace(0,5,number);
i=1;
aux1=0;
aux2=0;
list_changes=[];
for y0=y
    aux2=2*Omega([xmu,y0],xmu,sqrt((xmu-xmu)^2 + y0^2),sqrt((xmu-xmu+1)^2 + y0^2))-C1;
    if aux1*aux2<0
        list_changes=[list_changes, y0];
    end
    aux1=aux2;
end

ic4=bisection_method(xmu,list_changes(1),list_changes(1)-(y(2)-y(1)),xmu,C1,tol,100);
ic5=bisection_method(xmu,list_changes(2),list_changes(2)-(y(2)-y(1)),xmu,C1,tol,100);
orbit4=pseudoarc(xmu,ic4,xmu,C1,tol, max_iter,delta_s,2000);
orbit5=pseudoarc(xmu,ic5,xmu,C1,tol, max_iter,delta_s,2000);

figure
plot(orbit4(:,1),orbit4(:,2),'blue')
hold on;
plot(orbit5(:,1),orbit5(:,2),'blue')
scatter([xmu,xmu-1],0,'filled','red')
xlabel('x')
ylabel('y')
title('C1')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C2<C<C1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=(C2+C1)/2;
number=10000;
y=linspace(0,5,number);
i=1;
aux1=0;
aux2=0;
list_changes=[];
for y0=y
    aux2=2*Omega([xmu,y0],xmu,sqrt((xmu-xmu)^2 + y0^2),sqrt((xmu-xmu+1)^2 + y0^2))-C;
    if aux1*aux2<0
        list_changes=[list_changes, y0];
    end
    aux1=aux2;
end

ic6=bisection_method(xmu,list_changes(1),list_changes(1)-(y(2)-y(1)),xmu,C,tol,100);
ic7=bisection_method(xmu,list_changes(2),list_changes(2)-(y(2)-y(1)),xmu,C,tol,100);
orbit6=pseudoarc(xmu,ic6,xmu,C,tol, max_iter,delta_s,2000);
orbit7=pseudoarc(xmu,ic7,xmu,C,tol, max_iter,delta_s,2000);

figure
plot(orbit6(:,1),orbit6(:,2),'blue')
hold on;
plot(orbit7(:,1),orbit7(:,2),'blue')
scatter([xmu,xmu-1],0,'filled','red')
xlabel('x')
ylabel('y')
title('C2<C<C1')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
number=10000;
y=linspace(0,5,number);
i=1;
aux1=0;
aux2=0;
list_changes=[];
for y0=y
    aux2=2*Omega([xmu,y0],xmu,sqrt((xmu-xmu)^2 + y0^2),sqrt((xmu-xmu+1)^2 + y0^2))-C2;
    if aux1*aux2<0
        list_changes=[list_changes, y0];
    end
    aux1=aux2;
end

ic8=bisection_method(xmu,list_changes(1),list_changes(1)-(y(2)-y(1)),xmu,C2,tol,100);
orbit8=pseudoarc(xmu,ic8,xmu,C2,tol, max_iter,delta_s,3000);

figure
plot(orbit8(:,1),orbit8(:,2),'blue')
hold on;
scatter([xmu,xmu-1],0,'filled','red')
xlabel('x')
ylabel('y')
title('C2')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C3<C<C2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=(C2+C3)/2;
number=10000;
y=linspace(0,5,number);
i=1;
aux1=0;
aux2=0;
list_changes=[];
for y0=y
    aux2=2*Omega([xmu,y0],xmu,sqrt((xmu-xmu)^2 + y0^2),sqrt((xmu-xmu+1)^2 + y0^2))-C;
    if aux1*aux2<0
        list_changes=[list_changes, y0];
    end
    aux1=aux2;
end

ic9=bisection_method(xmu,list_changes(1),list_changes(1)-(y(2)-y(1)),xmu,C,tol,100);
orbit9=pseudoarc(xmu,ic9,xmu,C,tol, max_iter,delta_s,5000);
figure
plot(orbit9(:,1),orbit9(:,2),'blue')
hold on;
scatter([xmu,xmu-1],0,'filled','red')
xlabel('x')
ylabel('y')
title('C3<C<C2')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=C3;
number=10000;
y=linspace(0,5,number);
i=1;
aux1=0;
aux2=0;
list_changes=[];
for y0=y
    aux2=2*Omega([xmu,y0],xmu,sqrt((xmu-xmu)^2 + y0^2),sqrt((xmu-xmu+1)^2 + y0^2))-C;
    if aux1*aux2<0
        list_changes=[list_changes, y0];
    end
    aux1=aux2;
end

ic10=bisection_method(xmu,list_changes(1),list_changes(1)-(y(2)-y(1)),xmu,C,tol,100);
orbit10=pseudoarc(xmu,ic10,xmu,C,tol, max_iter,delta_s,5000);
figure
plot(orbit10(:,1),orbit10(:,2),'blue')
hold on;
scatter([xmu,xmu-1],0,'filled','red')
xlabel('x')
ylabel('y')
title('C3')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C4<C<C3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=(C3+C4)/2;
number=10000;
y=linspace(0,5,number);
i=1;
aux1=0;
aux2=0;
list_changes=[];
for y0=y
    aux2=2*Omega([xmu,y0],xmu,sqrt((xmu-xmu)^2 + y0^2),sqrt((xmu-xmu+1)^2 + y0^2))-C;
    if aux1*aux2<0
        list_changes=[list_changes, y0];
    end
    aux1=aux2;
end

y2=linspace(-5,0,number);
i=1;
aux1=0;
aux2=0;
for y0=y2
    aux2=2*Omega([xmu,y0],xmu,sqrt((xmu-xmu)^2 + y0^2),sqrt((xmu-xmu+1)^2 + y0^2))-C;
    if aux1*aux2<0
        list_changes=[list_changes, y0];
    end
    aux1=aux2;
end
ic11=bisection_method(xmu,list_changes(1),list_changes(1)-(y(2)-y(1)),xmu,C,tol,100);
orbit11=pseudoarc(xmu,ic11,xmu,C,tol, max_iter,delta_s,2000);
ic12=bisection_method(xmu,list_changes(end),list_changes(end)-(y2(2)-y2(1)),xmu,C,tol,100);
orbit12=pseudoarc(xmu,ic12,xmu,C,tol, max_iter,delta_s,2000);
figure
hold on;
plot(orbit11(:,1),orbit11(:,2),'blue')
hold on;
plot(orbit12(:,1),orbit12(:,2),'blue')
hold on;
scatter([xmu,xmu-1],0,'filled','red')
scatter([xmu-0.5,xmu-0.5],[sqrt(3)/2,-sqrt(3)/2],'filled','green')
xlabel('x')
ylabel('y')
title('C4<C<C3')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=C4;
number=10000;
y=linspace(0,5,number);
i=1;
aux1=0;
aux2=0;
list_changes=[];
for y0=y
    aux2=2*Omega([xmu,y0],xmu,sqrt((xmu-xmu)^2 + y0^2),sqrt((xmu-xmu+1)^2 + y0^2))-C;
    if aux1*aux2<0
        list_changes=[list_changes, y0];
    end
    aux1=aux2;
end
y2=linspace(-5,0,number);
i=1;
aux1=0;
aux2=0;
for y0=y2
    aux2=2*Omega([0,y0],xmu,sqrt((0-xmu)^2 + y0^2),sqrt((0-xmu+1)^2 + y0^2))-C;
    if aux1*aux2<0
        list_changes=[list_changes, y0];
    end
    aux1=aux2;
end
figure
hold on;
scatter([xmu,xmu-1],0,'filled','red')
%scatter([xmu-0.5,xmu-0.5],[sqrt(3)/2,-sqrt(3)/2],'filled','green')
xlabel('x')
ylabel('y')
title('C4')
%%
xmu = 0.1;
C = 4;
tol = 1e-13;
max_iter = 500;
delta_t = 1e-3; % Time step for integration
num_iterates = 250; % Number of iterates for PSP


x_vals = linspace(-3, 3, 10); % Vary x along x-axis
initial_conditions = [];
for x = x_vals
    % Solve for y' using energy relation
    r1 = sqrt((x - xmu)^2);
    r2 = sqrt((x - xmu + 1)^2);
    y_prime_squared = 2 * Omega([x, 0], xmu, r1, r2) - C;
    if y_prime_squared > 0
        y_prime = -sqrt(y_prime_squared); % Select negative y'
        initial_conditions = [initial_conditions; [x, 0, 0, y_prime]];
    end
end


all_iterates = cell(size(initial_conditions, 1), 1);
for i = 1:size(initial_conditions, 1)
    ic = initial_conditions(i, :);
    iterates = compute_psp(ic, xmu, C, tol, max_iter, delta_t, num_iterates);
    all_iterates{i} = iterates;
end


figure;
hold on;
for i = 1:length(all_iterates)
    iterates = all_iterates{i};
    plot(iterates(:, 1), iterates(:, 2), ".");
end
xlabel('x');
ylabel('y');
title('Poincare Section Plot for C = 4');
grid on;
hold off;

%%

function c=bisection_method(x,y1,y2,xmu,C,tol,max_counter)
    bool=0;
    if G(x,y1,xmu,C)==0
        c=x1;
    bool=1;
    end
    if G(x,y2,xmu,C)==0 && bool==0
        c=x2;
    bool=1;
    end
    counter=0;
    while abs(y1-y2)>tol  && bool==0 && counter<max_counter
        c=(y2+y1)/2;
        if G(x,c,xmu,C)==0 && bool==0
            bool=1;
        end
        if bool==0
        if G(x,y1,xmu,C)*G(x,c,xmu,C) <0
            y2=c;
        else
            y1=c;
        end
        end
        counter=counter+1;
    end
    c=(y2+y1)/2;
end


function orbit = pseudoarc(x0,y0, xmu,C,tol, max_iter,delta_s,n)
    % Initialize variables
    r1 = sqrt((x0-xmu)^2 + y0^2);
    r2 = sqrt((x0-xmu+1)^2 + y0^2);
    x0_prime_aux=-(2*Omegay([x0,y0],xmu,r1,r2))/(2*Omegax([x0,y0],xmu,r1,r2));
    y0_prime_aux=1;
    x0_prime=x0_prime_aux/(norm([x0_prime_aux,y0_prime_aux]));
    y0_prime=y0_prime_aux/(norm([x0_prime_aux,y0_prime_aux]));
    x1=x0+delta_s*x0_prime;
    y1=y0+delta_s*y0_prime;
    orbit=zeros(n,2);
    for i=1:n %Newton_for
    
    % Initialize variables

    if i~=1
    x0=x1;
    y0=y1;
    x_prime_aux=x0_prime;
    y_prime_aux=y0_prime;
    r1 = sqrt((x0-xmu)^2 + y0^2);
    r2 = sqrt((x0-xmu+1)^2 + y0^2);
    x0_prime_aux=-(2*Omegay([x0,y0],xmu,r1,r2))/(2*Omegax([x0,y0],xmu,r1,r2));
    y0_prime_aux=1;
    x0_prime=x0_prime_aux/(norm([x0_prime_aux,y0_prime_aux]));
    y0_prime=y0_prime_aux/(norm([x0_prime_aux,y0_prime_aux]));
    x0_prime=x0_prime;
    y0_prime=y0_prime;
    if dot([x0_prime,y0_prime],[x_prime_aux,y_prime_aux])<0
        x0_prime=-x0_prime;
        y0_prime=-y0_prime;
    end
    x1=x0+delta_s*x0_prime;
    y1=y0+delta_s*y0_prime;

    end

    iter = 0; % Iteration counter
    not_finished=true;
    while not_finished && iter < max_iter 
        % Evaluate F and J at current X
        FX = F(x1,y1,x0,y0,x0_prime,y0_prime,xmu,C,delta_s);
        JX = J(x1,y1,x0,y0,x0_prime,y0_prime,xmu,C,delta_s);

        % Solve the linear system J(X) * Delta = -F(X) for Delta
        Delta = -JX \ FX; % MATLAB's backslash operator for solving linear systems
        % Update the solution
        x1 = x1 + Delta(1);
        y1= y1+Delta(2);
        
        % Check for convergence
        if norm(FX, 'fro') < tol
            not_finished=false;
        end

        % Increment iteration counter
        iter = iter + 1;
    end
    % If max iterations are reached without convergence, warn the user
    orbit(i,:)=[x1,y1];
    end
end

function res=G(x,y,xmu,C)
    r1 = sqrt((x-xmu)^2 + y^2);
    r2 = sqrt((x-xmu+1)^2 + y^2);
    res=2*Omega([x y],xmu,r1,r2)-C;
end
function matrix = F(x1,y1,x0,y0,x0_prime,y0_prime,xmu,C,delta_s)
    matrix=[G(x1,y1,xmu,C);((x1-x0)*x0_prime+(y1-y0)*y0_prime)-delta_s];
end
function matrix=J(x1,y1,x0,y0,x0_prime,y0_prime,xmu,C,delta_s)
r1 = sqrt((x1-xmu)^2 + y1^2);
r2 = sqrt((x1-xmu+1)^2 + y1^2);
matrix=[2*Omegax([x1,y1],xmu,r1,r2),2*Omegay([x1,y1],xmu,r1,r2);x0_prime,y0_prime];
end




function [L1,C1,eigL1,L2,C2,eigL2,L3,C3,eigL3]=RTBP_eq_points(xmu,tol)
    %%%%%L1%%%%%
    xi=(xmu/(3*(1-xmu)))^(1/3);
    aux=0;
    while(abs(xi-aux)>tol)
        aux=xi;
        xi=F1(xi,xmu);
    end
    L1=[xmu-1+xi,0,0,0];    

    %%%%%L2%%%%%
    xi=(xmu/(3*(1-xmu)))^(1/3);
    aux=0;
    while(abs(xi-aux)>tol)
        aux=xi;
        xi=F2(xi,xmu);
    end
    L2=[xmu-1-xi,0,0,0];

    %%%%%L3%%%%%
    xi=1-(7*xmu)/(12);
    aux=0;
    while(abs(xi-aux)>tol)
        aux=xi;
        xi=F3(xi,xmu);
    end
    L3=[xmu+xi,0,0,0];
    xL1=L1(1);
    xL2=L2(1);
    xL3=L3(1);

    C1 = 2*Omega([xL1 0],xmu,abs(xL1-xmu),abs(xL1-xmu+1));
    C2 = 2*Omega([xL2 0],xmu,abs(xL2-xmu),abs(xL2-xmu+1));
    C3 = 2*Omega([xL3 0],xmu,abs(xL3-xmu),abs(xL3-xmu+1));

    DG1 = zeros(4,4); DG1(1,3)=1; DG1(2,4)=1; DG1(3,4)=2; DG1(4,3)=-2;
    DG2 = DG1; DG3 = DG1;
    DG1(3,1) = Omegaxx([xL1 0],xmu);
    DG1(3,2) = Omegaxy([xL1 0],xmu);
    DG1(4,1) = Omegaxy([xL1 0],xmu);
    DG1(4,2) = Omegayy([xL1 0],xmu);
    eigL1 = eig(DG1);
    DG2(3,1) = Omegaxx([xL2 0],xmu);
    DG2(3,2) = Omegaxy([xL2 0],xmu);
    DG2(4,1) = Omegaxy([xL2 0],xmu);
    DG2(4,2) = Omegayy([xL2 0],xmu);
    eigL2 = eig(DG2);
    DG3(3,1) = Omegaxx([xL3 0],xmu);
    DG3(3,2) = Omegaxy([xL3 0],xmu);
    DG3(4,1) = Omegaxy([xL3 0],xmu);
    DG3(4,2) = Omegayy([xL3 0],xmu);
    eigL3 = eig(DG3);
end

function res=Omegax(x,mu,r1,r2)
res=x(1) - ((1-mu)*(x(1)-mu)/(r1^3)) - mu*(x(1)-mu+1)/(r2^3);
end
function res=Omegay(x,mu,r1,r2)
res=+ x(2)*(1 - (1-mu)/(r1^3) - mu/(r2^3));
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
  res =  (mu-1)/((mu-x(1))^2+x(2)^2)^(3/2)-mu/((x(1)-mu+1)^2+x(2)^2)^(3/2)-(3*(2*mu-2*x(1))^2*(mu-1))/(4*((mu-x(1))^2+x(2)^2)^(5/2))+(3*mu*(2*x(1)-2*mu+2)^2)/(4*((x(1)-mu+1)^2+x(2)^2)^(5/2))+1;
end

function res = Omegayy(x,mu)
  res = (mu-1)/((mu-x(1))^2+x(2)^2)^(3/2)-mu/((x(1)-mu+1)^2+x(2)^2)^(3/2)+(3*mu*x(2)^2)/((x(1)-mu+1)^2+x(2)^2)^(5/2)-(3*x(2)^2*(mu-1))/((mu-x(1))^2+x(2)^2)^(5/2)+1;
end

function res = Omegaxy(x,mu)
  res = (3*x(2)*(2*mu-2*x(1))*(mu-1))/(2*((mu-x(1))^2+x(2)^2)^(5/2))+(3*mu*x(2)*(2*x(1)-2*mu+2))/(2*((x(1)-mu+1)^2+x(2)^2)^(5/2));
end


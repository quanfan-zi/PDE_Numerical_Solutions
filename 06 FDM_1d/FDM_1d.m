function FDM_1d_DN
%%
% y*=sin(pi/2*x)   in [0,1],
% -y''+pi^2/4*y=pi^2/2*sin(pi/2*x).
% y(0)=0, y'(1)=0.
tic
format long;
clear;clc;
close all;
pde.start_point=0;
pde.end_point=1;
% 边值条件的处理方法
% 1表示向后差商, 2表示中心差商, 3表示有限体积法
pde.boundary_method=2;

% 计算收敛阶
 [Max_Order, L2_Order,H1_Order]=convergenceOrder(pde)

% 画出函数图像 
% N=10;
% [Max_error,L2_error,H1_error]=run_main(N,pde)
toc
end

function [Max_error,L2_error,H1_error]=run_main(N,pde)
%%
boundary_method=pde.boundary_method; 
Length=pde.end_point-pde.start_point;
h=Length/N;
x=(pde.start_point:h:pde.end_point)';  

% 生成刚度矩阵和右端项
[A,F]=stiffnessMatrix(N,h,x);

% 处理边界条件
[A,F]=boundaryMatrix(A,F,x,N,h,boundary_method);

% 解方程组
u=A\F;

% 计算误差
[Max_error,L2_error,H1_error]=errorEstimate(u,x,N,h);

% 函数图像
plotFigure(x,u,pde)
end

function [A,F]=stiffnessMatrix(N,h,x)
%% Generate stiffness matrix and right hands
A=zeros(N+1,N+1);
F=zeros(N+1,1);

p=1; q=pi^2/4;
for i=2:N
    A(i,i-1)=-p/h;
    A(i,i)=2*p/h+h*q;
    A(i,i+1)=-p/h;
    F(i)=h*f(x(i)); 
end
end

function [A,F]=boundaryMatrix(A,F,x,N,h,boundary_method)
%% boundary condition
p=1; q=pi^2/4;
%Neumann边值条件
du_value=0;
if boundary_method==1          %向后差分
    A(N+1,N)=-1;
    A(N+1,N+1)=1;
    F(N+1)=h*du_value;              
elseif boundary_method==2   %中心差分    
    A(N+1,N)=-p/h;
    A(N+1,N+1)=p/h+q*h/2;
    F(N+1)=h*f(x(N+1))/2;
else                                            %有限体积法
    A(N+1,N)=-p/h;
    A(N+1,N+1)=p/h+q*h/2;
    F(N+1)=h*f(x(N+1)-h/4)/2;
end  

%Dirichelt边值条件
u_value=0;
A(1,1)=1;
F(1)=u_value;
end

function [Max_error,L2_error,H1_error]=errorEstimate(u,x,N,h)
%% Error estimate
u_true=true_sol(x);
u_err=u-u_true;

Max_error=max(abs(u_err));

L2_error=sum(u_err.^2)*h;
L2_error= L2_error^(1/2);

d_u_err=(u_err(2:N+1)-u_err(1:N))/h;
H1_error=sqrt(sum(d_u_err.^2)/h);
%H1_err_half=sum(d_u_err.^2)*h;
%H1_error=sqrt(H1_err_half+L2_error^2);
end

function plotFigure(x,u,pde)
    figure(1)
    u_true=true_sol(x);
    plot(x,u_true,'k-','linewidth',2)
    hold on;
    scatter(x,u,'filled','SizeData', 40);
    title('Numerical Solution (N=10)')
    hold off;
end

function y=true_sol(x)
%% Exact solution
y=sin(pi*x/2);
end

function y=f(x)
%% Right hand function 
y=pi^2/2*sin(pi*x/2);
end

function [Max_Order, L2_Order,H1_Order]=convergenceOrder(pde)
n=20;
L2_error_array=zeros(n,1);
H1_error_array=zeros(n,1);
Max_error_array=zeros(n,1);
L2_Order=zeros(n-1,1);
H1_Order=zeros(n-1,1);
Max_Order=zeros(n-1,1);
N_vector=10*(1:n);

for i=1:n
    [Max, L2, H1]=run_main(N_vector(i),pde);
    L2_error_array(i)=L2;
    H1_error_array(i)=H1;
    Max_error_array(i)=Max;
end
for i=1:n-1
    L2_Order(i)=-log2(L2_error_array(i)/L2_error_array(i+1))/log2(i/(i+1));
    H1_Order(i)=-log2(H1_error_array(i)/H1_error_array(i+1))/log2(i/(i+1));
    Max_Order(i)=-log2(Max_error_array(i)/Max_error_array(i+1))/log2(i/(i+1));
end


end

%% 测试方程
% u*=[exp(-pi^2*t)+t]*sin(pi*x)   in [0,1]*[0,1],
% \partial u/ \partial t = \partial^2 u/ \partial x^2+sin(pi*x)+pi^2*t*sin(pi*x)
% u(x,0)=sin(pi*x).
% u(0,t)=u(pi,t)=0

function FDM_parabolic
%% 参数函数, 画函数图像或求收敛阶
tic
format short;    clear;  clc;  close all;
pde.start_point=0;  pde.end_point=1;            % 横轴区间65741121
pde.start_time=0;  pde.end_time=1;               % 时间区间
pde.subdivision=640;                                         % 横轴剖分数
pde.r=1/2;                                                         % 网比r=at/h^2
% pde.method=1表示向前差分格式, 2表示向后差分格式, 3表示六点对称格式

% 单独画图像
% pde.method=1;
% [sDof, L2_error]=run_main(pde)

% 求收敛阶
convergenceOrder(pde)

toc
end

function [sDof, L2_error]=run_main(pde)
%% 求解函数值的主程序
method=pde.method;
N=pde.subdivision;
h=(pde.end_point-pde.start_point)/N;
t=pde.r*h^2;
M=round((pde.end_time-pde.start_time)/t);
r=pde.r;
x_array=pde.start_point:h:pde.end_point;
t_array=pde.start_time:t:pde.end_time;
sDof=(M+1)*(N+1);

% 生成解矩阵
u=SolutionInitial(M,N,x_array,t_array);

% 生成迭代矩阵
[AA, C]=MatrixGenerate(r,N,method);

% 求解函数值
for i=2:M+1
    F=RightFuction(x_array(2:N)',i*t*ones(N-1,1));
    tempu=C*u(i-1,2:N)'+t*AA*F;
    u(i,2:N)=tempu';
end

% 画图
%PlotFigure(x_array,t_array,u)

% 计算误差
L2_error=errorEstimate(u,h,x_array);
end

function u=SolutionInitial(M,N,x_array,t_array)
%% 初始化解矩阵
    u=zeros(M+1,N+1);
    u(1,:)=InitialFunction(x_array);                      % 初值
    u(:,1)=BoundaryFunction(t_array)';              % 边值
    u(:,N+1)=BoundaryFunction(t_array)';  
end

function [AA, C]=MatrixGenerate(r,N,method)
%% 生成迭代矩阵
%公式: U^(k+1)=C*U^k+tA^(-1)*F
S=toeplitz([0, 1, zeros(1, N-3)]);
if method == 1                                                 % 向前差分法
    C=(1-2*r)*eye(N-1)+r*S;
    AA=eye(N-1);
elseif method == 2                                           % 向后差分法
    C=((1+2*r)*eye(N-1)-r*S)^(-1);
    AA=C;
elseif method == 3                                             % 六点对称法
    AA=((1+r)*eye(N-1)-r/2*S)^(-1);
    C=AA*((1-r)*eye(N-1)+r/2*S);
end
end

function PlotFigure(x_array,t_array,u)
%% 作图
[X, T]=meshgrid(x_array,t_array);
subplot(1,2,1);
surf(X,T,u)
colormap(jet)
alpha(0.5)
view(-30,20)
colorbar
title('numerical solution');

subplot(1,2,2);
u_true=TrueSolution(X,T);
surf(X,T,u_true)
colormap(jet)
alpha(0.5)
view(-30,20)
colorbar
title('true solution');
end

function L2_error=errorEstimate(u,h,x_array)
%% 求t=1时的L2范数误差
u_num=u(end,:);
u_true=TrueSolution(x_array,ones(size(x_array)));
L2_error=sum((u_num-u_true).^2*h);
L2_error=sqrt(L2_error);
end

function u=InitialFunction(x)
    u=sin(pi*x);
end

function u=BoundaryFunction(t)
    u=0*t;
end

function u=RightFuction(x,t)
    u=sin(pi*x)+pi^2*t.*sin(pi*x);
end

function u =TrueSolution(X,T)
    u=(exp(-pi^2*T)+T).*sin(pi*X);
end

function convergenceOrder(pde)
N=[5,10,20,40,80,160,240,320];
sDof_array=zeros(1,8);
L2_error_array1=zeros(1,8);
L2_error_array2=zeros(1,8);
L2_error_array3=zeros(1,8);
L2_Order1=zeros(1,7);
L2_Order2=zeros(1,7);
L2_Order3=zeros(1,7);
for i=1:8
    pde.subdivision=N(i);
    
    pde.method=1;
    [sDof, L2_error]=run_main(pde);
    L2_error_array1(i)=L2_error;
    
    pde.method=2;
    [~, L2_error]=run_main(pde);
    L2_error_array2(i)=L2_error;
    
    pde.method=3;
    [~, L2_error]=run_main(pde);
    L2_error_array3(i)=L2_error;    
    sDof_array(i)=sDof;
end
for i=1:8-1
    fenmu=log2(N(i)/N(i+1));
    L2_Order1(i)=-log2(L2_error_array1(i)/L2_error_array1(i+1))/fenmu;
    L2_Order2(i)=-log2(L2_error_array2(i)/L2_error_array2(i+1))/fenmu;
    L2_Order3(i)=-log2(L2_error_array3(i)/L2_error_array3(i+1))/fenmu;
end
L2_Order1
L2_Order2
L2_Order3

figure(2)
loglog(N, L2_error_array1,'m-*','linewidth',1.5)
hold on
loglog(N, L2_error_array2,'-x','linewidth',1.5)
loglog(N, L2_error_array3,'b-d','linewidth',1.5)
loglog(N,0.5*N.^(-2),'k-d','linewidth',2)
hold off
legend('向前差分格式','向后差分格式','六点对称格式','2阶收敛速度','Location','NorthEast');
title('Error Estimation(loglog)','fontsize',14)
end






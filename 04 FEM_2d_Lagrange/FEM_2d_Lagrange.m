function FEM_2d_1p_sq_ND
%%
% u*=xy+sin(pi*x)sin(pi*y)   in [0,1]^2,
% -\Laplace u-2*pi^2*u=-2pi^2*xy
% u(0,y)=u(x,0)=0
% u_x(1,y)=y-pi*sin(pi*y)
% u_y(x,1)=x-pi*sin(pi*x)
tic
format long;
clear;clc;
close all;
pde.start_point_x=0;
pde.end_point_x=1;
pde.start_point_y=0;
pde.end_point_y=1;


% 计算收敛阶
 %n_iter=6;
 %[L2_convergenceOrder,H1_convergenceOrder]=convergenceOrder(n_iter,pde)

% 画出函数图像 
N=20
 [sDof, L2_error,H1_error, Mat_con]=run_main(N,pde)
toc
end

function [sDof, L2_error,H1_error, Mat_con]=run_main(N,pde)
%% intput N pde
%% output sDof,L2_error,H1_error, Mat_con
% 分别表示自由度, L2误差, H1误差和刚度矩阵条件数

% 处理边值条件的方法
boundary_method=1;  % 1  2  3
domainLenght=pde.end_point_x-pde.start_point_x;
h=domainLenght/N;

% 区域剖分
[p,e,t]=meshGeneration(pde,N);
sDof=size(p,2);
% sDof表示节点数量, 也就是自由度
% p是点, 每一列表示划分后点的x, y坐标
% t是三角形, 前三行表示三角形顶点在p中的指标(列数)
% e是边, 只包含边界上的边(因为处理边值条件与内部边无关), 前两行表示起点和终点在p中的指标 
%第三行表示小边属于那个大边, 其中第一条边x=0, 第二条边x=1, 第三条边y=0, 第四条边y=1

% 生成刚度矩阵
A=stiffnessMatrix(p,e,t,h);

% 生成右端项
F=rightHands(p,e,t,h);

% 处理边界条件
[A,F]=boundaryMatrix(A,F,p,e,boundary_method);
Mat_con=cond(A,2);

% solve A\F
u=solveAF(A,F,N,pde,boundary_method);

% 计算误差
[L2_error,H1_error]=errorEstimate(u,p,t,h);

% 画图
 plotFigure(N,pde,u)
end

function [p,e,t]=meshGeneration(pde,N)
%% Generate mesh
x_start=pde.start_point_x;
x_end=pde.end_point_x;
y_start=pde.start_point_y;
y_end=pde.end_point_y;

h_x=(x_end-x_start)/N;
h_y=(y_end-y_start)/N;

[X,Y] = meshgrid(x_start:h_x:x_end,y_start:h_y:y_end);
length=(N+1)*(N+1);
x=reshape(X,length,1);
y=reshape(Y,length,1);
p=[x,y]';

k=1:(N+1)*(N+1);
k=reshape(k,N+1,N+1);
k=k(1:N,1:N);
k=reshape(k,1,N*N);
t=[k;k+1;k+2+N;k+N+1];

loc_1=find(p(1,:)==x_start);
loc_2=find(p(1,:)==x_end);
loc_3=find(p(2,:)==y_start);
loc_4=find(p(2,:)==y_end);
e=[loc_1(1:N),loc_2(1:N),loc_3(1:N),loc_4(1:N);
    loc_1(2:N+1),loc_2(2:N+1),loc_3(2:N+1),loc_4(2:N+1);
    ones(1,N), 2*ones(1,N), 3*ones(1,N), 4*ones(1,N)];
% 第一条边x=0, 第二条边x=1, 第三条边y=0, 第四条边y=1

end

function A=stiffnessMatrix(p,~,t,h)
%% Generate stiffness matrix
sDof=size(p,2);     %节点数量
 A=sparse(sDof,sDof);
% A=zeros(sDof,sDof);
nel=size(t,2);         %小区间数量

% 单元刚度矩阵用mathematica提前算出来
a1=2/3-2*pi^2*h^2/9;
a2=-1/6-pi^2*h^2/9;
a3=-1/3-pi^2*h^2/18;
K=[a1,a2,a3,a2;
    a2,a1,a2,a3;
    a3,a2,a1,a2;
    a2,a3,a2,a1];
for iel=1:nel
    nd=t(1:4,iel);
    A(nd,nd)=A(nd,nd)+K;
end
end

function F=rightHands(p,~,t,h)
%% Generate right hands
sDof=size(p,2);
F=zeros(sDof,1);   %F的维数等于变量自由度
nel=size(t,2);         %剖分后矩形个数

[x,y]=get_integration_points([0,0], h);
b_value_00=basis00(x,y,h);
b_value_01=basis01(x,y,h);
b_value_10=basis10(x,y,h);
b_value_11=basis11(x,y,h);
for iel=1:nel
    % 找到小区间的四个顶点的指标
    nd=t(1:4,iel);
    a=p(:,nd(1));  %小区间左下角顶点坐标
    [x,y]=get_integration_points(a, h);
    f_value=rightFunction(x,y);
    f=zeros(4,1);
    f(1)=square_quadrature(f_value.*b_value_00,h);
    f(2)=square_quadrature(f_value.*b_value_01,h);
    f(3)=square_quadrature(f_value.*b_value_11,h);
    f(4)=square_quadrature(f_value.*b_value_10,h);
    F(nd)=F(nd)+f;
end
end

function [A,F]=boundaryMatrix(A,F,p,e,boundary_method)
%% 处理边值条件
values_1=zeros(5,1);
values_2=zeros(5,1);

% 处理Neumann边值条件
% 第四条边, y=1, 只和x有关
loc_neuEdges_4=find(e(3,:)==4);
neuEdges_4=e(:,loc_neuEdges_4);
nneu_4=size(loc_neuEdges_4,2);
for iedge=1:nneu_4
    x1=p(1,neuEdges_4(1,iedge));
    x2=p(1,neuEdges_4(2,iedge));
    min_x=min(x1,x2);
    max_x=max(x1,x2);
    x=get_integration_points_1d(min_x,max_x);
    for k=1:5
        % 赋值
        values_1(k)=nboundaryFunction_4(x(k), 1) * (x(k)-x2)/(x1-x2);
    end
% 更新右端向量
F(neuEdges_4(1,iedge))=F(neuEdges_4(1,iedge))+gaussint5(min_x,max_x,values_1);
    for k=1:5
          values_2(k)=nboundaryFunction_4(x(k), 1) * (x(k)-x1)/(x2-x1);
    end
F(neuEdges_4(2,iedge))=F(neuEdges_4(2,iedge))+gaussint5(min_x,max_x,values_2);
end

% 第二条边, x=1,只和y有关
loc_neuEdges_2=find(e(3,:)==2);
neuEdges_2=e(:,loc_neuEdges_2);
nneu_2=size(loc_neuEdges_2,2);
for iedge=1:nneu_2
    y1=p(2,neuEdges_2(1,iedge));
    y2=p(2,neuEdges_2(2,iedge));
    min_y=min(y1,y2);
    max_y=max(y1,y2);
    y=get_integration_points_1d(min_y,max_y);
    for k=1:5
        % 赋值
        values_1(k)=nboundaryFunction_2(1,y(k)) * (y(k)-y2)/(y1-y2);
    end
% 更新右端向量
F(neuEdges_2(1,iedge))=F(neuEdges_2(1,iedge))+gaussint5(min_y,max_y,values_1);
    for k=1:5
          values_2(k)=nboundaryFunction_2(1,y(k)) * (y(k)-y1)/(y2-y1);
    end
F(neuEdges_2(2,iedge))=F(neuEdges_2(2,iedge))+gaussint5(min_y,max_y,values_2);
end

% 处理第一边值条件
% 找到x=0或者y=0的节点
loc_d=find(p(1,:)==0|p(2,:)==0);
if boundary_method==1 
    x=p(1,loc_d);
    y=p(2,loc_d);
    A(loc_d,:)=0;
    A(loc_d,loc_d)=eye(size(loc_d,2));    %3*9变为3*3
    F(loc_d)=boundaryFunction(x,y)';
end
F;
end

function u=solveAF(A,F,N,pde,boundary_method)
%% solve A\F
if boundary_method==1||boundary_method==3
    u=A\F;
else
    u=zeros(N+1,1);
    u(1:N)=A\F;
    u(N+1)=pde.right_boundary_Dirichlet;
end
end

function [L2_error,H1_error]=errorEstimate(u,p,t,h)
%% Error estimate
% L2=sqrt(/iint(u-u_h)^2dxdy)

L2_error=0;
H1_half=0;

[x,y]=get_integration_points([0,0], h);
b_value_00=basis00(x,y,h);
b_value_01=basis01(x,y,h);
b_value_10=basis10(x,y,h);
b_value_11=basis11(x,y,h);
% 如果想提升程序运行速度， 应该将导数节点处的数值也写在循环外面
% 但是这样太冗长太丑.
for i=1:size(t,2)
    ssquare=t(1:4,i);
    a=p(:,ssquare(1));
    [x,y]=get_integration_points(a,h);
    % L2误差
    values_numerial=b_value_00.*u(ssquare(1))+b_value_01.*u(ssquare(2))+...
        b_value_11.*u(ssquare(3))+b_value_10.*u(ssquare(4));
    values_true=true_sol(x,y);
    values=(values_numerial-values_true).^2;
    L2_error=L2_error+square_quadrature(values,h);
    % H1误差
    values_numerial_gradient=basis00_gradient(x,y,h).*u(ssquare(1))+...
        basis01_gradient(x,y,h).*u(ssquare(2))+...
        basis11_gradient(x,y,h).*u(ssquare(3))+basis10_gradient(x,y,h).*u(ssquare(4));
    exactValues_gradient=ture_sol_gradient(x,y);
    values_gradient=(values_numerial_gradient-exactValues_gradient).^2;
    H1_half=H1_half+square_quadrature(values_gradient(1,:),h)+...
                                    square_quadrature(values_gradient(2,:),h);
end
L2_error=sqrt(L2_error);
H1_error=sqrt(L2_error^2+H1_half);
end

function plotFigure(N,pde,u)
%% Plot figure
h=(pde.end_point_x-pde.start_point_x)/N;
x=pde.start_point_x:h:pde.end_point_x;
y=pde.start_point_y:h:pde.end_point_y;
[X, Y]=meshgrid(x,y);

% subplot(1,3,1);
% U=reshape(u,N+1,N+1);
% surf(X,Y,U)
% colormap(jet)
% alpha(0.5)
% view(-15,20)
% title('numerical solution(N=5)');

subplot(1,2,1);
U=reshape(u,N+1,N+1);
surf(X,Y,U)
colormap(jet)
alpha(0.5)
view(-15,20)
colorbar
title('numerical solution');

subplot(1,2,2);
%pdemesh(p,e,t,z);
Z=X.*Y+sin(pi*X).*sin(pi*Y);
surf(X,Y,Z)
colormap(jet)
alpha(0.5)
view(-15,20)
% grid off
colorbar
title('exact solution');
end

function z=true_sol(x,y)
%% Exact solution
z=x.*y+sin(pi*x).*sin(pi*y);
end

function z=ture_sol_gradient(x,y)
%% Derivatives of exaxt solution
z_x=y+pi*cos(pi*x).*sin(pi*y);
z_y=x+pi*cos(pi*y).*sin(pi*x);
z=[z_x;z_y];
end

function z=rightFunction(x,y)
%% 偏微分方程右端函数
% z=2*pi*pi*sin(pi.*x).*sin(pi.*y);
z=-2*pi*pi*x.*y;
end

function z=boundaryFunction(~,~)
%% Boundary function
% z=sin(pi.*x).*sin(pi.*y);
z=0;
end

function z=nboundaryFunction_4(x,~)
%% Boundary function
z=x-pi*sin(pi*x);
end

function z=nboundaryFunction_2(~,y)
%% Boundary function
% z=sin(pi.*x).*sin(pi.*y);
z=y-pi*sin(pi*y);
end

function z=basis00(x, y, h)
%% 区间[0,h]^2上点(0,0)的基函数
z=(x/h-1).*(y/h-1);
end

function z=basis01(x, y, h)
%% 区间[0,h]^2上点(0,h)的基函数
z=(1-x/h).*(y/h);
end

function z=basis10(x, y, h)
%% 区间[0,h]^2上点(h,0)的基函数
z=(x/h).*(1-y/h);
end

function z=basis11(x, y, h)
%% 区间[0,h]^2上点(h,h)的基函数
z=(x/h).*(y/h);
end

function z=basis00_gradient(x, y, h)
%% 区间[0,h]^2上点(0,0)的基函数
z_x=-1/h*(1-y/h);
z_y=-1/h*(1-x/h);
z=[z_x; z_y];
end

function z=basis01_gradient(x, y, h)
%% 区间[0,h]^2上点(0,0)的基函数
z_x=-y/h^2;
z_y=1/h*(1-x/h);
z=[z_x; z_y];
end

function z=basis10_gradient(x, y, h)
%% 区间[0,h]^2上点(0,0)的基函数
z_x=1/h*(1-y/h);
z_y=-1/h^2*x;
z=[z_x; z_y];
end

function z=basis11_gradient(x, y, h)
%% 区间[0,h]^2上点(0,0)的基函数
z_x=1/h^2*y;
z_y=1/h^2*x;
z=[z_x; z_y];
end

function [x,y]=get_integration_points(point,h)
%% 把[-1,1]*[-1,1]区间上的求积结点线性变换到指定区间上
r=sqrt(15)/5;
x=[-r,0,r,-r,0,r,-r,0,r];  %values为行向量
y=[-r,-r,-r,0,0,0,r,r,r];
x=h/2*(x+1)+point(1);
y=h/2*(y+1)+point(2);
end

function y = square_quadrature(values,h)
%% Quadrature formula
weights = reshape([5/9,8/9,5/9]'*[5/9,8/9,5/9],9,1); %weights列向量
y=h^2/4*(values*weights);
end

function x=get_integration_points_1d(a,b)
%% Generate intergration points(1D)
location=[0.046910077030668 0.230765344947158 0.5  0.769234655052842 0.953089922969332];
for i=1:5
x(i)=a+location(i)*(b-a);
end
end

function y=gaussint5(a,b,values)
weights=[0.236926885056189 0.478628670499366 0.568888888888889 0.478628670499366 0.236926885056189];
y=(weights*values)*((b-a)/2);
end

function [L2_convergenceOrder,H1_convergenceOrder]=convergenceOrder(N,pde)
%% Get integration points
L2_error_array=zeros(N,1);
H1_error_array=zeros(N,1);
Mat_con_array=zeros(N,1);
num_freedom_array=zeros(N,1);
L2_convergenceOrder=zeros(N-1,1);
H1_convergenceOrder=zeros(N-1,1);
Mat_convergenceOrder=zeros(N-1,1);
N_vector=5*2.^(0:N-1);

for i=1:N
    [sDof, L2_error, H1_error, Mat_con]=run_main(N_vector(i),pde);
    num_freedom_array(i)=sDof;
    L2_error_array(i)=L2_error;
    H1_error_array(i)=H1_error;
    Mat_con_array(i)=Mat_con;
end
for i=1:N-1
    L2_convergenceOrder(i)=log2(L2_error_array(i)/L2_error_array(i+1));
    H1_convergenceOrder(i)=-log2(H1_error_array(i)/H1_error_array(i+1));
    Mat_convergenceOrder(i)=-log2(Mat_con_array(i)/Mat_con_array(i+1));
end
figure(2)
%
subplot(1,3,1)
%plot(num_freedom_array,H1_error_array,'-o',num_freedom_array,L2_error_array,'-*')
plot(num_freedom_array,L2_error_array,'-*')
% hleg1=legend('H1-norm','L2-norm','Location','NorthEast');
title('Error Estimation')
subplot(1,3,2)
loglog(num_freedom_array,L2_error_array,'m-*')
%loglog(num_freedom_array,H1_error_array,'-o',num_freedom_array,L2_error_array,'-*')
% hleg1=legend('H1-norm','L2-norm','Location','NorthEast');
title('Error Estimation(loglog)')
subplot(1,3,3)
plot(num_freedom_array(1:N-1),L2_convergenceOrder,'r-*')
legend('L2-norm','Location','NorthEast');
%plot(num_freedom_array(1:N-1),H1_convergenceOrder,'-o',num_freedom_array(1:N-1),L2_convergenceOrder,'-*')
%hleg1=legend('H1-norm','L2-norm','Location','NorthEast');
title('Convergence Order')
%plot(num_freedom,H1_error_array,'-o',num_freedom,exp(-1*log(num_freedom)+1.5),'',num_freedom,L2_error_array,'-*',num_freedom,exp(-2*log(num_freedom)+1),'-.')
%xlabel('degree of freedom in log scale');
%ylabel('error in log scale');



figure(3) %
subplot(1,2,1)
loglog(num_freedom_array, Mat_con_array, '-*')
title('Matrix Condition','FontSize', 14)
subplot(1,2,2)
plot(num_freedom_array(1:N-1), Mat_convergenceOrder, 'm-*')
title('Divergence Order','FontSize', 14)

%axis([1e2,1e4,1e-4,1e1])
end
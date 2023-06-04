function FEM_1d_2p_for_ND
%%
%-y"+pi^2/4*y=(pi^2)/2*sin(pi*x/2)   in [0,1],
%y=sin(pi*x/2)
%y(0)=0, y'(0)=0;
%
tic
%format long;
pde.left_boundary_Dirichlet=0;
pde.right_boundary_Neumann=0;
pde.start_point=0;
pde.end_point=1;

 % 计算收敛阶
n_iter=20;
[L2_convergenceOrder,H1_convergenceOrder]=convergenceOrder(n_iter,pde)

% 画出函数图像 
 %N=10;
 %[L2_error,H1_error, Mat_con]=run_main(N,pde)
toc
end

function [L2_error,H1_error, Mat_con]=run_main(N,pde)
%%
% 处理边值条件的方法
boundary_method=1;  % 1  2  3
Length=pde.end_point-pde.start_point;
h=Length/N;
x=(pde.start_point:h:pde.end_point)';  %x是区间剖分, 不考虑半节点

% [-1,1]区间上五点Guass积分公式的取点和权重
positions=[-0.906179845938664,-0.538469310105683,0,0.538469310105683,0.906179845938664];
weights=[0.236926885056189;0.478628670499366;0.568888888888889;0.478628670499366;0.236926885056189];

% 生成刚度矩阵
A=stiffnessMatrix(N,h,positions,weights);
Mat_con=cond(A,2);

% 生成右端项
F=rightHands(x,N,h,positions,weights);

% 处理边界条件
[A,F]=boundaryMatrix(A,F,x,N,h,pde,boundary_method);

% solve A\F
u=solveAF(A,F,N,pde,boundary_method);

% 计算误差
[L2_error,H1_error]=errorEstimate(u,x,N,h,positions,weights);

% plot figure
xx=(pde.start_point:h/2:pde.end_point)';
plotFigure(xx,u,pde)
end

function A=stiffnessMatrix(N,h,positions,weights)
%% Generate stiffness matrix
A=zeros(2*N+1,2*N+1);

% 边值问题中的系数
p1=1; p2=pi^2/4;

% 用数值积分计算每个小区间上基函数的积分
% 找到对应函数值
inter_point=get_integration_points(0, h, positions);
f_l_values=basis_l(inter_point, 0, h);
f_r_values=basis_r(inter_point, 0, h);
f_m_values=basis_m(inter_point, 0, h);
f_l_der_values=basis_l_der(inter_point, 0, h);
f_r_der_values=basis_r_der(inter_point, 0, h);
f_m_der_values=basis_m_der(inter_point, 0, h);
% 计算高斯积分值
temp=p1*f_l_der_values.^2+p2*f_l_values.^2;
a1=Gauss_Quadrature_formula(h, temp, weights);
temp=p1*f_l_der_values.*f_r_der_values+p2*f_l_values.*f_r_values;
a2=Gauss_Quadrature_formula(h, temp, weights);
temp=p1*f_m_der_values.^2+p2*f_m_values.^2;
a3=Gauss_Quadrature_formula(h, temp, weights);
temp=p1*f_l_der_values.*f_m_der_values+p2*f_l_values.*f_m_values;
a4=Gauss_Quadrature_formula(h, temp, weights);

% K表示单元刚度矩阵
K=[a1,a4,a2;
      a4,a3,a4;
      a2,a4,a1];
for i=1:N
    A(2*i-1:2*i+1,2*i-1:2*i+1)=A(2*i-1:2*i+1,2*i-1:2*i+1)+K;
end
end

function F=rightHands(x,N,h,positions,weights)
%% Generate right hands
F=zeros(2*N+1,1);
inter_point=get_integration_points(0, h, positions);
b_values_l=basis_l(inter_point, 0, h);
b_values_r=basis_r(inter_point, 0, h);
b_values_m=basis_m(inter_point, 0, h);
for i=1:N
    inter_point=get_integration_points(x(i),x(i+1),positions);
    f_values=f(inter_point);
    F(2*i-1)=F(2*i-1)+Gauss_Quadrature_formula(h,f_values.*b_values_r,weights);
    F(2*i)=F(2*i)+Gauss_Quadrature_formula(h,f_values.*b_values_m,weights); 
    F(2*i+1)=F(2*i+1)+Gauss_Quadrature_formula(h,f_values.*b_values_l,weights);
end
end

function [A,F]=boundaryMatrix(A,F,x,N,h,pde,boundary_method)
%% boundary condition
if boundary_method==1
    A(1,:)=0;
    A(1,1)=1;
    F(1)=pde.left_boundary_Dirichlet;
    %F(N+1)=F(N+1)-pde.right_boundary_Neumann*basis_function_l(pde.end_point,x(N),h);
elseif boundary_method==2
    A(N+1,:)=[];
    A(:,N+1)=[];
    F(N+1)=[];
    F(1)=F(1)-pde.left_boundary_Neumann*basis_function_r(pde.start_point,x(2),h);
else
    A(N+1,N+1)=10^15;
    F(N+1)=10^15*pde.right_boundary_Dirichlet;
    F(1)=F(1)-pde.left_boundary_Neumann*basis_function_r(pde.start_point,x(2),h);
end  %neumann边值条件只有一种处理方式, 可以写在if外面
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

function [L2_error,H1_error]=errorEstimate(u,x,N,h,positions,weights)
%% Error estimate
L2_error=0;
H1_error=0;
inter_point=get_integration_points(0, h, positions);
b_values_l=basis_l(inter_point, 0, h);
b_values_r=basis_r(inter_point, 0, h);
b_values_m=basis_m(inter_point, 0, h);
b_values_l_der=basis_l_der(inter_point, 0, h);
b_values_r_der=basis_r_der(inter_point, 0, h);
b_values_m_der=basis_m_der(inter_point, 0, h);
%figure(4)
%hold on;
for i=1:N
    inter_point=get_integration_points(x(i),x(i+1),positions);
    value_ture=ture_sol(inter_point);
    value_ture_derivatives=ture_sol_der(inter_point);
    value_num=u(2*i-1)*b_values_r+u(2*i)*b_values_m+u(2*i+1)*b_values_l;
    value_num_derivatives=u(2*i-1)*b_values_r_der+u(2*i)*b_values_m_der+u(2*i+1)*b_values_l_der;
    values=(value_num-value_ture).^2;
    values_derivatives=(value_num_derivatives-value_ture_derivatives).^2;
    L2_error=L2_error+Gauss_Quadrature_formula(h,values,weights);
    % H1模直接用导数L2模即可
    H1_error=H1_error+Gauss_Quadrature_formula(h,values_derivatives,weights);
    
    %scatter(inter_point,value_num_derivatives,'filled','b');
end
L2_error=sqrt(L2_error);
H1_error=sqrt(H1_error);
end

function plotFigure(x,u,pde)
%% Plot figure
figure(1)
x_ture=pde.start_point:0.01:pde.end_point;
ture_u=ture_sol(x_ture);
plot(x_ture,ture_u,'k','linewidth',2);
hold on;
scatter(x,u,'filled','SizeData', 50);
title('Numerical Solution (N=10)')
end

function y=ture_sol(x)
%% Exact solution
y=sin(pi*x/2);
end

function y=ture_sol_der(x)
%% Derivatives of exaxt solution
y=pi/2*cos(pi*x/2);
end

function y=f(x)
%% Right hand function 
y=pi^2/2*sin(pi*x/2);
end

function  y=basis_l(x, left, right)
%% Left of basis function
%区间[left, right]上，以left和left+0.5h为零点的基函数
h=right-left;
y=2*((x-left)./h).*((x-left)./h-0.5);
end

function y=basis_l_der(x, left, right)
%% Derivatives of left of basis function
h=right-left;
y=4*(x-left)./h^2-1/h;
end

function  y=basis_r(x, left, right)
%% Right of basis function
%区间[left, right]上，以left+0.5h和left+h为零点的基函数
h=right-left;
y=2*((x-left)./h-1).*((x-left)./h-0.5);
end

function y=basis_r_der(x, left, right)
%% Derivatives of right of basis function
h=right-left;
y=4*(x-left)./h^2-3/h;
end

function  y=basis_m(x, left, right)
%% Right of basis function
%区间[left, right]上，以left和right为零点的基函数
h=right-left;
% y=4*((x-left)./h).*(1-(x-left)./h);
y=4*((x-left)./h).*((right-x)./h);
end

function y=basis_m_der(x, left, right)
%% Derivatives of middle of basis function
h=right-left;
y=4*(1/h-2*(x-left)./h^2);
end

function y=Gauss_Quadrature_formula(h,values,weights)
%% Gauss Quadrature formula
y=0.5*h*values*weights;
end

function y=get_integration_points(x1,x2,positions)
%% Get integration points
h=x2-x1;
y=h*0.5*positions+repmat((x1+x2)*0.5,1,5);
end

function [L2_convergenceOrder,H1_convergenceOrder]=convergenceOrder(N,pde)
%% Get integration points
L2_error_array=zeros(N,1);
H1_error_array=zeros(N,1);
Mat_con_array=zeros(N,1);
L2_convergenceOrder=zeros(N-1,1);
H1_convergenceOrder=zeros(N-1,1);
Mat_convergenceOrder=zeros(N-1,1);
N_vector=10*(1:N);

num_freedom=2*10*(1:N-1);
for i=1:N
    [L2_error,H1_error, Mat_con]=run_main(N_vector(i),pde);
    L2_error_array(i)=L2_error;
    H1_error_array(i)=H1_error;
    Mat_con_array(i)=Mat_con;
end
for i=1:N-1
    L2_convergenceOrder(i)=-log2(L2_error_array(i)/L2_error_array(i+1))/log2(i/(i+1));
    H1_convergenceOrder(i)=-log2(H1_error_array(i)/H1_error_array(i+1))/log2(i/(i+1));
    Mat_convergenceOrder(i)=log2(Mat_con_array(i)/Mat_con_array(i+1))/log2(i/(i+1));
end
figure(2)
%
subplot(1,3,1)
plot(N_vector,H1_error_array,'-o',N_vector,L2_error_array,'-*')
hleg1=legend('H1-norm','L2-norm','Location','NorthEast');
title('Error Estimation')
subplot(1,3,2)
loglog(N_vector,H1_error_array,'-o',N_vector,L2_error_array,'-*')
hleg1=legend('H1-norm','L2-norm','Location','NorthEast');
title('Error Estimation(loglog)')
subplot(1,3,3)
plot(num_freedom,H1_convergenceOrder,'-o',num_freedom,L2_convergenceOrder,'-*')
hleg1=legend('H1-norm','L2-norm','Location','NorthEast');
title('Convergence Order')
%plot(num_freedom,H1_error_array,'-o',num_freedom,exp(-1*log(num_freedom)+1.5),'',num_freedom,L2_error_array,'-*',num_freedom,exp(-2*log(num_freedom)+1),'-.')
%xlabel('degree of freedom in log scale');
%ylabel('error in log scale');



figure(3) %
subplot(1,2,1)
loglog(N_vector, Mat_con_array, '-*')
title('Matrix Condition','FontSize', 14)
subplot(1,2,2)
plot(num_freedom, Mat_convergenceOrder, 'm-*')
title('Divergence Order','FontSize', 14)

%axis([1e2,1e4,1e-4,1e1])
end
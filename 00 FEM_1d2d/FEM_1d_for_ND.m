function FEM_1d_for_ND
%%
%-u"=sinx   in [0,2*pi],
%u'(0)=1, u(2*pi)=0;
%
tic
format long;
pde.left_boundary_Neumann=1;
pde.right_boundary_Dirichlet=0;
pde.start_point=0;
pde.end_point=2*pi;

 %n_iter=13;
 %[L2_convergenceOrder,H1_convergenceOrder]=convergenceOrder(n_iter,pde)

N=1000;
[L2_error,H1_error]=run_main(N,pde)
toc
end

function [L2_error,H1_error]=run_main(N,pde)
%%
% 处理边值条件的方法
boundary_method=1;  % 1  2  3
Length=pde.end_point-pde.start_point;
h=Length/N;
x=(pde.start_point:h:pde.end_point)';

% [-1,1]区间上五点Guass积分公式的取点和权重
positions=[-0.906179845938664,-0.538469310105683,0,0.538469310105683,0.906179845938664];
weights=[0.236926885056189;0.478628670499366;0.568888888888889;0.478628670499366;0.236926885056189];

% 生成刚度矩阵
A=stiffnessMatrix(N,h);

% 生成右端项
F=rightHands(x,N,h,positions,weights);

% 处理边界条件
[A,F]=boundaryMatrix(A,F,x,N,h,pde,boundary_method);

% solve A\F
u=solveAF(A,F,N,pde,boundary_method);

% 计算误差
[L2_error,H1_error]=errorEstimate(u,x,N,h,positions,weights);

% plot figure
plotFigure(x,u,pde)
end

function A=stiffnessMatrix(N,h)
%% Generate stiffness matrix
A=sparse(N+1,N+1);
% 因为函数简单, 直接计算出数值即可
for i=1:N
    A(i,i)=A(i,i)+1/h;  
    A(i,i+1)=A(i,i+1)-1/h;
    A(i+1,i)=A(i+1,i)-1/h;
    A(i+1,i+1)=A(i+1,i+1)+1/h;   
end
end

function F=rightHands(x,N,h,positions,weights)
%% Generate right hands
F=zeros(N+1,1);
for i=1:N
    inter_point=get_integration_points(x(i),x(i+1),positions);
    f_values=f(inter_point);
    b_values_l=basis_function_l(inter_point,x(i),h);
    b_values_r=basis_function_r(inter_point,x(i+1),h);
    F(i)=F(i)+Gauss_Quadrature_formula(h,f_values.*b_values_r,weights);
    F(i+1)=F(i+1)+Gauss_Quadrature_formula(h,f_values.*b_values_l,weights);   
end
end

function [A,F]=boundaryMatrix(A,F,x,N,h,pde,boundary_method)
%% boundary condition
if boundary_method==1
    A(N+1,:)=0;
    A(N+1,N+1)=1;
    F(N+1)=pde.right_boundary_Dirichlet;
    F(1)=F(1)-pde.left_boundary_Neumann*basis_function_r(pde.start_point,x(2),h);
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
i1=1:N;
j1=2:N+1;
inter_point=get_integration_points(x(i1),x(j1),positions);
value_num=basis_function_l(inter_point,repmat(x(i1),1,length(positions)),h).*repmat(u(j1),1,length(positions))+...
    basis_function_r(inter_point,repmat(x(j1),1,length(positions)),h).*repmat(u(i1),1,length(positions));

value_num_derivatives=(1/h).*ones(size(inter_point)).*repmat(u(j1)-u(i1),1,length(positions));
value_ture=ture_sol(inter_point);
value_ture_derivatives=ture_sol_der(inter_point);
values=(value_num-value_ture).^2;
values_derivatives=(value_num_derivatives-value_ture_derivatives).^2;
L2_error=sqrt(sum(Gauss_Quadrature_formula(h,values,weights)));
H1_error=sqrt(sum(Gauss_Quadrature_formula(h,values_derivatives,weights)));  %H1半模
end

function plotFigure(x,u,pde)
%% Plot figure
x_ture=pde.start_point:0.01:pde.end_point;
ture_u=ture_sol(x_ture);
plot(x,u,'.',x_ture,ture_u)
end

function y=ture_sol(x)
%% Exact solution
y=sin(x);
end

function y=ture_sol_der(x)
%% Derivatives of exaxt solution
y=cos(x);
end

function y=f(x)
%% Right hand function 
y=sin(x);
end

function  y=basis_function_l(x,intervel_left,h)
%% Left of basis function
y=(x-intervel_left)./h;
end

function  y=basis_function_r(x,intervel_right,h)
%% Right of basis function
y=-(x-intervel_right)./h;
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
%% Compute the convergence order
L2_error_array=zeros(N,1);
H1_error_array=zeros(N,1);
L2_convergenceOrder=zeros(N-1,1);
H1_convergenceOrder=zeros(N-1,1);
N_vector=2.^(0:N-1).*10;

num_freedom=N_vector;
for i=1:N
    [L2_error,H1_error]=run_main(N_vector(i),pde);
    L2_error_array(i)=L2_error;
    H1_error_array(i)=H1_error;
end
for i=1:N-1
    L2_convergenceOrder(i)=log2(L2_error_array(i)/L2_error_array(i+1));
    H1_convergenceOrder(i)=log2(H1_error_array(i)/H1_error_array(i+1));
end
figure(2)
loglog(num_freedom,H1_error_array,'-o',num_freedom,exp(-1*log(num_freedom)+1.5),'',...
    num_freedom,L2_error_array,'-*',num_freedom,exp(-2*log(num_freedom)+1),'-.')
xlabel('degree of freedom in log scale');
ylabel('error in log scale');
hleg1=legend('H1-norm','slope of order 1','L2-norm','slope of order 2','Location','NorthEast');
title('Convergence Order')
%axis([1e2,1e4,1e-4,1e1])
end
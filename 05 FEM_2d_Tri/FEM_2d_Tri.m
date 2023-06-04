function FEM_2d_tri_ND_final
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

% ����������
 n_iter=10;
 [L2_convergenceOrder,H1_convergenceOrder]=convergenceOrder(n_iter,pde)

% ��������ͼ�� 
% N=10;
% [sDof, L2_error,H1_error]=run_main(N,pde)
toc
end

function [sDof, L2_error,H1_error]=run_main(N,pde)
%% intput N pde
%% output sDof,L2_error,H1_error, Mat_con
% �ֱ��ʾ���ɶ�, L2���, H1���͸նȾ���������

% �����ֵ�����ķ���
boundary_method=1;  % 1  2  3
domainLenght=pde.end_point_x-pde.start_point_x;
hmax=domainLenght/N;

% �����ʷ�
[p,e,t]=meshGeneration(pde,hmax);
sDof=size(p,2);

% ���ɸնȾ���
A=stiffnessMatrix(p,e,t);

% �����Ҷ���
F=rightHands(p,e,t);

% ����߽�����
[A,F]=boundaryMatrix(A,F,p,e,boundary_method);

% solve A\F
u=A\F;

% �������
[L2_error,H1_error]=errorEstimate(u,p,t);

% ��ͼ
plotFigure(p,e,t,u)
end

function [p,e,t]=meshGeneration(~,hmax)
%% Generate mesh
g=[ 2     2     2     2
    0     1     1    0
     1     1    0    0
    0    0     1     1
    0     1     1    0
     1     1     1     1
     0     0     0     0];
[p,e,t]=initmesh(g,'Hmax',hmax);
%figure(1)
%pdemesh(p,e,t,'NodeLabels', 'on');
% ��һ����y=0, �ڶ�����x=1, ��������y=1, ��������x=0
end

function A=stiffnessMatrix(p,~,t)
%% Generate stiffness matrix
sDof=size(p,2);     %�ڵ�����
% A=sparse(sDof,sDof);
A=zeros(sDof,sDof);
nel=size(t,2);         %С��������

for iel=1:nel
    nd=t(1:3,iel);
    a=p(:,nd(1));
    b=p(:,nd(2));
    c=p(:,nd(3));
    k=elK(a,b,c);
    A(nd,nd)=A(nd,nd)+k;
end
end

function F=rightHands(p,~,t)
%% Generate right hands
sDof=size(p,2);
F=zeros(sDof,1);   %F��ά�����ڱ������ɶ�
nel=size(t,2);         %�ʷֺ������θ���
for iel=1:nel
    nd=t(1:3,iel);
    a=p(:,nd(1));
    b=p(:,nd(2));
    c=p(:,nd(3));
    f=elF(a,b,c);
    F(nd)=F(nd)+f;
end
end   

function [A,F]=boundaryMatrix(A,F,p,e,boundary_method)
%% �����ֵ����
values_1=[0, 0];
values_2=[0, 0];
% ����Neumann��ֵ����
% ��������, y=1, ֻ��x�й�
loc_neuEdges_4=find(e(5,:)==3);
neuEdges_4=e(:,loc_neuEdges_4);
nneu_4=size(loc_neuEdges_4,2);
for iedge=1:nneu_4
    x1=p(1,neuEdges_4(1,iedge));
    x2=p(1,neuEdges_4(2,iedge));
    min_x=min(x1,x2);
    max_x=max(x1,x2);
    x=get_integration_points_1d(min_x,max_x);
    for k=1:2
        % ��ֵ
        values_1(k)=nboundaryFunction_4(x(k), 1) * (x(k)-x2)/(x1-x2);
    end
% �����Ҷ�����
gauss1=gaussint5(min_x,max_x,values_1);
F(neuEdges_4(1,iedge))=F(neuEdges_4(1,iedge))+gauss1;
    for k=1:2
          values_2(k)=nboundaryFunction_4(x(k), 1) * (x(k)-x1)/(x2-x1);
    end
    gauss2=gaussint5(min_x,max_x,values_2);
F(neuEdges_4(2,iedge))=F(neuEdges_4(2,iedge))+gauss2;
end

% �ڶ�����, x=1,ֻ��y�й�
loc_neuEdges_2=find(e(5,:)==2);
neuEdges_2=e(:,loc_neuEdges_2);
nneu_2=size(loc_neuEdges_2,2);
for iedge=1:nneu_2
    y1=p(2,neuEdges_2(1,iedge));
    y2=p(2,neuEdges_2(2,iedge));
    min_y=min(y1,y2);
    max_y=max(y1,y2);
    y=get_integration_points_1d(min_y,max_y);
    for k=1:2
        % ��ֵ
        values_1(k)=nboundaryFunction_2(1,y(k)) * (y(k)-y2)/(y1-y2);
    end
% �����Ҷ�����
F(neuEdges_2(1,iedge))=F(neuEdges_2(1,iedge))+gaussint5(min_y,max_y,values_1);
    for k=1:2
          values_2(k)=nboundaryFunction_2(1,y(k)) * (y(k)-y1)/(y2-y1);
    end
F(neuEdges_2(2,iedge))=F(neuEdges_2(2,iedge))+gaussint5(min_y,max_y,values_2);
end

% �����һ��ֵ����
% �ҵ�x=0����y=0�Ľڵ�
loc_d=find(p(1,:)==0|p(2,:)==0);
if boundary_method==1 
    x=p(1,loc_d);
    y=p(2,loc_d);
    A(loc_d,:)=0;
    A(loc_d,loc_d)=eye(size(loc_d,2));   
    F(loc_d)=boundaryFunction(x,y)';
end
end

function [L2_error,H1_error]=errorEstimate(u,p,t)
%% Error estimate
% L2=sqrt(/iint(u-u_h)^2dxdy)

%% Error estimate
% L2=sqrt(/iint(u-u_h)^2dxdy)

L2_error=0;
H1_half=0;
for i=1:size(t,2)
    triangle=t(1:3,i);
    a=p(:,triangle(1));
    b=p(:,triangle(2));
    c=p(:,triangle(3));
    [x,y]=get_integration_points(a,b,c);
    values_numerial=tri_basis(a,b,c,a,x,y).*u(triangle(1))+...
        tri_basis(a,b,c,b,x,y).*u(triangle(2))+...
        tri_basis(a,b,c,c,x,y).*u(triangle(3));
    values_true=true_sol(x,y);
    values=(values_numerial-values_true).^2;
    L2_error=L2_error+triangle_quadrature(a,b,c,values);
    
    values_numerial_gradient=tri_basis_gradient(a,b,c,a,x,y).*u(triangle(1))+...
        tri_basis_gradient(a,b,c,b,x,y).*u(triangle(2))+...
        tri_basis_gradient(a,b,c,c,x,y).*u(triangle(3));
    exactValues_gradient=ture_sol_gradient(x,y);
    values_gradient=sum((values_numerial_gradient-exactValues_gradient).^2,2);
    H1_half=H1_half+triangle_quadrature(a,b,c,values_gradient);    
end
L2_error=sqrt(L2_error);
H1_error=sqrt(L2_error^2+H1_half);
end

function plotFigure(p,e,t,u)
%% Plot figure
x=p(1,:)';
y=p(2,:)';
z=true_sol(x,y);

figure(2)
subplot(1,2,1);
%pdemesh(p,e,t,u);
pdeplot(p,e,t,"XYData",u,"ZData",u, ...
                  "FaceAlpha",0.5, ...
                  "ColorMap","jet", ...
                  "Mesh","on");
title('numerical solution');
subplot(1,2,2);
%pdemesh(p,e,t,z);
pdeplot(p,e,t,"XYData",z,"ZData",z, ...
                  "FaceAlpha",0.5, ...
                  "ColorMap","jet", ...
                  "Mesh","on");
title('exact solution');
end

function k_matrix=elK(a,b,c)
%% Generate local stiffness matrix
% ��basis_j'*��basis_k-2pi^2*basis_j.*basis_k
% ����������ȡ�߸���, ��ֵ���ֵõ�aij^k
k_matrix=zeros(3,3);
w=[a b c];
[x,y]=get_integration_points(a,b,c);
for j=1:3
    gradient_values_j=tri_basis_gradient(a,b,c,w(:,j),x,y);     
    %gradient_values 2*7�ľ���, ��һ��y������, �ڶ���x����
    for k=1:3
        gradient_values_k=tri_basis_gradient(a,b,c,w(:,k),x,y);
        s1=sum(gradient_values_j.*gradient_values_k,2);
        s2=-2*pi^2*tri_basis(a,b,c,w(:,j),x,y).*tri_basis(a,b,c,w(:,k),x,y);
        k_matrix(j,k)=triangle_quadrature(a,b,c,s1+s2);
    end
end
end

function k_f=elF(a,b,c)
%% Local right hands
%/iint fv dxdy
 k_f=zeros(3,1);
w=[a b c];
[x,y]=get_integration_points(a,b,c);
for j=1:3
    basis_values_j=tri_basis(a,b,c,w(:,j),x,y);
    right_values=rightFunction(x,y);
    values=sum(basis_values_j.*right_values,2);
    k_f(j)=triangle_quadrature(a,b,c,values);
end
end

function z=true_sol(x,y)
%% Exact solution
z=x.*y+sin(pi*x).*sin(pi*y);
end

function z=ture_sol_gradient(x,y)
%% Derivatives of exaxt solution
z_x=y+pi*cos(pi*x).*sin(pi*y);
z_y=x+pi*cos(pi*y).*sin(pi*x);
z=[z_x,z_y];
end

function z=rightFunction(x,y)
%% ƫ΢�ַ����Ҷ˺���
z=-2*pi*pi*x.*y;
end

function z=boundaryFunction(~,~)
%% Boundary function
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

function z=tri_basis(a,b,c,point,x,y)
%% Basis function
x1=a(1); y1=a(2); 
x2=b(1); y2=b(2);
x3=c(1); y3=c(2);
S=0.5*(x2*y3+x1*y2+x3*y1-x2*y1-x1*y3-x3*y2);
if point==a
    Sxy=0.5*(x2*y3+x*y2+x3*y-x2*y-x*y3-x3*y2);
elseif point==b
    Sxy=0.5*(x*y3+x1*y+x3*y1-x*y1-x1*y3-x3*y);
else
    Sxy=0.5*(x2*y+x1*y2+x*y1-x2*y1-x1*y-x*y2);
end
z=Sxy/S;
end

function z=tri_basis_gradient(a,b,c,point,x,y)
%% Gradient function
x1=a(1); y1=a(2); 
x2=b(1); y2=b(2);
x3=c(1); y3=c(2);
S=0.5*(x2*y3+x1*y2+x3*y1-x2*y1-x1*y3-x3*y2);
z=zeros((length(x)+length(y))/2,2);    %(length(x)+length(y))/2 ����ǿ��֤, ʹ��length(x)����
w=[a b c];
x_vector=zeros(2,1);
y_vector=zeros(2,1);
if point==a
    num=1;
elseif point==b
    num=2;
else
    num=3;
end
k=mod(num,3)+1;
x_vector(1)=w(1,k);
y_vector(1)=w(2,k);
k=k+1;
if(k>3)
    k=1;
end
x_vector(2)=w(1,k);
y_vector(2)=w(2,k);
z(:,1)=(y_vector(1)-y_vector(2))/S/2;
z(:,2)=(x_vector(2)-x_vector(1))/S/2;
end

function [x,y]=get_integration_points(a,b,c)
%% Generate intergration points
x_temp=[a(1);b(1);c(1)];
y_temp=[a(2);b(2);c(2)];
x=[x_temp;sum(x_temp)/3;sum(x_temp([1,2]))/2;sum(x_temp([2,3]))/2;sum(x_temp([1,3]))/2];
y=[y_temp;sum(y_temp)/3;sum(y_temp([1,2]))/2;sum(y_temp([2,3]))/2;sum(y_temp([1,3]))/2];
end

function x=get_integration_points_1d(a,b)
%% Generate intergration points(1D)
% location=[0.046910077030668 0.230765344947158 0.5  0.769234655052842 0.953089922969332];
% for i=1:5
% x(i)=a+location(i)*(b-a);
% end
x=[(b-a)*(1/2-3^(1/2)/6)+a,(b-a)*(1/2+3^(1/2)/6)+a];
end

function y=gaussint5(a,b,values)
% weights=[0.236926885056189 1.47863 0.568888888888889 0.478628670499366 0.236926885056189];
weights=[1 1];
y=(weights*values')*((b-a)/2);
end

function y = triangle_quadrature(a,b,c,values)
%% Quadrature formula
weights = [1/20;1/20;1/20;9/20;2/15;2/15;2/15];
x1=a(1); y1=a(2); 
x2=b(1); y2=b(2);
x3=c(1); y3=c(2);
S=0.5*(x2*y3+x1*y2+x3*y1-x2*y1-x1*y3-x3*y2);
y=S*(weights'*values);
end

function [L2_convergenceOrder,H1_convergenceOrder]=convergenceOrder(n,pde)
%% Get integration points
L2_error_array=zeros(n,1);
H1_error_array=zeros(n,1);
num_freedom_array=zeros(n,1);
L2_convergenceOrder=zeros(n-1,1);
H1_convergenceOrder=zeros(n-1,1);
N_vector=5*(1:n);

for i=1:n
    [sDof, L2_error, H1_error]=run_main(N_vector(i),pde);
    num_freedom_array(i)=sDof;
    L2_error_array(i)=L2_error;
    H1_error_array(i)=H1_error;
end
for i=1:n-1
    L2_convergenceOrder(i)=-log2(L2_error_array(i)/L2_error_array(i+1))/log2(i/(i+1));
    H1_convergenceOrder(i)=-log2(H1_error_array(i)/H1_error_array(i+1))/log2(i/(i+1));
end

figure(3)
subplot(1,3,1)
plot(num_freedom_array,H1_error_array,'-o',num_freedom_array,L2_error_array,'-*')
legend('H1-norm','L2-norm','Location','NorthEast');
title('Error Estimation')
subplot(1,3,2)
loglog(num_freedom_array,H1_error_array,'-o',num_freedom_array,L2_error_array,'-*')
legend('H1-norm','L2-norm','Location','NorthEast');
title('Error Estimation(loglog)')
subplot(1,3,3)
plot(num_freedom_array(1:n-1),H1_convergenceOrder,'-o',num_freedom_array(1:n-1),L2_convergenceOrder,'-*')
legend('H1-norm','L2-norm','Location','NorthEast');
title('Convergence Order')
end


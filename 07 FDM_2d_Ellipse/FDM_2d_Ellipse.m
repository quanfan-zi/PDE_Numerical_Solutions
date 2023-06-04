function FDM_2d_DN
%%
% u*=xy+sin(pi*x)sin(pi*y)   in [0,1]^2,
% -\Laplace u-2*pi^2*u=-2pi^2*xy
% u(0,y)=u(x,0)=0
% u_x(1,y)=y-pi*sin(pi*y)
% u_y(x,1)=x-pi*sin(pi*x)
tic
format short;
clear;clc;
close all;
pde.start_point_x=0;
pde.end_point_x=1;
pde.start_point_y=0;
pde.end_point_y=1;
% �����ֵ�����ķ���, 1��ʾ���̷�, 2��ʾ�������
pde.boundary_method=2;

% ����������
 n_iter=10;
 [Max_convergenceOrder,L2_convergenceOrder]=convergenceOrder(n_iter,pde)

% ��������ͼ�� 
% N=8;
% [sDof, Max_error, L2_error]=run_main(N,pde)
toc
end

function  [sDof, Max_error, L2_error]=run_main(N,pde)
%% output sDof,L2_error,H1_error, Mat_con

% �����ֵ�����ķ���
% 1��ʾ���̷�, 2��ʾ�������
boundary_method=pde.boundary_method;
domainLenght=pde.end_point_x-pde.start_point_x;
h=domainLenght/N;

% �����ʷ�
[p,e,t]=meshGeneration(pde,N);
sDof=size(p,2);
% sDof��ʾ�ڵ�����, Ҳ�������ɶ�
% eֻ�����߽��ϵı�, �����б�ʾС�������ĸ����, ���е�һ����x=0, �ڶ�����x=1, ��������y=0, ��������y=1

% ���ɸնȾ�����Ҷ���
[A, F]=stiffnessMatrix(p,e,t,N,h);

% ����߽�����
[A, F]=boundaryMatrix(A,F,p,N,h,boundary_method);

% solve A\F
u=A\F;

% �������
[Max_error, L2_error]=errorEstimate(u,p,t,h);

% ��ͼ
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
% ��һ����x=0, �ڶ�����x=1, ��������y=0, ��������y=1
end

function [A, F]=stiffnessMatrix(p,~,~,N,h)
%% Generate stiffness matrix
sDof=size(p,2);     
% A=sparse(sDof,sDof);
A=zeros(sDof,sDof);
F=zeros(sDof,1);
q=2*pi^2;  % ƫ΢�ַ��̲���
elK=[1/h^2, 1/h^2, q-4/h^2, 1/h^2, 1/h^2];  % ��Ԫ�նȾ���?

%��ʵ��߿���ֱ���þ������д����, ��������...
for i=2:N
    for j=2:N
        point_num=(i-1)*(N+1)+j;
        x_temp=p(1,point_num);
        y_temp=p(2,point_num);
        F(point_num)=rightFunction(x_temp,y_temp);
        points_num=[point_num-(N+1),point_num-1,point_num,point_num+1,point_num+(N+1)];
        A(point_num, points_num)=elK;
    end
end
end

function [A,F]=boundaryMatrix(A,F,p,N,h,boundary_method)
%% �����ֵ����
q=2*pi^2;  % ƫ΢�ַ��̲���

% ����Neumann��ֵ����
% boundary_method=1Ϊ΢�̷�, 2Ϊ���������
if boundary_method==1
    elk=[1/h^2, q-3/h^2, 1/h^2, 1/h^2];
    for i=2:N
        % ��������y=1
         point_num=(i)*(N+1);
         points_num=[point_num-(N+1),point_num,point_num+(N+1), point_num-1];
         x_temp=p(1,point_num);
         F(point_num)=rightFunction(x_temp,1)-1/h*nboundaryFunction_4(x_temp,1);
         A(point_num,points_num)=elk;
         
         % �ڶ�����x=1
         point_num=N*(N+1)+i;
         points_num=[point_num-1,point_num,point_num+1, point_num-(N+1)];
         y_temp=p(2,point_num);
         F(point_num)=rightFunction(1,y_temp)-1/h*nboundaryFunction_2(1,y_temp);
         A(point_num,points_num)=elk;
    end
    % ���ϽǶ���
    elk=[1/h^2, q-2/h^2,1/h^2];
    point_num=(N+1)^2;
    points_num=[point_num-(N+1),point_num,point_num-1];
    A(point_num,points_num)=elk;
    DFvalue=1/h*(nboundaryFunction_4(1,1)+nboundaryFunction_2(1,1));
    F(point_num)=rightFunction(1,y_temp)-DFvalue;
    
else
    elk=[1/h^2/2, q/2-2/h^2, 1/h^2/2, 1/h^2];
    for i=2:N
        % ��������y=1
        point_num=(i)*(N+1);
         points_num=[point_num-(N+1),point_num,point_num+(N+1), point_num-1];
         x_temp=p(1,point_num);
         F(point_num)=rightFunction(x_temp,1)/2-1/h*nboundaryFunction_4(x_temp,1);
         A(point_num,points_num)=elk;
         
         % �ڶ�����x=1
         point_num=N*(N+1)+i;
         points_num=[point_num-1,point_num,point_num+1, point_num-(N+1)];
         y_temp=p(2,point_num);
         F(point_num)=rightFunction(1,y_temp)/2-1/h*nboundaryFunction_2(1,y_temp);
         A(point_num,points_num)=elk;
    end
    % ���ϽǶ���
    elk=[1/h^2/2, q/4-1/h^2, 1/h^2/2];
    point_num=(N+1)^2;
    points_num=[point_num-(N+1),point_num,point_num-1];
    A(point_num,points_num)=elk;
    DFvalue=1/(2*h)*(nboundaryFunction_4(1,1)+nboundaryFunction_2(1,1));
    F(point_num)=rightFunction(1,y_temp)/4-DFvalue;
end

% ����Dirichlet��ֵ����
% �ҵ�x=0����y=0�Ľڵ�
loc_d=find(p(1,:)==0|p(2,:)==0);
A(loc_d,loc_d)=eye(size(loc_d,2));
end

function [Max_error, L2_error]=errorEstimate(u,p,t,h)
%% �������ģ��L2ģ
x=p(1,:)';
y=p(2,:)';
u_true=true_sol(x,y);
u_err=u-u_true;

Max_error=max(abs(u_err));

L2_error=sum(u_err.^2)*h;
L2_error= L2_error^(1/2);
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
Z=X.*Y+sin(pi*X).*sin(pi*Y);
surf(X,Y,Z)
colormap(jet)
alpha(0.5)
view(-15,20)
colorbar
title('exact solution');
end

function z=true_sol(x,y)
%% Exact solution
z=x.*y+sin(pi*x).*sin(pi*y);
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

function z=rightFunction(x,y)
%% ƫ΢�ַ����Ҷ˺���
z=2*pi*pi*x.*y;
end

function [Max_Order, L2_Order]=convergenceOrder(n,pde)
L2_error_array=zeros(n,1);
Max_error_array=zeros(n,1);
sDof_array=zeros(n,1);
L2_Order=zeros(n-1,1);
Max_Order=zeros(n-1,1);
N_vector=10*(1:n);

for i=1:n
    [sDof, Max_error, L2_error]=run_main(N_vector(i),pde);
    L2_error_array(i)=L2_error;
    Max_error_array(i)=Max_error;
    sDof_array(i)=sDof;
end
for i=1:n-1
    L2_Order(i)=-log2(L2_error_array(i)/L2_error_array(i+1))/log2(i/(i+1));
    Max_Order(i)=-log2(Max_error_array(i)/Max_error_array(i+1))/log2(i/(i+1));
end

figure(2)
loglog(sDof_array, L2_error_array,'m-*','linewidth',2)
hold on
loglog(sDof_array, Max_error_array,'-x','linewidth',2)
loglog(sDof_array,sDof_array.^(-1),'k-d','linewidth',2)
loglog(sDof_array,sDof_array.^(-0.5),'b-d','linewidth',2)
hold off
%loglog(num_freedom_array,H1_error_array,'-o',num_freedom_array,L2_error_array,'-*')
%legend('0-norm','1-norm','C-norm','y=x^{-2}/200','Location','NorthEast');
legend('0-norm','C-norm','2�������ٶ�','1�������ٶ�','Location','NorthEast');
title('Error Estimation(loglog)','fontsize',14)

end

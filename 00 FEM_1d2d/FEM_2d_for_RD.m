function FEM_2d_for_RD
%% 
% u= cos(pi.*x).*cos(pi.*y);
% -▽*(e^c(x.^2+y.^2)▽u)+u=f(x,y), c=2.1
% top and down  Dirichlet boundary condtion
% left and right Robin boundary condtion
tic
format long;
pde.start_point_x=-1;
pde.end_point_x=1;
pde.start_point_y=-1;
pde.end_point_y=1;

%N=5;    %收敛阶
%[L2_convergenceOrder,H1_convergenceOrder]=convergenceOrder(pde,N)
% N_vector=5*2.^(0:N-1); 跑一次大概需要6-10秒

 M=30;
 [sDof,L2_error,H1_error]=run_main(M,pde)

toc
end

function [sDof,L2_error,H1_error]=run_main(M,pde)
%% intput M pde 维数
%% output sDof,L2_error,H1_error 维数

%使用第一种方法处理Dirichlet边值条件
domainLenght=pde.end_point_x-pde.start_point_x;
hmax=domainLenght/M;

% Generate mesh
[p,e,t]=meshGeneration(hmax);
% p是点, 每一列表示划分后点的x, y坐标
% t是三角形, 前三行表示三角形顶点在p中的指标(列数)
% e是边, 只包含边界上的边(因为处理边值条件与内部边无关), 前两行表示起点和终点在p中的指标 
sDof=size(p,2);   %统计自由度

% Generate stiffness matrix
A=stiffnessMatrix(p,e,t);

% Generate right hands
F=rightHands(p,e,t);

% boundary condition
% 处理Dirichlet边值条件只用了第一种方法
 [A,F]=boundaryMatrix(A,F,p,e);

% solve A\F
u=A\F;

% error estimate
[L2_error,H1_error]=errorEstimate(u,p,t);

% plot figure
plotFigure(p,e,t,u);

end

function [p,e,t]=meshGeneration(hmax)
%% Generate mesh
g=[ 2     2     2     2
    -1     1     1    -1
     1     1    -1    -1
    -1    -1     1     1
    -1     1     1    -1
     1     1     1     1
     0     0     0     0];
[p,e,t]=initmesh(g,'Hmax',hmax);
end

function A=stiffnessMatrix(p,e,t)
%% Generate stiffness matrix
sDof=size(p,2);
A=sparse(sDof,sDof);
nel=size(t,2);          %nel表示网格剖分后三角形的个数
for iel=1:nel
    nd=t(1:3,iel);      %nd表示三角形的三个顶点
    a=p(:,nd(1));
    b=p(:,nd(2));
    c=p(:,nd(3));
    k=elK(a,b,c);
    A(nd,nd)=A(nd,nd)+k;
end

%加上robin边值条件的部分
%第二条边上
loc_robEdges_2=find(p(1,e(1,:))==1&p(1,e(2,:))==1);   
%通过find函数找到e中起点和终点都在x=1上的线段, 这样就在第二条边上了
robEdges_2=e(:,loc_robEdges_2);
nrob_2=size(robEdges_2,2);
for iedge=1:nrob_2
    %每条小边上的积分之和区间长度有关
    nd=robEdges_2(1:2,iedge);         %区间的端点指标
    y1=p(2,nd(1)); y2=p(2,nd(2));
    edgelength=abs(y1-y2);              %区间的长度
    k=[2, 1; 1, 2]*edgelength/6;
    %k=[2, 1; 1, 2]./(6*edgelength);
    A(nd,nd)=A(nd,nd)+k;
end

%第四条边上
loc_robEdges_4=find(p(1,e(1,:))==-1&p(1,e(2,:))==-1);   
robEdges_4=e(:,loc_robEdges_4);
nrob_4=size(robEdges_4,2);
for iedge=1:nrob_4
    nd=robEdges_4(1:2,iedge);         %区间的端点指标
    y1=p(2,nd(1)); y2=p(2,nd(2));
    edgelength=abs(y1-y2);              %区间的长度
    k=[2, 1; 1, 2]*edgelength/6;
    A(nd,nd)=A(nd,nd)+k;
end
end

function F=rightHands(p,~,t)
%% Generate right hands
sDof=size(p,2);
F=zeros(sDof,1);   %F的维数等于变量自由度
nel=size(t,2);         %剖分后三角形个数
for iel=1:nel
    nd=t(1:3,iel);
    a=p(:,nd(1));
    b=p(:,nd(2));
    c=p(:,nd(3));
    f=elF(a,b,c);
    F(nd)=F(nd)+f;
end
end    

function [A,F]=boundaryMatrix(A,F,p,e)
%% boundary condition
values_1=zeros(5,1);
values_2=zeros(5,1);
loc_neuEdges_2=find(p(1,e(1,:))==1&p(1,e(2,:))==1);   %修改
neuEdges_2=e(:,loc_neuEdges_2);
nneu_2=size(neuEdges_2,2);
for iedge=1:nneu_2                          %第二条边
    y1=p(2,neuEdges_2(1,iedge));
    y2=p(2,neuEdges_2(2,iedge));
    min_y=min(y1,y2);
    max_y=max(y1,y2);
    y=get_integration_points_1d(min_y,max_y);
    for k=1:5
          values_1(k)=rboundaryFunction_2(1,y(k))*(y(k)-y2)/(y1-y2);
    end
    F(neuEdges_2(1,iedge))=F(neuEdges_2(1,iedge))+gaussint5(min_y,max_y,values_1);
    for k=1:5
          values_2(k)=rboundaryFunction_2(1,y(k))*(y(k)-y1)/(y2-y1);
    end
    F(neuEdges_2(2,iedge))=F(neuEdges_2(2,iedge))+gaussint5(min_y,max_y,values_2);
end

loc_neuEdges_4=find(p(1,e(1,:))==-1&p(1,e(2,:))==-1);    %修改
neuEdges_4=e(:,loc_neuEdges_4);
nneu_4=size(neuEdges_4,2);
for iedge=1:nneu_4                     %第四条边
    y1=p(2,neuEdges_4(1,iedge));
    y2=p(2,neuEdges_4(2,iedge));
    min_y=min(y1,y2);
    max_y=max(y1,y2);
    y=get_integration_points_1d(min_y,max_y);
    for k=1:5
          values_1(k)=rboundaryFunction_4(-1,y(k))*(y(k)-y2)/(y1-y2);
    end
    F(neuEdges_4(1,iedge))=F(neuEdges_4(1,iedge))+gaussint5(min_y,max_y,values_1);
    for k=1:5
          values_2(k)=rboundaryFunction_4(-1,y(k))*(y(k)-y1)/(y2-y1);
    end
    F(neuEdges_4(2,iedge))=F(neuEdges_4(2,iedge))+gaussint5(min_y,max_y,values_2);  
end

% 处理Dirichlet边值条件
loc_d=find(p(2,:)==-1|p(2,:)==1);
    x=p(1,loc_d);
    y=p(2,loc_d);
    A(loc_d,:)=0;
    A(loc_d,loc_d)=eye(size(loc_d,2));    
    F(loc_d)=dboundaryFunction(x,y)';
end

function [L2_error,H1_error]=errorEstimate(u,p,t)
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
    values_true=cos(pi*x).*cos(pi*y);
    values=(values_numerial-values_true).^2;
    L2_error=L2_error+triangle_quadrature(a,b,c,values);
    
    values_numerial_gradient=tri_basis_gradient(a,b,c,a,x,y).*u(triangle(1))+...
        tri_basis_gradient(a,b,c,b,x,y).*u(triangle(2))+...
        tri_basis_gradient(a,b,c,c,x,y).*u(triangle(3));
    exactValues_gradient=[-pi.*sin(pi.*x).*cos(pi.*y),-pi.*cos(pi.*x).*sin(pi.*y)];
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
z=cos(pi*x).*cos(pi*y);
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
%exp(2.1*(x.^2+y.^2)).*▽basis_j'*▽basis_k+basis_j.*basis_k
%对上述函数取七个点, 数值积分得到aij^k
k_matrix=zeros(3,3);
w=[a b c];
[x,y]=get_integration_points(a,b,c);
for j=1:3
    gradient_values_j=tri_basis_gradient(a,b,c,w(:,j),x,y);     
    %gradient_values 2*7的矩阵, 第一行y方向导数, 第二行x方向
    for k=1:3
        gradient_values_k=tri_basis_gradient(a,b,c,w(:,k),x,y);
        gradient_values=sum(gradient_values_j.*gradient_values_k,2);
        s1=gradient_values.*exp(2.1*(x.^2+y.^2));
        s2=tri_basis(a,b,c,w(:,j),x,y).*tri_basis(a,b,c,w(:,k),x,y);
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

function z=rightFunction(x,y)
%% 偏微分方程右端项
s1=(pi*pi*cos(pi*x)+2*2.1*pi.*x.*sin(pi*x)).*exp(2.1*(x.^2+y.^2)).*cos(pi*y);
s2=(pi*pi*cos(pi*y)+2*2.1*pi.*y.*sin(pi*y)).*exp(2.1*(x.^2+y.^2)).*cos(pi*x);
s3=cos(pi.*x).*cos(pi.*y);
z=s1+s2+s3;
end

function z=dboundaryFunction(x,y)
%% Dirichlet边值条件
%输入和输出都是行向量
z=cos(pi.*x).*cos(pi.*y);
end

function z=rboundaryFunction_2(x,y)
%% Robin边值条件, 在第二条边
z=-exp(2.1)*pi*sin(pi*x).*cos(pi*y)+cos(pi.*x).*cos(pi.*y);
end

function z=rboundaryFunction_4(x,y)
%% Robin边值条件, 在第四条边
z=exp(2.1)*pi*sin(pi*x).*cos(pi*y)+cos(pi.*x).*cos(pi.*y);
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
z=zeros((length(x)+length(y))/2,2);    %(length(x)+length(y))/2 纯属强迫证, 使用length(x)即可
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
location=[0.046910077030668 0.230765344947158 0.5  0.769234655052842 0.953089922969332];
x=zeros(5,1);
for i=1:5
x(i)=a+location(i)*(b-a);
end
end

function y=gaussint5(a,b,values)
weights=[0.236926885056189 0.478628670499366 0.568888888888889 0.478628670499366 0.236926885056189];
y=(weights*values)*((b-a)/2);
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

function [L2_convergenceOrder,H1_convergenceOrder]=convergenceOrder(pde,N)
%% Compute the convergence order
L2_error_array=zeros(N,1);
H1_error_array=zeros(N,1);
num_freedom_array=zeros(N,1);
L2_convergenceOrder=zeros(N-1,1);
H1_convergenceOrder=zeros(N-1,1);
N_vector=5*2.^(0:N-1);

for i=1:N
    [num_freedom,L2_error,H1_error]=run_main(N_vector(i),pde);
    num_freedom_array(i)=num_freedom;
    L2_error_array(i)=L2_error;
    H1_error_array(i)=H1_error;
end
for i=1:N-1
    L2_convergenceOrder(i)=log2(L2_error_array(i)/L2_error_array(i+1));
    H1_convergenceOrder(i)=log2(H1_error_array(i)/H1_error_array(i+1));
end
figure(2)
loglog(num_freedom_array,H1_error_array,'-o',num_freedom_array,exp(-1*log(sqrt(num_freedom_array))+1.9),'',...
    num_freedom_array,L2_error_array,'-*',num_freedom_array,exp(-2*log(sqrt(num_freedom_array))+1.5),'-.')
xlabel('degree of freedom in log scale');
ylabel('error in log scale');
legend('H1-norm','slope of order 1','L2-norm','slope of order 2','Location','SouthWest');
title('Convergence Order')
%axis([1e2,1e5,1e-4,1e1])
end
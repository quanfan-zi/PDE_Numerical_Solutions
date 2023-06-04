%% 测试方程
% u*=sin(pi*x)*cos(pi*y)*exp(-pi^2*t/8),    0 < x, y < 1, t > 0
% u_t=(u_xx+u_yy)/16,
% u(x,y,0)=sin(pi*x)*cos(pi*y)
% u(0,y,t)=u(1,y,t)=0
% u_y(x,0,t)=u_y(x,1,t)=0

function FDM_parabolic_pro
%% 参数函数, 画函数图像或求收敛阶
tic 
format short;    clear;  clc;  close all;
pde.start_point=0;  pde.end_point=1;            % 横轴区间
pde.start_time=0;  pde.end_time=1;               % 时间区间
pde.subdivision_xy=40;                                    % 横轴剖分数
pde.r=1;                                                             % 网比r=at/h^2
% pde.method=1表示ADI 格式, 2表示预估校正法, 3表示LOD 格式

% 画t=1时函数图像
pde.method=1;
run_main(pde)

toc
end




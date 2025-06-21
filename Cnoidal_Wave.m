clc; clear; close all;
set(0, 'DefaultAxesFontName', '宋体');
set(0, 'DefaultTextFontName', '宋体');
% 输入参数
H = 0.125;   % 波高 (m)
h = 0.4;     % 水深 (m)
T = 2;       % 周期 (s)
g=9.81;

% m_vals = linspace(0.01, 0.999, 1000);
% f_vals = zeros(size(m_vals));
% for i = 1:length(m_vals)
%     m = m_vals(i);
%     K = ellipticK(m);
%     E = ellipticE(m);
%     lambda = sqrt((16*h*m)/(3*H)) *h* K;
%     c = sqrt(g*h*(1 + H/(m*h)*(2 - m - 3*E/K)));
%     f_vals(i) = lambda/c - T;
% end
% figure(1);
% plot(m_vals, f_vals);
% xlabel('m'); 
% ylabel('f(m) = λ/c - T');
% grid on; 
% yline(0, '--r');

% 初始猜测
m0 = 0.999999999999;          % for T=5s  m=0.9998779         for T=2s  m=0.8747414
tol = 1e-12;  % 收敛容差
% 步骤1: 通过牛顿迭代法求解 Jacobi 参数 m (根据 Cho, 2003 公式(7,9,10))
% 公式来源: Cho, Y.-S. (2003). A note on estimation of the Jacobian elliptic
%           parameter in cnoidal wave theory. Ocean Engineering, 30, 1915-1922.
m = m0;
for iter = 1:100
    if m <= 0 || m >= 1
        fprintf('m went out of bounds: m = %.5f\n', m);
    end
    K = ellipticK(m); E = ellipticE(m);
    % T^2 = (16 * h^3 * m^2 * K^2) / (3 * g * H * [m * h + H * (2 - m - 3 * E/K)])
    Y = m*h + 2*H - m*H - 3*H*(E/K) - (16*h^3 * m^2 * K^2)/(3*g*H*T^2);
    m1 = 1-m;
    dYdm = h - H + (3*H)/(2*m*m1*K^2)*(m1*K^2 + E^2 - 2*m1*K*E)...
        - (16*h^3 * m * K)/(3*g* m1*H*T^2 )*(m1*K +E);
    
    m = m - Y / (dYdm);
    % fprintf('Estimated m = %.10f after %d iterations\n', m, iter);
    if abs(Y / (dYdm)) < tol
        break
    end
end
fprintf('Estimated m = %.10f after %d iterations\n', m, iter);

% m=0.9998779
%% 步骤2: 自由面公式 z(x,t) 公式2
% Compute trough elevation eta2 (so mean level = 0)
K = ellipticK(m); E = ellipticE(m);
c = sqrt(g*h*(1 + H/(m*h)*(2 - m - 3*E/K)));
lambda=c*T
lambda = sqrt((16*h*m)/(3*H)) *h* K
x = linspace(0, 20*lambda, 50000);

eta2 = -H*(E/K - (1-m))/m;
% Compute cn^2 at each x

[sn, cn, dn] = ellipj(2*K*(x)/lambda, m);
eta = eta2 + H*(cn.^2);
figure(2);
plot(x, eta, 'LineWidth',1.5);
xlabel('x (m)'); ylabel('\eta(x,0) (m)');
title('Cnoidal Wave Free Surface (t=0)');
grid on;

%%  步骤3: Wave‐Maker Displacement (RK4 Integration) 公式(3)
zt = h + eta2;   % trough height above bottom

dt = 0.01; tmax = 80;
t = 0:dt:tmax; xi = zeros(size(t));
for n = 1:length(t)-1
    tn = t(n); xin = xi(n);
    phi1 = 2*K*(xin - c*tn)/lambda;
    [~, cn1, ~] = ellipj(phi1, m);
    f1 = c*(1 - h/(zt + H*cn1.^2));

    phi2 = 2*K*((xin+0.5*dt*f1) - c*(tn+0.5*dt))/lambda;
    [~, cn2, ~] = ellipj(phi2, m);
    f2 = c*(1 - h/(zt + H*(cn2.^(2))));

    phi3 = 2*K*((xin+0.5*dt*f2) - c*(tn+0.5*dt))/lambda;
    [~, cn3, ~] = ellipj(phi3, m);
    f3 = c*(1 - h/(zt + H*(cn3.^(2))));

    phi4 = 2*K*((xin+dt*f3) - c*(tn+dt))/lambda;
    [~, cn4, ~] = ellipj(phi4, m);
    f4 = c*(1 - h/(zt + H*(cn4.^(2))));
    xi(n+1) = xin + (dt/6)*(f1 + 2*f2 + 2*f3 + f4);
end
figure(3);
plot(t, xi, 'LineWidth',1.5);hold on;

t_ramp = 1.5*T;         % 启动时间 [s]
% 应用启动渐进 (ramp-up)
ramp = (1 - cos(pi * t / t_ramp)) / 2;
ramp(t > t_ramp) = 1;
xii = xi .* ramp;
plot(t, xii, 'LineWidth',1.5);

xlabel('t (s)'); ylabel('Wave-maker displacement x (m)');
title('Wave-maker Displacement vs Time');
grid on;
legend('解析位移','带ramp-up位移');
% 保存摇板位移数据
output_data = [t' xii'];
writematrix(output_data, 'paddle_displacement_cnoidal.csv');
fprintf('\n摇板位移数据已保存至: paddle_displacement_cnoidal.csv\n');
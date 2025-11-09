% 模型混合函数
function [x_pre, P] = Model_mix_SingerAP_CS_jerk(x1, x2,  P1, P2,  u)
    % 综合状态估计
    x_pre = x1 * u(1) + x2 * u(2) ;
    
    % 综合协方差估计
    P = u(1)*(P1 + (x1 - x_pre)*(x1 - x_pre)') + ...
        u(2)*(P2 + (x2 - x_pre)*(x2 - x_pre)');
end
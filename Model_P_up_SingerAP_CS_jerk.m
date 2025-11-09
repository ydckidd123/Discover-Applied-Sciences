function [u,L] = Model_P_up_SingerAP_CS_jerk(r1, r2,  S1, S2, c_j)
       % 计算似然函数（多元高斯分布）
    L1 = exp(-0.5 * r1' / S1 * r1) / sqrt(det(2*pi*S1));
    L2 = exp(-0.5 * r2' / S2 * r2) / sqrt(det(2*pi*S2));
    
    
    % 归一化似然
    L_sum = L1 + L2 ;
    L1 = L1 / L_sum;
    L2 = L2 / L_sum;
    L=[L1;L2];
    
    % 更新模型概率
    u = [L1; L2] .* c_j;
    u = u / sum(u); % 归一化
end
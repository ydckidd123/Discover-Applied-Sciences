% 标准IMM-SingerAPCS-jerk-CKF算法
function [u_IMM]=IMM_CKF(x0_filter,z,R,time,T)
% 设置Singer模型参数
alpha1 =1/2;           % 机动频率（机动时间常数的倒数） 
sigma_a1 = 1;           % 加速度标准差
q_s1 = 2 * alpha1 * sigma_a1^2; % 连续时间过程噪声强度

% 计算Singer模型的一维状态转移矩阵
F1d1 = [1, T, (alpha1*T-1+exp(-alpha1*T))/(alpha1^2),0;
        0, 1, (1-exp(-alpha1*T))/alpha1,0;
        0, 0, exp(-alpha1*T),0;
        0, 0, 0, 0];
% 二维Singer模型
F_singer1 = blkdiag(F1d1, F1d1);

% 计算Singer模型的一维过程噪声协方差
q111 = (2*alpha1^3*T^3 - 6*alpha1^2*T^2 + 6*alpha1*T + 12*alpha1*T*exp(-alpha1*T) ...
       - 3*alpha1^2*T^2*exp(-alpha1*T) - 3*exp(-2*alpha1*T) + 3) / (6*alpha1^5);
q121 = (alpha1^2*T^2 - 2*alpha1*T + 2 - 2*exp(-alpha1*T) ...
       + 2*alpha1*T*exp(-alpha1*T) + alpha1^2*T^2*exp(-alpha1*T) + exp(-2*alpha1*T) - 1) / (2*alpha1^4);
q131 = (1 - 2*exp(-alpha1*T) + 2*alpha1*T*exp(-alpha1*T) + exp(-2*alpha1*T)) / (2*alpha1^3);
q221 = (2*alpha1*T - 3 + 4*exp(-alpha1*T) - exp(-2*alpha1*T)) / (2*alpha1^3);
q231 = (1 - 2*alpha1*T*exp(-alpha1*T) - exp(-2*alpha1*T)) / (2*alpha1^2);
q331 = (1 - exp(-2*alpha1*T)) / (2*alpha1);

Q1d1 = q_s1 * [q111, q121, q131, 0;
             q121, q221, q231, 0;
             q131, q231, q331, 0;
             0,   0,   0,   0];
% 二维过程噪声协方差
Q_singer1 = blkdiag(Q1d1, Q1d1);

% APCS-Jerk 模型参数
alpha_cs1 = 2;       % 机动频率 (1/秒)
alpha_cs01 = 2;       % 机动频率 (1/秒)
j_max1 = 5;            % 最大加加速度 (m/s^3)
j_max01 = 5;            % 最大加加速度 (m/s^3)

% Jerk模型状态转移矩阵 (8维状态: [x, vx, ax, jx, y, vy, ay, jy])
F_jerk1 = [1, T, T^2/2, T^3/6, 0, 0, 0, 0;
          0, 1, T,     T^2/2, 0, 0, 0, 0;
          0, 0, 1,     T,     0, 0, 0, 0;
          0, 0, 0,     1,     0, 0, 0, 0;
          0, 0, 0,     0,     1, T, T^2/2, T^3/6;
          0, 0, 0,     0,     0, 1, T,     T^2/2;
          0, 0, 0,     0,     0, 0, 1,     T;
          0, 0, 0,     0,     0, 0, 0,     1];
% Jerk模型过程噪声离散化矩阵 (CS模型基础矩阵)
q111=1/(2*alpha_cs1^7)*(alpha_cs1^5*T^5/10-alpha_cs1^4*T^4/2+4*alpha_cs1^3*T^3/3-2*alpha_cs1^2*T^2+2*alpha_cs1*T-3+4*exp(-alpha_cs1*T)+2*alpha_cs1^2*T^2*exp(-alpha_cs1*T)-exp(-2*alpha_cs1*T));
q121=1/(2*alpha_cs1^6)*(1-2*alpha_cs1*T+2*alpha_cs1^2*T^2-alpha_cs1^3*T^3+alpha_cs1^4*T^4/4+exp(-2*alpha_cs1*T)+2*alpha_cs1*T*exp(-alpha_cs1*T)-2*exp(-alpha_cs1*T)-alpha_cs1^2*T^2*exp(-alpha_cs1*T));
q131=1/(2*alpha_cs1^5)*(2*alpha_cs1*T-alpha_cs1^2*T^2-alpha_cs1^3*T^3/3-3-2*exp(-2*alpha_cs1*T)+4*exp(-alpha_cs1*T)+alpha_cs1^2*T^2*exp(-alpha_cs1*T));
q141=1/(2*alpha_cs1^4)*(1+exp(-2*alpha_cs1*T)-2*exp(-alpha_cs1*T)-alpha_cs1^2*T^2*exp(-alpha_cs1*T));
q221=1/(2*alpha_cs1^5)*(1-exp(-2*alpha_cs1*T)+2*alpha_cs1^3*T^3/3+2*alpha_cs1*T-2*alpha_cs1^2*T^2-4*alpha_cs1*T*exp(-alpha_cs1*T));
q231=1/(2*alpha_cs1^4)*(1+alpha_cs1^2*T^2-2*alpha_cs1*T+2*alpha_cs1*T*exp(-alpha_cs1*T)+exp(-2*alpha_cs1*T)-2*exp(-alpha_cs1*T));
q241=1/(2*alpha_cs1^3)*(1-exp(-2*alpha_cs1*T)-2*alpha_cs1*T*exp(-alpha_cs1*T));
q331=1/(2*alpha_cs1^3)*(4*exp(-alpha_cs1*T)-exp(-2*alpha_cs1*T)+2*alpha_cs1*T-3);
q341=1/(2*alpha_cs1^2)*(1-2*exp(-alpha_cs1*T)+exp(-2*alpha_cs1*T));
q441=1/(2*alpha_cs1)*(1-exp(-2*alpha_cs1*T));

% 基础矩阵 (不含噪声强度)
Q_base01 = [q111, q121, q131, q141;
           q121, q221, q231, q241;
           q131, q231, q331, q341;
           q141, q241, q341, q441];
G1 = [T^2/2, 0;
     T, 0;
     1, 0;
     1, 0;
     0, T^2/2;
     0, T;
     0, 1;
     0, 1];  

%IMM-CKF初始化
X_IMM1 = zeros(8, time); 
P_IMM1 = zeros(8, 8, time); 
pij11 = [0.7, 0.3;
       0.3, 0.7];    
u_IMM = zeros(2, time);   % 模型概率
u_IMM(:,1) = [0.6, 0.4]'; 
P011 = diag([1000, 500, 100,1, 1000, 500, 100,1]); % 初始协方差

% 各模型滤波器初始化
X_Singer1 = x0_filter; X_APCSJerk1 = x0_filter; 

P_Singer1 = P011; P_APCSJerk1 = P011; 

X_IMM1(:,1) = x0_filter; 
P_IMM1(:,:,1) = P011;
% 记录各模型各时刻状态
X_Singer_01 = zeros(8, time); X_APCSJerk_01 = zeros(8, time);
X_Singer_01(:,1) = x0_filter; X_APCSJerk_01(:,1) = x0_filter; 
d=zeros(time,1);
D_APCS_Jerk1=zeros(time,1);
%% IMM-CKF迭代
for t1=1:time-1
    % 第一步: 交互混合
    c_j1 = pij11' * u_IMM(:,t1); % 归一化常数
    
    % 计算混合概率
    ui11 = (1/c_j1(1)) * pij11(:,1) .* u_IMM(:,t1);
    ui21 = (1/c_j1(2)) * pij11(:,2) .* u_IMM(:,t1);
    
   % 计算混合后的状态和协方差
    x011 = X_Singer1*ui11(1) + X_APCSJerk1*ui11(2) ;
    x021 = X_Singer1*ui21(1) + X_APCSJerk1*ui21(2) ;
    
    P011 = (P_Singer1 + (X_Singer1- x011)*(X_Singer1 - x011)')*ui11(1) + ...
          (P_APCSJerk1 + (X_APCSJerk1 - x011)*(X_APCSJerk1 - x011)')*ui11(2);
    
    P021 = (P_Singer1 + (X_Singer1 - x021)*(X_Singer1 - x021)')*ui21(1) + ...
          (P_APCSJerk1 + (X_APCSJerk1 - x021)*(X_APCSJerk1 - x021)')*ui21(2) ;

    % 第二步: 模型条件滤波（使用CKF代替UKF）
    [X_Singer1, P_Singer1, r_Singer1, S_Singer1] = CKF_Singer(x011, P011, z(:,t1+1), F_singer1,  Q_singer1, R);
    [X_APCSJerk1, P_APCSJerk1, r_APCSJerk1, S_APCSJerk1] = CKF_APCSJerk(x021, P021, z(:,t1+1), F_jerk1,  Q_base01, R,alpha_cs1,j_max1);
    d(t1)=r_APCSJerk1'/S_APCSJerk1*r_APCSJerk1;

   
    if t1>3
        D_APCS_Jerk1(t1)=(d(t1)+d(t1-1)+d(t1-2))/3;
    else
        D_APCS_Jerk1(t1)=d(t1);
    end

    if t1>3
       if D_APCS_Jerk1(t1)>9.21
          f1=exp(0.75*(D_APCS_Jerk1(t1)-D_APCS_Jerk1(t1-1)));
       else
          f1=1;
       end
    else
        f1=1;
    end

    % N=4
%     if t1>4
%         D_APCS_Jerk1(t1)=(d(t1)+d(t1-1)+d(t1-2)+d(t1-3))/4;
%     else
%         D_APCS_Jerk1(t1)=d(t1);
%     end
% 
%     if t1>3
%        if D_APCS_Jerk1(t1)>9.21
%           f1=exp(0.75*(D_APCS_Jerk1(t1)-D_APCS_Jerk1(t1-1)));
%        else
%           f1=1;
%        end
%     else
%         f1=1;
%     end

    % N=5
    if t1>5
        D_APCS_Jerk1(t1)=(d(t1)+d(t1-1)+d(t1-2)+d(t1-3)+d(t1-4))/5;
    else
        D_APCS_Jerk1(t1)=d(t1);
    end

    if t1>3
       if D_APCS_Jerk1(t1)>9.21
          f1=exp(0.75*(D_APCS_Jerk1(t1)-D_APCS_Jerk1(t1-1)));
       else
          f1=1;
       end
    else
        f1=1;
    end


   

    if f1>5
        f1=5;
    end
    j_max1=f1*j_max01;
    alpha_cs1=f1*alpha_cs01;


    % 第三步: 模型概率更新
    [u_IMM(:,t1+1)] = Model_P_up_SingerAP_CS_jerk(r_Singer1, r_APCSJerk1, S_Singer1, S_APCSJerk1, c_j1);

    % 第四步: 估计融合
    [X_IMM1(:,t1+1), P_IMM1(:,:,t1+1)] = Model_mix_SingerAP_CS_jerk(X_Singer1, X_APCSJerk1, P_Singer1, P_APCSJerk1, u_IMM(:,t1+1));

    % 保存各模型状态
    X_Singer_01(:, t1+1) = X_Singer1; 
    X_APCSJerk_01(:, t1+1) = X_APCSJerk1; 
    
end


end






end
%% ATPM-PIMM-SingerAPCS-jerk-RACKF
function[X_IMM,u_IMM,JC]=APIMM_RACKF(x0_filter,z,R,time,T,JC)
% ATPM预处理
K1=[0;0];
[u_IMM_CKF]=IMM_CKF(x0_filter,z,R,time,T);
index_match=zeros(1,time);
Lrl=zeros(1,1);


% 设置Singer模型参数
%T=0.5;
alpha =1/2;           % 机动频率（机动时间常数的倒数） 
sigma_a = 1;           % 加速度标准差
q_s = 2 * alpha * sigma_a^2; % 连续时间过程噪声强度

% 计算Singer模型的一维状态转移矩阵
F1d = [1, T, (alpha*T-1+exp(-alpha*T))/(alpha^2),0;
        0, 1, (1-exp(-alpha*T))/alpha,0;
        0, 0, exp(-alpha*T),0;
        0, 0, 0, 0];
% 二维Singer模型
F_singer = blkdiag(F1d, F1d);

% 计算Singer模型的一维过程噪声协方差
q11 = (1-exp(-2*alpha*T)+2*alpha*T+2*alpha^3*T^3/3-2*alpha^2*T^2-4*alpha*T*exp(-alpha*T)) / (2*alpha^5);
q12 = (alpha^2*T^2 - 2*alpha*T + 1 - 2*exp(-alpha*T) ...
       + 2*alpha*T*exp(-alpha*T)  + exp(-2*alpha*T) ) / (2*alpha^4);
q13 = (1 - 2*alpha*T*exp(-alpha*T) -2*alpha*T*exp(-alpha*T)) / (2*alpha^3);
q22 = (2*alpha*T - 3 + 4*exp(-alpha*T) - exp(-2*alpha*T)) / (2*alpha^3);
q23 = (1 - 2*exp(-alpha*T) + exp(-2*alpha*T)) / (2*alpha^2);
q33 = (1 - exp(-2*alpha*T)) / (2*alpha);

Q1d = q_s * [q11, q12, q13, 0;
             q12, q22, q23, 0;
             q13, q23, q33, 0;
             0,   0,   0,   0];
% 二维过程噪声协方差
Q_singer = blkdiag(Q1d, Q1d);

% APCS-Jerk 模型参数
alpha_cs = 4;       % 机动频率 (1/秒)
alpha_cs0 = 4;       % 机动频率 (1/秒)
j_max =10;            % 最大加加速度 (m/s^3)
j_max0 = 10;            % 最大加加速度 (m/s^3)

% Jerk模型状态转移矩阵 (8维状态: [x, vx, ax, jx, y, vy, ay, jy])
F_jerk = [1, T, T^2/2, T^3/6, 0, 0, 0, 0;
          0, 1, T,     T^2/2, 0, 0, 0, 0;
          0, 0, 1,     T,     0, 0, 0, 0;
          0, 0, 0,     1,     0, 0, 0, 0;
          0, 0, 0,     0,     1, T, T^2/2, T^3/6;
          0, 0, 0,     0,     0, 1, T,     T^2/2;
          0, 0, 0,     0,     0, 0, 1,     T;
          0, 0, 0,     0,     0, 0, 0,     1];
% Jerk模型过程噪声离散化矩阵 (CS模型基础矩阵)
q11=1/(2*alpha_cs^7)*(alpha_cs^5*T^5/10-alpha_cs^4*T^4/2+4*alpha_cs^3*T^3/3-2*alpha_cs^2*T^2+2*alpha_cs*T-3+4*exp(-alpha_cs*T)+2*alpha_cs^2*T^2*exp(-alpha_cs*T)-exp(-2*alpha_cs*T));
q12=1/(2*alpha_cs^6)*(1-2*alpha_cs*T+2*alpha_cs^2*T^2-alpha_cs^3*T^3+alpha_cs^4*T^4/4+exp(-2*alpha_cs*T)+2*alpha_cs*T*exp(-alpha_cs*T)-2*exp(-alpha_cs*T)-alpha_cs^2*T^2*exp(-alpha_cs*T));
q13=1/(2*alpha_cs^5)*(2*alpha_cs*T-alpha_cs^2*T^2-alpha_cs^3*T^3/3-3-2*exp(-2*alpha_cs*T)+4*exp(-alpha_cs*T)+alpha_cs^2*T^2*exp(-alpha_cs*T));
q14=1/(2*alpha_cs^4)*(1+exp(-2*alpha_cs*T)-2*exp(-alpha_cs*T)-alpha_cs^2*T^2*exp(-alpha_cs*T));
q22=1/(2*alpha_cs^5)*(1-exp(-2*alpha_cs*T)+2*alpha_cs^3*T^3/3+2*alpha_cs*T-2*alpha_cs^2*T^2-4*alpha_cs*T*exp(-alpha_cs*T));
q23=1/(2*alpha_cs^4)*(1+alpha_cs^2*T^2-2*alpha_cs*T+2*alpha_cs*T*exp(-alpha_cs*T)+exp(-2*alpha_cs*T)-2*exp(-alpha_cs*T));
q24=1/(2*alpha_cs^3)*(1-exp(-2*alpha_cs*T)-2*alpha_cs*T*exp(-alpha_cs*T));
q33=1/(2*alpha_cs^3)*(4*exp(-alpha_cs*T)-exp(-2*alpha_cs*T)+2*alpha_cs*T-3);
q34=1/(2*alpha_cs^2)*(1-2*exp(-alpha_cs*T)+exp(-2*alpha_cs*T));
q44=1/(2*alpha_cs)*(1-exp(-2*alpha_cs*T));

% 基础矩阵 (不含噪声强度)
Q_base0 = [q11, q12, q13, q14;
           q12, q22, q23, q24;
           q13, q23, q33, q34;
           q14, q24, q34, q44];
G = [T^2/2, 0;
     T, 0;
     1, 0;
     1, 0;
     0, T^2/2;
     0, T;
     0, 1;
     0, 1];  

%IMM-CKF初始化
X_IMM = zeros(8, time); 
P_IMM = zeros(8, 8, time); 
pij = [0.7, 0.3;
       0.3, 0.7];    
pij_Singer= [0.4,0.3,0.3;
             0.3,0.4,0.3;
             0.3,0.3,0.4];
pij_APCSjerk= [0.4,0.3,0.3;
               0.3,0.4,0.3;
               0.3,0.3,0.4];
u_IMM = zeros(2, time);   % 模型概率
u_Singer = zeros(3, time);
u_APCSjerk= zeros(3, time);

u_IMM(:,1) = [0.6, 0.4]'; 
u_Singer(:,1) =  [0.4,0.3,0.3]';
u_APCSjerk(:,1) =  [0.4,0.3,0.3]';

P0 = diag([1000, 500, 100,1, 1000, 500, 100,1]); % 初始协方差
% 各模型滤波器初始化
X_Singer = x0_filter; X_APCSJerk = x0_filter; 
X_Singer_CKF = x0_filter; X_APCSJerk_CKF = x0_filter; 
X_Singer_RCKF = x0_filter; X_APCSJerk_RCKF = x0_filter; 
X_Singer_ACKF = x0_filter; X_APCSJerk_ACKF = x0_filter; 

P_Singer = P0; P_APCSJerk = P0;
P_Singer_CKF = P0; P_APCSJerk_CKF = P0; 
P_Singer_RCKF = P0; P_APCSJerk_RCKF = P0;
P_Singer_ACKF = P0; P_APCSJerk_ACKF = P0;

X_IMM(:,1) = x0_filter; 
P_IMM(:,:,1) = P0;
% 记录各模型各时刻状态
X_Singer_0 = zeros(8, time); X_APCSJerk_0 = zeros(8, time);
X_Singer_0(:,1) = x0_filter; X_APCSJerk_0(:,1) = x0_filter; 
d1=zeros(time,1);
D_APCS_Jerk=zeros(time,1);
D_RCKF_Singer=zeros(time,1);D_ACKF_Singer=zeros(time,1);
D_RCKF_APCSJerk=zeros(time,1);D_ACKF_APCSJerk=zeros(time,1);
%% IMM-CKF迭代
for t=1:time-1
    % 第一步: 交互混合
    c_j = pij' * u_IMM(:,t); % 归一化常数
    
    % 计算混合概率
    ui1 = (1/c_j(1)) * pij(:,1) .* u_IMM(:,t);
    ui2 = (1/c_j(2)) * pij(:,2) .* u_IMM(:,t);
    
   % 计算混合后的状态和协方差
    x01_CKF = X_Singer_CKF*ui1(1) + X_APCSJerk_CKF*ui1(2) ;
    x01_RCKF = X_Singer_RCKF*ui1(1) + X_APCSJerk_RCKF*ui1(2) ;
    x01_ACKF = X_Singer_ACKF*ui1(1) + X_APCSJerk_ACKF*ui1(2) ;
    
    x02_CKF = X_Singer_CKF*ui2(1) + X_APCSJerk_CKF*ui2(2) ;
    x02_RCKF = X_Singer_RCKF*ui2(1) + X_APCSJerk_RCKF*ui2(2) ;
    x02_ACKF = X_Singer_ACKF*ui2(1) + X_APCSJerk_ACKF*ui2(2) ;

    P01_CKF = (P_Singer_CKF + (X_Singer_CKF - x01_CKF)*(X_Singer_CKF - x01_CKF)')*ui1(1) + ...
          (P_APCSJerk_CKF + (X_APCSJerk_CKF - x01_CKF)*(X_APCSJerk_CKF - x01_CKF)')*ui1(2);
    P01_RCKF = (P_Singer_RCKF + (X_Singer_RCKF - x01_RCKF)*(X_Singer_RCKF - x01_RCKF)')*ui1(1) + ...
          (P_APCSJerk_RCKF + (X_APCSJerk_RCKF - x01_RCKF)*(X_APCSJerk_RCKF - x01_RCKF)')*ui1(2);
    P01_ACKF = (P_Singer_ACKF + (X_Singer_ACKF - x01_ACKF)*(X_Singer_ACKF - x01_ACKF)')*ui1(1) + ...
          (P_APCSJerk_ACKF + (X_APCSJerk_ACKF - x01_ACKF)*(X_APCSJerk_ACKF - x01_ACKF)')*ui1(2);
    
    P02_CKF = (P_Singer_CKF + (X_Singer_CKF - x02_CKF)*(X_Singer_CKF - x02_CKF)')*ui2(1) + ...
          (P_APCSJerk_CKF + (X_APCSJerk_CKF - x02_CKF)*(X_APCSJerk_CKF - x02_CKF)')*ui2(2) ;
    P02_RCKF = (P_Singer_RCKF + (X_Singer_RCKF - x02_RCKF)*(X_Singer_RCKF - x02_RCKF)')*ui2(1) + ...
          (P_APCSJerk_RCKF + (X_APCSJerk_RCKF - x02_RCKF)*(X_APCSJerk_RCKF - x02_RCKF)')*ui2(2) ;
    P02_ACKF = (P_Singer_ACKF + (X_Singer_ACKF - x02_ACKF)*(X_Singer_ACKF - x02_ACKF)')*ui2(1) + ...
          (P_APCSJerk_ACKF + (X_APCSJerk_ACKF - x02_ACKF)*(X_APCSJerk_ACKF - x02_ACKF)')*ui2(2) ;

    % 第二步: 模型条件滤波
    %Singer
    [X_Singer, P_Singer, r_Singer, S_Singer,u_Singer(:,t+1),X_Singer_CKF,X_Singer_RCKF,X_Singer_ACKF,P_Singer_CKF,P_Singer_RCKF,P_Singer_ACKF,pij_Singer,D_RCKF_Singer,D_ACKF_Singer] =...
        RACKF_Singer(x01_CKF,x01_RCKF,x01_ACKF, P01_CKF,P01_RCKF,P01_ACKF, z(:,t+1), F_singer,  Q_singer, R,u_Singer(:,t),pij_Singer,t,D_RCKF_Singer,D_ACKF_Singer);
    

    %AP-CS-jerk
    [X_APCSJerk, P_APCSJerk, r_APCSJerk, S_APCSJerk,u_APCSjerk(:,t+1),X_APCSJerk_CKF,X_APCSJerk_RCKF,X_APCSJerk_ACKF,P_APCSJerk_CKF,P_APCSJerk_RCKF,P_APCSJerk_ACKF,pij_APCSjerk,D_RCKF_APCSJerk,D_ACKF_APCSJerk] =...
        RACKF_APCSJerk(x02_CKF,x02_RCKF,x02_ACKF, P02_CKF,P02_RCKF,P02_ACKF, z(:,t+1), F_jerk,  Q_base0, R,u_APCSjerk(:,t),pij_APCSjerk,alpha_cs,j_max,t,D_RCKF_APCSJerk,D_ACKF_APCSJerk);
    
    d1(t)=r_APCSJerk'/S_APCSJerk*r_APCSJerk;

   

    % N=3
    if t>3
        D_APCS_Jerk(t)=(d1(t)+d1(t-1)+d1(t-2))/3;
    else
        D_APCS_Jerk(t)=d1(t);
    end

    if t>3
       if D_APCS_Jerk(t)>5.991
          f=exp(0.75*(D_APCS_Jerk(t)-D_APCS_Jerk(t-1)));
          %f=exp(1*(d1(t)-d1(t-1)));
           %fprintf('%f',f);
       else
          f=1;
       end
    else
        f=1;
     end

   


    if f>10
        f=10;
    end
    j_max=f*j_max0;
    alpha_cs=f*alpha_cs0;


    % 第三步: 模型概率更新
    [u_IMM(:,t+1),L] = Model_P_up_SingerAP_CS_jerk(r_Singer, r_APCSJerk, S_Singer, S_APCSJerk, c_j);
    

    %      % 第四步： 转移概率TPM修正
    %定义阈值Th
    Th=0.9 ;
    
        index_match(t)=find(L==max(L), 1, 'first');
    
    if t>1
        i=1;
        for j=1:2
            if j~=index_match(t-1)
                Lrl(i)=L(index_match(t-1))/L(j);
                i=i+1;
            end
        end
        if min(Lrl)<=Th
            a=1;
            JC=JC+1;
        else
            a=0;
        end
        %实时更新Markov矩阵
        pij1=zeros(2,2);
        fij=zeros(4,4);
        %构造相同时刻的修正函数
        for i=1:2
            for j=1:2
                fij(i,j)=exp(0.1*(u_IMM(j,t+1)-u_IMM(i,t+1)));
            end
        end
        for j=1:2
            K1(j)=exp((1-a)*(u_IMM(j,t+1)-u_IMM(j,t))+a*(u_IMM_CKF(j,t+1)-u_IMM_CKF(j,t)));
            for i=1:2
                pij1(i,j)=K1(j)*pij(i,j); %修正
                %pij1(i,j)=fij(i,j)*K1(j)*pij(i,j); %修正 planB
            end
        end
        
        %归一化
        pij1sum=sum(pij1,2);
        for i=1:2
            for j=1:2
                pij(i,j)=pij1(i,j)/pij1sum(i);
            end
        end
        %设置转移矩阵的上下限，以保证主对角元占优
        T_th=0.95;
        T_tl=0.8;
        %考虑主对角元
        for i=1:2
            if pij(i,i)>T_th
                pij(i,i)=T_th;
                for j=1:2
                    if j~=i
                        pij(i,j)=(1-T_th)*pij1(i,j)/(pij1sum(i)-pij1(i,i));
                    end
                end
            elseif pij(i,i)<T_tl
                pij(i,i)=T_tl;
                for j=1:2
                    if j~=i
                        pij(i,j)=(1-T_tl)*pij1(i,j)/(pij1sum(i)-pij1(i,i));
                    end
                end
            end



         end



    end

    % 第五步: 估计融合
    [X_IMM(:,t+1), P_IMM(:,:,t+1)] = Model_mix_SingerAP_CS_jerk(X_Singer, X_APCSJerk, P_Singer, P_APCSJerk, u_IMM(:,t+1));

    % 保存各模型状态
    X_Singer_0(:, t+1) = X_Singer; 
    X_APCSJerk_0(:, t+1) = X_APCSJerk; 
    
end
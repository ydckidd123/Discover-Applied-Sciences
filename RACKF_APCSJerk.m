function [X, P, e_IMM, S_IMM, u_IMM, X_CKF,X_RCKF,X_ACKF,P_CKF,P_RCKF,P_ACKF,pij,D_RCKF_APCSJerk,D_ACKF_APCSJerk] = RACKF_APCSJerk(X_CKF,X_RCKF,X_ACKF,P_CKF,P_RCKF,P_ACKF, Z, A,  Q, R, u_IMM,pij,alpha_cs,j_max,t,D_RCKF_APCSJerk,D_ACKF_APCSJerk)
 n = length(X_ACKF); % 状态维度
    m = length(Z);         % 量测维度

    %记录初始概率
    u=u_IMM;
    % 第一步: 交互混合
    c_j = pij' * u_IMM; % 归一化常数
    % 计算混合概率
    ui1 = (1/c_j(1)) * pij(:,1) .* u_IMM;
    ui2 = (1/c_j(2)) * pij(:,2) .* u_IMM;
    ui3 = (1/c_j(3)) * pij(:,3) .* u_IMM;
    
    % 计算混合后的状态和协方差
    x01 = X_CKF*ui1(1) + X_ACKF*ui1(2) + X_RCKF*ui1(3);
    x02 = X_CKF*ui2(1) + X_ACKF*ui2(2) + X_RCKF*ui2(3);
    x03 = X_CKF*ui3(1) + X_ACKF*ui3(2) + X_RCKF*ui3(3);
    

    
    P01 = (P_CKF + (X_CKF - x01)*(X_CKF - x01)')*ui1(1) + ...
          (P_ACKF + (X_ACKF - x01)*(X_ACKF - x01)')*ui1(2) + ...
          (P_RCKF + (X_RCKF - x01)*(X_RCKF - x01)')*ui1(3);
    
    P02 = (P_CKF + (X_CKF - x02)*(X_CKF - x02)')*ui2(1) + ...
          (P_ACKF + (X_ACKF - x02)*(X_ACKF - x02)')*ui2(2) + ...
          (P_RCKF + (X_RCKF - x02)*(X_RCKF - x02)')*ui2(3);

    P03 = (P_CKF + (X_CKF - x03)*(X_CKF - x03)')*ui3(1) + ...
          (P_ACKF + (X_ACKF - x03)*(X_ACKF - x03)')*ui3(2) + ...
          (P_RCKF + (X_RCKF - x03)*(X_RCKF - x03)')*ui3(3);
    
    


    % 第二步: 模型条件滤波
    [X_CKF, P_CKF, r_CKF, S_CKF] = CKF_APCSJerk(x01, P01, Z, A,  Q, R,alpha_cs,j_max);
    [X_RCKF, P_RCKF, r_RCKF, S_RCKF,D_RCKF_APCSJerk] = RCKF_APCSJerk(x03, P03, Z, A,  Q, R,alpha_cs,j_max,t,D_RCKF_APCSJerk);
    [X_ACKF, P_ACKF, r_ACKF, S_ACKF,D_ACKF_APCSJerk] = ACKF_APCSJerk(x02, P02, Z, A,  Q, R,alpha_cs,j_max,t,D_ACKF_APCSJerk);
     % 第三步: 模型概率更新
    [u_IMM,L] = Model_P_up( r_CKF,r_RCKF, r_ACKF,S_CKF,  S_RCKF, S_ACKF, c_j);
    
    %设置转移矩阵的上下限，以保证主对角元占优
    T_th=0.7;
    T_tl=0.4;
    K=[0;0];
   
    %实时更新Markov矩阵
    pij1=zeros(2,2);
    for j=1:2
        K(j)=exp(u_IMM(j)-u(j));
        for i=1:2
        pij1(i,j)=K(j)*pij(i,j); %修正
        end
    end
     pij1sum=sum(pij1,2);

      for i=1:2
        for j=1:2
            pij(i,j)=pij1(i,j)/pij1sum(i);
        end
      end
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
     % 第四步: 估计融合
    [X, P, e_IMM, S_IMM] = Model_mix(X_CKF,X_RCKF,X_ACKF, P_CKF,P_RCKF,P_ACKF, r_CKF,r_RCKF,r_ACKF,S_CKF,S_RCKF,S_ACKF, u_IMM);

    %% 模型概率更新函数
function [u,L] = Model_P_up(r1,r2,r3,   S1,S2,S3,   c_j)
    % 计算似然函数（多元高斯分布）
    L1 = exp(-0.5 * r1' / S1 * r1) / sqrt(det(2*pi*S1));
    L2 = exp(-0.5 * r2' / S2 * r2) / sqrt(det(2*pi*S2));
    L3 = exp(-0.5 * r3' / S2 * r3) / sqrt(det(2*pi*S3));
    
    
    % 归一化似然
     L_sum = L1*c_j(1) + L2*c_j(2) + L3*c_j(3);
    % L_sum = L1 + L2 + L3+L4;

    L1 = L1 / L_sum;
    L2 = L2 / L_sum;
    L3 = L3 / L_sum;
    
    L=[L1; L2;L3];
    % 更新模型概率
    u = [L1; L2;L3] .* c_j;
    u = u / sum(u); % 归一化
end

%% 模型混合函数
function [x_pre, P, e, S] = Model_mix(x1,x2,x3, P1,P2,P3, r1,r2,r3, S1,S2,S3, u)
    % 综合状态估计
    x_pre = x1 * u(1) + x2 * u(2) + x3 * u(3) ;
    
    % 综合协方差估计
    P = u(1)*(P1 + (x1 - x_pre)*(x1 - x_pre)') + ...
        u(2)*(P2 + (x2 - x_pre)*(x2 - x_pre)')+ ...
        u(3)*(P3 + (x3 - x_pre)*(x3 - x_pre)');
   
    % 综合新息估计
    e = r1 * u(1) + r2 * u(2) + r3 * u(3) ;

    % 综合新息协方差估计
    S = u(1)*(S1 + (r1 - e)*(r1 - e)') + ...
        u(2)*(S2 + (r2 - e)*(r2 - e)')+ ...
        u(3)*(S3 + (r3 - e)*(r3 - e)');

% S = u(1)*S1 + u(2)*S2 + u(3)*S3 + ...
%           (u(1)*(r1*r1') + u(2)*(r2*r2') + u(3)*(r3*r3')) - e*e';

end

end
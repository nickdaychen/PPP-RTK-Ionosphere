function [ xk,pk,kof ] = kalman_filtering(sno,sls,pno,meas,x1k,p1k,sit,bp,freq,options)

c = 299792458;
%观测方程应为4n
mno = 4*sno;

Hk = zeros(mno,pno);
Zk = zeros(mno,  1);
Ck = zeros(mno,  1);
glok = [1 -4 5 6 1 -4 5 6 -2 -7 0 -1 -2 -7 0 -1 4 -3 3 2 4 -3 3 2 0 0];

Rk = eye(mno);
for k=1:sno
    %从model矩阵中获得卫星xyz坐标
    sat = meas((4*k)-3,8:10);
    
    rho = norm(sat - x1k(1:3,1)');
    
    Jx = (x1k(1,1) - sat(1,1))/rho;
    Jy = (x1k(2,1) - sat(1,2))/rho;
    Jz = (x1k(3,1) - sat(1,3))/rho;
    Jt = 1;
    Jw = meas((4*k - 3),28);
    Jtn= meas((4*k - 3),29);   
    Jte= meas((4*k - 3),30);
    %Jn = 1;
    Jr = 1;
    Je = 1;
    Jc = 1;
    %cdy
    Jk = 1;
    Jk1 = -1;
    %波长
    %v1=c/f1;
    %v2=c/f2;



%     s = (2*k - 1);
%     f = (2*k);
    %cdy
    s = (4*k - 3); 
    f = (4*k-2);
    m = (4*k-1); 
    g = 4*k;



%cdy
    Hk(s:g,1) = Jx;
    Hk(s:g,2) = Jy;
    Hk(s:g,3) = Jz;
    Hk(s:g,4) = Jt;
    Hk(s:g,5) = Jw;
%cdy
    Hk(s,(k+bp))= Jk;
    Hk(m,(k+bp))= Jk1;
%   Hk(f,(bp+k)) = Jn;
    switch options.TroGrad
        case 0
            if sit==3 || sit==4 || sit==5
                if sls(k)<33
                    Hk(s:g,6) = 0;
                else
                    Hk(s:g,6) = 1;
                end
            elseif sit==6 || sit==7
                if sls(k)<59
                    Hk(s:g,6) = 0;
                else
                    Hk(s:g,6) = 1;
                end
            elseif sit==8 || sit==9
                if sls(k)<33
                    Hk(s:g,6) = 0;
                    Hk(s:g,7) = 0;
                elseif sls(k)<59
                    Hk(s:g,6) = 1;
                    Hk(s:g,7) = 0;
                else
                    Hk(s:g,6) = 0;
                    Hk(s:g,7) = 1;
                end
            elseif sit==10
                if sls(k)<33
                    Hk(s:g,6) = 0;
                    Hk(s:g,7) = 0;
                elseif sls(k)<89
                    Hk(s:g,6) = 1;
                    Hk(s:g,7) = 0;
                else
                    Hk(s:g,6) = 0;
                    Hk(s:g,7) = 1;
                end
            elseif sit==11
                if sls(k)<59
                    Hk(s:g,6) = 0;
                    Hk(s:g,7) = 0;
                elseif sls(k)<89
                    Hk(s:g,6) = 1;
                    Hk(s:g,7) = 0;
                else
                    Hk(s:g,6) = 0;
                    Hk(s:g,7) = 1;
                end
            elseif sit==12
                if sls(k)<33
                    Hk(s:g,6) = 0;
                    Hk(s:g,7) = 0;
                    Hk(s:g,8) = 0;
                elseif sls(k)<59
                    Hk(s:g,6) = 1;
                    Hk(s:g,7) = 0;
                    Hk(s:g,8) = 0;
                elseif sls(k)<89
                    Hk(s:g,6) = 0;
                    Hk(s:g,7) = 1;
                    Hk(s:g,8) = 0;
                elseif sls(k)<106
                    Hk(s:g,6) = 0;
                    Hk(s:g,7) = 0;
                    Hk(s:g,8) = 1;
                end
            end
        case 1
            Hk(s:g,6) = Jtn;
            Hk(s:g,7) = Jte;
            if sit==3 || sit==4 || sit==5
                if sls(k)<33
                    Hk(s:g,8) = 0;
                else
                    Hk(s:g,8) = 1;
                end
            elseif sit==6 || sit==7
                if sls(k)<59
                    Hk(s:g,8) = 0;
                else
                    Hk(s:g,8) = 1;
                end
            elseif sit==8 || sit==9
                if sls(k)<33
                    Hk(s:g,8) = 0;
                    Hk(s:g,9) = 0;
                elseif sls(k)<59
                    Hk(s:g,8) = 1;
                    Hk(s:g,9) = 0;
                else
                    Hk(s:g,8) = 0;
                    Hk(s:g,9) = 1;
                end
            elseif sit==10
                if sls(k)<33
                    Hk(s:g,8) = 0;
                    Hk(s:g,9) = 0;
                elseif sls(k)<89
                    Hk(s:g,8) = 1;
                    Hk(s:g,9) = 0;
                else
                    Hk(s:g,8) = 0;
                    Hk(s:g,9) = 1;
                end
            elseif sit==11
                if sls(k)<59
                    Hk(s:g,8) = 0;
                    Hk(s:g,9) = 0;
                elseif sls(k)<89
                    Hk(s:g,8) = 1;
                    Hk(s:g,9) = 0;
                else
                    Hk(s:g,8) = 0;
                    Hk(s:g,9) = 1;
                end
            elseif sit==12
                if sls(k)<33
                    Hk(s:g,8) = 0;
                    Hk(s:g,9) = 0;
                    Hk(s:g,10) = 0;
                elseif sls(k)<59
                    Hk(s:g,8) = 1;
                    Hk(s:g,9) = 0;
                    Hk(s:g,10) = 0;
                elseif sls(k)<89
                    Hk(s:g,8) = 0;
                    Hk(s:g,9) = 1;
                    Hk(s:g,10) = 0;
                elseif sls(k)<106
                    Hk(s:g,8) = 0;
                    Hk(s:g,9) = 0;
                    Hk(s:g,10) = 1;
                end
            end
    end
    % fill the measurement vector
    % iono-free measurement
    % 应从obs矩阵中获取观测值，且以4个为一组填入观测方程Zk中
    %cdy
      Zk(s,1) = meas((4*k-3),6); 
      Zk(f,1) = meas((4*k-2),6);
      Zk(m,1) = meas((4*k-1),6);  
      Zk(g,1) = meas(4*k,6);

     if sls(k)<33
%         Zk(s,1) = i_free(meas((4*k - 3),6),meas((4*k - 2),6),0);
%         Zk(f,1) = i_free(meas((4*k - 1),6),meas((4*k    ),6),0);
          f1 = 10.23*10^6*154; 
          f2 = 10.23*10^6*120;    
     elseif sls(k)<59
%         Zk(s,1) = i_free(meas((4*k - 3),6),meas((4*k - 2),6),1);
%         Zk(f,1) = i_free(meas((4*k - 1),6),meas((4*k    ),6),1);
          f1 = (1602 + 0.5625*glok(sls(k)-32))*10^6; 
          f2 = (1246 + 0.4375*glok(sls(k)-32))*10^6;  
     elseif sls(k)<89
%         Zk(s,1) = i_free(meas((4*k - 3),6),meas((4*k - 2),6),2);
%         Zk(f,1) = i_free(meas((4*k - 1),6),meas((4*k    ),6),2);
          f1 = 10.23*10^6*154; 
          f2 = 10.23*10^6*115;  
     elseif sls(k)<106
%         Zk(s,1) = i_free(meas((4*k - 3),6),meas((4*k - 2),6),3);
%         Zk(f,1) = i_free(meas((4*k - 1),6),meas((4*k    ),6),3);
          f1 = 10.23*10^6*152.6; 
          f2 = 10.23*10^6*118;  
     end

    Jq = ((f1^2)/(f2^2));
    Jq1 = -((f1^2)/(f2^2));
    
    v1=c/f1;
    v2=c/f2;


    Hk(f,(k+bp))= Jq;
    Hk(g,(k+bp))= Jq1;
    Hk(m,(k+bp+sno)) = v1;
    Hk(g,(k+bp+2*sno)) = v2;



    % fill the computed vector
    % iono-free corrections for code and phase observations
    %此处应去掉无电离层函数，而是将误差校正量直接代入
    p1c = (meas((4*k - 3),7));
    p2c = (meas((4*k - 2),7));
    l1c = (meas((4*k - 1),7));
    l2c = (meas((4*k    ),7));
%     if sls(k)<33
%         pc = i_free(p1c,p2c,0);
%         lc = i_free(l1c,l2c,0);
%     elseif sls(k)<59
%         pc = i_free(p1c,p2c,1);
%         lc = i_free(l1c,l2c,1);
%     elseif sls(k)<89
%         pc = i_free(p1c,p2c,2);
%         lc = i_free(l1c,l2c,2);
%     elseif sls(k)<106
%         pc = i_free(p1c,p2c,3);
%         lc = i_free(l1c,l2c,3);
%     end
    %此处应替换为非差非组合模型，其余参数保持不变，而且为四个一组
    switch options.TroGrad
        case 0
            if sit==1
                Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jk*x1k(bp+k));
                Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jq*x1k(bp+k));
                Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jk1*x1k(bp+k))+(v1*x1k(bp+sno+k));
                Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jq1*x1k(bp+k))+(v2*x1k(bp+k+2*sno));

            elseif sit==2
%                   Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5));
%                   Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jn*x1k(bp+k));
                Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jk*x1k(bp+k));
                Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jq*x1k(bp+k));
                Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jk1*x1k(bp+k))+(v1*x1k(bp+sno+k));
                Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jq1*x1k(bp+k))+(v2*x1k(bp+k+2*sno));

            elseif sit==3 || sit==4 || sit==5
                if sls(k)<33
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5));
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jn*x1k(bp+k));

                      Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jk*x1k(bp+k));
                      Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jq*x1k(bp+k));
                      Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                      Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno)); 
                else
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6));
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jn*x1k(bp+k));

                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno)); 
                end
            elseif sit==6 || sit==7
                if sls(k)<59
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5));
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jn*x1k(bp+k));

                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno)); 
                else
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6));
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jn*x1k(bp+k));
                    
                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno)); 
                end
            elseif sit==8 || sit==9
                if sls(k)<33
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5));
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jn*x1k(bp+k));

                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno)); 
                elseif sls(k)<59
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6));
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jn*x1k(bp+k));

                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno)); 
                else
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(7));
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(7)) + (Jn*x1k(bp+k));

                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(7)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(7)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(7)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(7)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno)); 
                end
            elseif sit==10
                if sls(k)<33
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5));
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jn*x1k(bp+k));

                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno)); 
                elseif sls(k)<89
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6));
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jn*x1k(bp+k));

                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno)); 
                else
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(7));
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(7)) + (Jn*x1k(bp+k));

                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(7)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(7)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(7)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(7)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno)); 
                end
            elseif sit==11
                if sls(k)<59
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5));
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jn*x1k(bp+k));
                    
                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno)); 


                elseif sls(k)<89
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6));
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jn*x1k(bp+k));

                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno)); 
                else
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(7));
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(7)) + (Jn*x1k(bp+k));

                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(7)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(7)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(7)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(7)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno)); 
                end
            elseif sit == 12
                if sls(k)<33
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5));
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jn*x1k(bp+k));

                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno)); 
                elseif sls(k)<59
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6));
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jn*x1k(bp+k));

                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jr*x1k(6)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno)); 
                elseif sls(k)<89
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Je*x1k(7));
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Je*x1k(7)) + (Jn*x1k(bp+k));

                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Je*x1k(7)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Je*x1k(7)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Je*x1k(7)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Je*x1k(7)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno)); 
                elseif sls(k)<106
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jc*x1k(8));
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jc*x1k(8)) + (Jn*x1k(bp+k));

                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jc*x1k(8)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jc*x1k(8)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jc*x1k(8)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jc*x1k(8)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno)); 
                end
            end
        case 1
            if sit==1
%                 Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7));
%                 Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jn*x1k(bp+k));

                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno)); 

            elseif sit==2
%                 Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) ;
%                 Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jn*x1k(bp+k));

                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno)); 

            elseif sit==3 || sit==4 || sit==5
                if sls(k)<33
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) ;
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jn*x1k(bp+k));

                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno)); 

                else
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8));
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jn*x1k(bp+k));

                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno)); 

                end
            elseif sit==6 || sit==7
                if sls(k)<59
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) ;
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jn*x1k(bp+k));

                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno)); 
                else
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8));
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jn*x1k(bp+k));

                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno)); 

                end
            elseif sit==8 || sit==9
                if sls(k)<33
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) ;
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jn*x1k(bp+k));

                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno)); 
                elseif sls(k)<59
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8));
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jn*x1k(bp+k));

                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno));
                else
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(9));
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(9)) + (Jn*x1k(bp+k));

                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(9)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(9)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(9)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(9)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno));
                end
            elseif sit==10
                if sls(k)<33
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) ;
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jn*x1k(bp+k));

                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno)); 
                elseif sls(k)<89
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8));
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jn*x1k(bp+k));

                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno));
                else
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(9));
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(9)) + (Jn*x1k(bp+k));

                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(9)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(9)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(9)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(9)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno));
                end
            elseif sit==11
                if sls(k)<59
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) ;
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jn*x1k(bp+k));

                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno)); 
                elseif sls(k)<89
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8));
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jn*x1k(bp+k));

                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno));

                else
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(9));
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(9)) + (Jn*x1k(bp+k));

                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(9)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(9)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(9)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(9)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno));

                end
            elseif sit == 12
                if sls(k)<33
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) ;
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jn*x1k(bp+k));

                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno)); 
                elseif sls(k)<59
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8));
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jn*x1k(bp+k));

                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jr*x1k(8)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno));
                elseif sls(k)<89
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Je*x1k(9));
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Je*x1k(9)) + (Jn*x1k(bp+k));

                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Je*x1k(9)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Je*x1k(9)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Je*x1k(9)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Je*x1k(9)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno));
                elseif sls(k)<106
%                     Ck(s,1) = rho + pc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jc*x1k(10));
%                     Ck(f,1) = rho + lc + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jc*x1k(10)) + (Jn*x1k(bp+k));

                    Ck(s,1) = rho + p1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jc*x1k(10)) + (Jk*x1k(bp+k));
                    Ck(f,1) = rho + p2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jc*x1k(10)) + (Jq*x1k(bp+k));
                    Ck(m,1) = rho + l1c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jc*x1k(10)) + (Jk1*x1k(bp+k)) + (v1*x1k(bp+sno+k));
                    Ck(g,1) = rho + l2c + (Jt*x1k(4)) + (Jw*x1k(5)) + (Jtn*x1k(6)) + (Jte*x1k(7)) + (Jc*x1k(10)) + (Jq1*x1k(bp+k)) + (v2*x1k(bp+k+2*sno));
                end
            end
    end
    %基于高度角的随机模型确定测量噪声协方差矩阵Rk（4n*4n），且Rk为对角矩阵，随机模型或许可以得到优化
    if strcmp(options.WeMethod,'Elevation Dependent')
        
        if sls(k)<33
            f1 = freq(sls(k),1); 
            f2 = freq(sls(k),2);
            
%             ab1 = ((f1^2)/(f1^2 - f2^2))^2;
%             ab2 = ((f2^2)/(f1^2 - f2^2))^2;
            elv = meas((4*k - 1),26);
            
            kp = 1; kl = 1;
%             code_var = ((options.CodeStd *kp)^2)*(ab1+ab2);
%             phas_var = ((options.PhaseStd*kl)^2)*(ab1+ab2);
            
            code_var = (options.CodeStd *kp)^2;
            phas_var = (options.PhaseStd*kl)^2;
            
            Rk(s,s) = code_var/(sind(elv));
            Rk(f,f) = code_var/(sind(elv));            
            Rk(m,m) = phas_var/(sind(elv));
            Rk(g,g) = phas_var/(sind(elv));
            
        elseif sls(k)<59
            f1 = freq(sls(k),1); 
            f2 = freq(sls(k),2); 
            
%             ab1 = ((f1^2)/(f1^2 - f2^2))^2;
%             ab2 = ((f2^2)/(f1^2 - f2^2))^2;
            elv = meas((4*k - 1),26);
            
            if sit==2 || sit==6 || sit==7 || sit==11
                kp = 1; kl = 1;
            else
                kp = 2; kl = 1;
            end
%             code_var = ((options.CodeStd*kp)^2)*(ab1+ab2);
%             phas_var = ((options.PhaseStd*kl)^2)*(ab1+ab2);
            
%             Rk(s,s) = code_var/(sind(elv));
%             Rk(f,f) = phas_var/(sind(elv));

            code_var = (options.CodeStd*kp)^2;
            phas_var = (options.PhaseStd*kl)^2;
            
            Rk(s,s) = code_var/(sind(elv));
            Rk(f,f) = code_var/(sind(elv));            
            Rk(m,m) = phas_var/(sind(elv));
            Rk(g,g) = phas_var/(sind(elv));
            
        elseif sls(k)<89
            f1 = freq(sls(k),1); 
            f2 = freq(sls(k),2); 
            
%              ab1 = ((f1^2)/(f1^2 - f2^2))^2;
%              ab2 = ((f2^2)/(f1^2 - f2^2))^2;
            elv = meas((4*k - 1),26);
            
            kp = 2; kl = 2;
%             code_var = ((options.CodeStd *kp)^2)*(ab1+ab2);
%             phas_var = ((options.PhaseStd*kl)^2)*(ab1+ab2);
%             Rk(s,s) = code_var/(sind(elv));
%             Rk(f,f) = phas_var/(sind(elv));
            code_var = (options.CodeStd *kp)^2;
            phas_var = (options.PhaseStd*kl)^2;
            
            Rk(s,s) = code_var/(sind(elv));
            Rk(f,f) = code_var/(sind(elv));            
            Rk(m,m) = phas_var/(sind(elv));
            Rk(g,g) = phas_var/(sind(elv));
            
        elseif sls(k)<106
            f1 = freq(sls(k),1); 
            f2 = freq(sls(k),2);
            
%             ab1 = ((f1^2)/(f1^2 - f2^2))^2;
%             ab2 = ((f2^2)/(f1^2 - f2^2))^2;
            elv = meas((4*k - 1),26);
            
            kp = 2; kl = 2;
%             code_var = ((options.CodeStd *kp)^2)*(ab1+ab2);
%             phas_var = ((options.PhaseStd*kl)^2)*(ab1+ab2);
%             Rk(s,s) = code_var/(sind(elv));
%             Rk(f,f) = phas_var/(sind(elv));
            code_var = (options.CodeStd *kp)^2;
            phas_var = (options.PhaseStd*kl)^2;
            
            Rk(s,s) = code_var/(sind(elv));
            Rk(f,f) = code_var/(sind(elv));            
            Rk(m,m) = phas_var/(sind(elv));
            Rk(g,g) = phas_var/(sind(elv));
        end
    end
end
% 该代码使用该函数检查变量 Zk 和 Ck 中是否存在任何 NaN（非数字）值isnan。
% 如果有任何NaN 值，代码会找到它们出现的行索引，并从 Zk、Ck、Hk、Rk 和 Rk 的相应列中删除这些行。
if any(isnan(Zk)) || any(isnan(Ck))
    [a,~] = find(isnan(Zk));
    for b=a
        Zk(b,:) = [];
        Ck(b,:) = [];
        Hk(b,:) = [];
        Rk(b,:) = [];
        Rk(:,b) = [];
    end
end


sres0 = zeros(size(Zk,1),1);
while 1
    %Vk为观测噪声   Rk为其协方差
    Vk   = Zk - Ck; 
    %预测状态值的方差，即先验估计方差Sk
    Sk = (Hk*p1k*Hk') + Rk;
    %自适应卡尔曼滤波因子af计算
    abf = sum(Vk.^2)/trace(Sk);
    c0 = 2.5; c1 = 6.5;
    if abf>c1
        af = 10^10;
    elseif abf>c0
        af = (c0/abs(abf))*((c1 - abs(abf))/(c1 - c0))^2;
    else
        af = 1;
    end
    Sk = ((1/af).*(Hk*p1k*Hk')) + Rk;  
    %计算Kk卡尔曼增益
    %计算后验方差
    %计算当前历元基础估计值的方差
    % rcond 1范数条件数倒数估计值
    % pinv：摩尔彭若思伪逆矩阵
    if rcond(Sk)>1*10^-15         
        Kk = ((1/af).*(p1k*Hk'))/(Sk);
    else
        Kk = ((1/af).*(p1k*Hk'))*pinv(Sk);
    end
    dx   = Kk*Vk;
    %状态更新
    xk = x1k + dx;                
    tnk = (eye(pno) - Kk*Hk);
    %更新滤波估值协方差矩阵，即卡尔曼增益计算后的后验方差
    pk = tnk*p1k*tnk' + Kk*Rk*Kk';
    
    if rcond(Rk)>1*10^-15
        kof = pinv(Hk'*(Rk\Hk));      
    else
        kof = pinv(Hk'*(pinv(Rk)*Hk));      
    end
    kof = kof(1:5,1:5);
    res  = (Ck + (Hk*dx)) - Zk;   
    vres = abs(Rk - (Hk*pk*Hk'));      
    sres = zeros(size(res,1),1);
    for si = 1:size(res,1)
        sres(si,1) = abs(res(si,1))/sqrt(vres(si,si));
    end
    % 然后代码检查 sres 和 sres0 之间的差异是否大于 0.1。
    % 如果是，则代码更新测量噪声协方差矩阵Rk 以减少异常值对估计的影响
    dres = abs(sres - sres0);
    if any(dres>0.1)
        mm = find(abs(sres) == max(abs(sres)));
        k0 = 2.5; k1 = 6.5;
        if sres(mm,1)>k1
            sm = 1*10^-10;
            Rk(mm,mm) = Rk(mm,mm)/sm;
        elseif sres(mm,1)>k0
            sm = (k0/abs(sres(mm,1)))*((k1 - abs(sres(mm,1)))/(k1 - k0))^2;
            Rk(mm,mm) = Rk(mm,mm)/sm;
        end
        sres0 = sres;
    else
        break
    end
end
end

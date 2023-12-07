function [xs,kofs,pks] = MGNSS_filter(model,data,options)


[satlist] = dtr_satlist(data.obs);

[satno] = dtr_satno(data.obs);

%历元个数
rn = size(data.obs.st,1);

[freq,~] = frequencies;

[bp,~,sit] = dtr_sys(options);

xs  = zeros(bp+315,rn);
pks = zeros(bp+315,bp+315,rn);
kofs = zeros(5,5,rn);

n = 1;
% glok = [1 -4 5 6 1 -4 5 6 -2 -7 0 -1 -2 -7 0 -1 4 -3 3 2 4 -3 3 2 0 0];
  


for i=1:rn
    
    sls = satlist{i};
     % 可用卫星的编号（如2,5,7,85,97）
    sno = satno(i);
   % 当前历元可用卫星个数
    pno = bp + 3*sno;
    % 真正待估参数，需要修改3倍
    nk = n + (4*sno) - 1;
    meas = model(n:nk,:);
    
    n  = nk + 1;
    
    Q  = zeros(pno);
    %加入电离层延迟
    if options.ProMod == 0
        Q(1,1) = (options.NosPos*10^(options.NosPos2))*(data.inf.time.int);
        % int： 采样间隔
        Q(2,2) = (options.NosPos*10^(options.NosPos2))*(data.inf.time.int);
        Q(3,3) = (options.NosPos*10^(options.NosPos2))*(data.inf.time.int);
    end
    Q(4,4) = (options.NosClk*10^(options.NosClk2))*(data.inf.time.int);
    Q(5,5) = (options.NosTrop*10^(options.NosTrop2))*(data.inf.time.int);
    
    switch options.TroGrad
        case 0
            if bp==6
                Q(6,6) = (options.NosSTD*10^(options.NosSTD2))*(data.inf.time.int);
            elseif bp==7
                Q(6,6) = (options.NosSTD*10^(options.NosSTD2))*(data.inf.time.int);
                Q(7,7) = (options.NosSTD*10^(options.NosSTD2))*(data.inf.time.int);
            elseif bp==8
                Q(6,6) = (options.NosSTD*10^(options.NosSTD2))*(data.inf.time.int);
                Q(7,7) = (options.NosSTD*10^(options.NosSTD2))*(data.inf.time.int);
                Q(8,8) = (options.NosSTD*10^(options.NosSTD2))*(data.inf.time.int);
            end
        case 1
            Q(6,6) = (1*10^(-12))*(data.inf.time.int);
            Q(7,7) = (1*10^(-12))*(data.inf.time.int);
            if bp==8
                Q(8,8) = (options.NosSTD*10^(options.NosSTD2))*(data.inf.time.int);
            elseif bp==9
                Q(8,8) = (options.NosSTD*10^(options.NosSTD2))*(data.inf.time.int);
                Q(9,9) = (options.NosSTD*10^(options.NosSTD2))*(data.inf.time.int);
            elseif bp==10
                Q(8,8) = (options.NosSTD*10^(options.NosSTD2))*(data.inf.time.int);
                Q(9,9) = (options.NosSTD*10^(options.NosSTD2))*(data.inf.time.int);
                Q(10,10) = (options.NosSTD*10^(options.NosSTD2))*(data.inf.time.int);
            end
    end

    for ls = (bp+1) : (bp+sno)
        Q(ls,ls) = 10000;
    end

    for ls1 = (bp+sno+1) : pno
        Q(ls1,ls1) = 1000000000;
    end
    
    F = eye(pno);
    %    修改，增加电离层方差100
    if i == 1 && options.InMethod == 0
        
        x1k = zeros(pno,1);

        %添加电离层初值   
%           for ori = (bp+1) : (bp+sno)
%               if x1k(ori,1) == 0
%               x1k(ori,1) = 0.1;
%               end
%           end
        
        p1k = zeros(pno);
        p1k(1,1) = (options.IntPos*10^(options.IntPos2))^2;
        p1k(2,2) = (options.IntPos*10^(options.IntPos2))^2;
        p1k(3,3) = (options.IntPos*10^(options.IntPos2))^2;
        p1k(4,4) = (options.IntClk*10^(options.IntClk2))^2;
        p1k(5,5) = (options.IntTrop*10^(options.IntTrop2))^2;
        % int 系数
        switch options.TroGrad
            case 0
                if bp==6
                    p1k(6,6) = (options.IntSTD*10^(options.IntSTD2))^2;
                elseif bp==7
                    p1k(6,6) = (options.IntSTD*10^(options.IntSTD2))^2;
                    p1k(7,7) = (options.IntSTD*10^(options.IntSTD2))^2;
                elseif bp==8
                    p1k(6,6) = (options.IntSTD*10^(options.IntSTD2))^2;
                    p1k(7,7) = (options.IntSTD*10^(options.IntSTD2))^2;
                    p1k(8,8) = (options.IntSTD*10^(options.IntSTD2))^2;
                end
            case 1
                p1k(6,6) = (0*10^(1))^2;
                p1k(7,7) = (0*10^(1))^2;
                if bp==8
                    p1k(8,8) = (options.IntSTD*10^(options.IntSTD2))^2;
                elseif bp==9
                    p1k(8,8) = (options.IntSTD*10^(options.IntSTD2))^2;
                    p1k(9,9) = (options.IntSTD*10^(options.IntSTD2))^2;
                elseif bp==10
                    p1k(8,8) = (options.IntSTD*10^(options.IntSTD2))^2;
                    p1k(9,9) = (options.IntSTD*10^(options.IntSTD2))^2;
                    p1k(10,10) = (options.IntSTD*10^(options.IntSTD2))^2;
                end
        end
        %amb：整周模糊度系数

        %cdy

%         for u=(bp+1):pno
%             p1k(u,u) = (options.IntAmb*10^(options.IntAmb2))^2;
%         end
        for v = (bp+1):(bp+sno)
            p1k(v,v) = 1000000;
        
        end

        for u=(bp+sno+1):pno
            p1k(u,u) = (options.IntAmb*10^(options.IntAmb2))^2;
        end

    elseif i ~= 1
        %将上一个历元的基础参数和其他参数存入x1k中
        %修改数据存储逻辑后，结果有提升
        x1k = zeros(pno,1);
        x1k(1:bp,1) = xs(1:bp,i-1);
        %cdy
        for k=1:3*sno
            if k < (1+sno)
            snm = sls(k);
            x1k(bp+k,1) = xs(bp+snm,i-1);
            elseif k < (2*sno+1)
            snm = sls(k-sno) + 105;
            x1k(bp+k,1) = xs(bp+snm,i-1);
            elseif k < (3*sno+1)
            snm = sls(k-2*sno) + 210;
            x1k(bp+k,1) = xs(bp+snm,i-1);
            end
        end
        %F为状态转移矩阵
        x1k = F*x1k;
        
        p1k = zeros(pno);
        %获取p1k第一二维度的长度
        for r=1:size(p1k,1)
            for c=1:size(p1k,2)

           if r<(bp+1) && c<(bp+1)
               p1k(r,c) = ps(r,c);
           elseif r <(bp+1) && c > bp && c < (bp+sno+1)
               f2 = sls(c-bp);
               p1k(r,c) = ps(r,f2+bp);
           elseif r < (bp+1) && c > (bp+sno) && c < (bp+2*sno+1)
               f2 = sls(c-sno-bp);
               p1k(r,c) = ps(r,f2+bp+sno);
           elseif r < (bp+1) && c > (bp+2*sno) && c < (bp+3*sno+1)
               f2 = sls(c-2*sno-bp);
               p1k(r,c) = ps(r,f2+bp+2*sno);


           elseif r > bp &&  r < (bp+sno+1) && c<(bp+1)
                f1 = sls(r-bp);
                p1k(r,c) = ps(f1+bp,c);
           elseif r > bp &&  r < (bp+sno+1) && c > bp && c < (bp+sno+1)
               f1 = sls(r-bp);
               f2 = sls(c-bp);
               p1k(r,c) = ps(f1+bp,f2+bp);
           elseif r > bp &&  r < (bp+sno+1) && c > (bp+sno) && c < (bp+2*sno+1)
               f1 = sls(r-bp);
               f2 = sls(c-sno-bp);
               p1k(r,c) = ps(f1+bp,f2+bp+sno);
           elseif r > bp &&  r < (bp+sno+1) && c > (bp+2*sno) && c < (bp+3*sno+1)
               f1 = sls(r-bp);
               f2 = sls(c-2*sno-bp);
               p1k(r,c) = ps(f1+bp,f2+bp+2*sno);


           elseif r > (bp+sno) &&  r < (bp+2*sno+1) && c<(bp+1)
               f1 = sls(r-bp-sno);
                p1k(r,c) = ps(f1+bp+sno,c);
           elseif r > (bp+sno) &&  r < (bp+2*sno+1) && c > bp && c < (bp+sno+1)
               f1 = sls(r-bp-sno);
               f2 = sls(c-bp);
               p1k(r,c) = ps(f1+bp+sno,f2+bp);
           elseif r > (bp+sno) &&  r < (bp+2*sno+1) && c > (bp+sno) && c < (bp+2*sno+1)
               f1 = sls(r-bp-sno);
               f2 = sls(c-sno-bp);
               p1k(r,c) = ps(f1+bp+sno,f2+bp+sno);
           elseif r > (bp+sno) &&  r < (bp+2*sno+1) && c > (bp+2*sno) && c < (bp+3*sno+1)
               f1 = sls(r-bp-sno);
               f2 = sls(c-2*sno-bp);
               p1k(r,c) = ps(f1+bp+sno,f2+bp+2*sno);
           
           elseif r > (bp+2*sno) &&  r < (bp+3*sno+1) && c<(bp+1)
               f1 = sls(r-bp-2*sno);
               p1k(r,c) = ps(f1+bp+2*sno,c);
           elseif r > (bp+2*sno) &&  r < (bp+3*sno+1) && c > bp && c < (bp+sno+1)
               f1 = sls(r-bp-2*sno);
               f2 = sls(c-bp);
               p1k(r,c) = ps(f1+bp+2*sno,f2+bp);
           elseif r > (bp+2*sno) &&  r < (bp+3*sno+1) && c > (bp+sno) && c < (bp+2*sno+1)
               f1 = sls(r-bp-2*sno);
               f2 = sls(c-sno-bp);
               p1k(r,c) = ps(f1+bp+2*sno,f2+bp+sno);
           else
               f1 = sls(r-bp-2*sno);
               f2 = sls(c-2*sno-bp);
               p1k(r,c) = ps(f1+bp+2*sno,f2+bp+2*sno) ;
%                 if r<(bp+1) && c<(bp+1)
%                     p1k(r,c) = ps(r,c);
%                 elseif r<(bp+1) && c>bp
%                     sn = sls(c-bp);
%                     p1k(r,c) = ps(r,sn+bp);
%                 elseif r>bp && c<(bp+1)
%                     sn = sls(r-bp);
%                     p1k(r,c) = ps(sn+bp,c);
%                 else
%                     f1 = sls(r-bp);
%                     f2 = sls(c-bp);
%                     p1k(r,c) = ps(f1+bp,f2+bp);
%                 end
           end

            end
        end
        
        for k=1:3*sno
            if k < (sno+1)
                snm = sls(k)+bp;
                if ps(snm,snm)==0
                p1k(k+bp,k+bp) = 1000000;
                end
            elseif k < (2*sno+1)
                snm = sls(k-sno) + bp + 105;
                if ps(snm,snm)==0
                p1k(k+bp,k+bp) = (options.IntAmb*10^(options.IntAmb2))^2;
                end
            else
                snm = sls(k-2*sno) + bp + 210;
                if ps(snm,snm)==0
                p1k(k+bp,k+bp) = (options.IntAmb*10^(options.IntAmb2))^2;
                end
            end
        end
        
        p1k = F*p1k*F' + Q;
    end
    
    if i == 1
        if options.InMethod == 1
            [ xk,pk,kof ] = least_sqr(sno,sls,pno,meas,sit,bp,freq,options);
        else
            % initial point
            if strcmp(options.ApMethod,'RINEX') % first choice is from RINEX
                x1k(1:3,1) = data.inf.rec.pos';
            elseif strcmp(options.ApMethod,'Specify')% second choice is from User specific
                x1k(1:3,1) = [options.AprioriX options.AprioriY options.AprioriZ];
            end
            [ xk,pk,kof ] = kalman_filtering(sno,sls,pno,meas,x1k,p1k,sit,bp,freq,options);
        end
    else
        [ xk,pk,kof ] = kalman_filtering(sno,sls,pno,meas,x1k,p1k,sit,bp,freq,options);
    end 
    
    kofs(:,:,i) = kof(:,:,1);
    
    xs(1:bp,i) = xk(1:bp,1);

    %cdy
%     for k = 1:sno
%        snm = sls(k) + bp;
%        xs(snm,i) = xk(bp+k,1);
%     end
    for k = 1:3*sno
        if k < sno + 1
            snm = sls(k) + bp;
%             if sls(k) < 33
%             f1 = 10.23*10^6*154; 
%             f2 = 10.23*10^6*120;
%             fdcb = ((f2^2)/(f1^2-f2^2));
%             xs(snm,i) = xk(bp+k,1) + fdcb*(dcb(sls(k),1));
%             elseif sls(k) < 59
%             f1 = (1602 + 0.5625*glok(sls(k)-32))*10^6; 
%             f2 = (1246 + 0.4375*glok(sls(k)-32))*10^6;
%             fdcb = ((f2^2)/(f1^2-f2^2));
%             xs(snm,i) = xk(bp+k,1) + fdcb*(dcb(sls(k),1));
%             else
            xs(snm,i) = xk(bp+k,1);
%             end
        elseif k < 2*sno + 1
            snm = sls(k-sno) + bp + 105;
            xs(snm,i) = xk(bp+k,1);
        else
            snm =  sls(k-2*sno) + bp + 210;
            xs(snm,i) = xk(bp+k,1);
        end
    end


    %cdy
    ps = zeros(bp + 315);
    for r=1:size(pk,1)
       for c=1:size(pk,2)
%            if r<(bp+1) && c<(bp+1)
%                ps(r,c) = pk(r,c);
%            elseif r<(bp+1) && c>bp
%                sn = sls(c-bp);
%                ps(r,sn+bp) = pk(r,c);
%            elseif r>bp && c<(bp+1)
%                sn = sls(r-bp);
%                ps(sn+bp,c) = pk(r,c);
%            else
%                sn1 = sls(r-bp); sn2 = sls(c-bp);
%                ps(sn1+bp,sn2+bp) = pk(r,c);
%            end

           if r<(bp+1) && c<(bp+1)
                ps(r,c) = pk(r,c);
           elseif r <(bp+1) && c > bp && c < (bp+sno+1)
               sn = sls(c-bp);
               ps(r,sn+bp) = pk(r,c);
           elseif r < (bp+1) && c > (bp+sno) && c < (bp+2*sno+1)
               sn = sls(c-sno-bp);
               ps(r,sn+bp+sno) = pk(r,c);
           elseif r < (bp+1) && c > (bp+2*sno) && c < (bp+3*sno+1)
               sn = sls(c-2*sno-bp);
               ps(r,sn+bp+2*sno) = pk(r,c);


           elseif r > bp &&  r < (bp+sno+1) && c<(bp+1)
               sn = sls(r-bp);
               ps(sn+bp,c) = pk(r,c);
           elseif r > bp &&  r < (bp+sno+1) && c > bp && c < (bp+sno+1)
               sn1 = sls(r-bp);
               sn2 = sls(c-bp);
               ps(sn1+bp,sn2+bp) = pk(r,c);
           elseif r > bp &&  r < (bp+sno+1) && c > (bp+sno) && c < (bp+2*sno+1)
               sn1 = sls(r-bp);
               sn2 = sls(c-sno-bp);
               ps(sn1+bp,sn2+bp+sno) = pk(r,c);
           elseif r > bp &&  r < (bp+sno+1) && c > (bp+2*sno) && c < (bp+3*sno+1)
               sn1 = sls(r-bp);
               sn2 = sls(c-2*sno-bp);
               ps(sn1+bp,sn2+bp+2*sno) = pk(r,c);


           elseif r > (bp+sno) &&  r < (bp+2*sno+1) && c<(bp+1)
               sn = sls(r-bp-sno);
               ps(sn+bp+sno,c) = pk(r,c);
           elseif r > (bp+sno) &&  r < (bp+2*sno+1) && c > bp && c < (bp+sno+1)
               sn1 = sls(r-bp-sno);
               sn2 = sls(c-bp);
               ps(sn1+bp+sno,sn2+bp) = pk(r,c);
           elseif r > (bp+sno) &&  r < (bp+2*sno+1) && c > (bp+sno) && c < (bp+2*sno+1)
               sn1 = sls(r-bp-sno);
               sn2 = sls(c-sno-bp);
               ps(sn1+bp+sno,sn2+bp+sno) = pk(r,c);
           elseif r > (bp+sno) &&  r < (bp+2*sno+1) && c > (bp+2*sno) && c < (bp+3*sno+1)
               sn1 = sls(r-bp-sno);
               sn2 = sls(c-2*sno-bp);
               ps(sn1+bp+sno,sn2+bp+2*sno) = pk(r,c);
           
           elseif r > (bp+2*sno) &&  r < (bp+3*sno+1) && c<(bp+1)
               sn = sls(r-bp-2*sno);
               ps(sn+bp+2*sno,c) = pk(r,c);
           elseif r > (bp+2*sno) &&  r < (bp+3*sno+1) && c > bp && c < (bp+sno+1)
               sn1 = sls(r-bp-2*sno);
               sn2 = sls(c-bp);
               ps(sn1+bp+2*sno,sn2+bp) = pk(r,c);
           elseif r > (bp+2*sno) &&  r < (bp+3*sno+1) && c > (bp+sno) && c < (bp+2*sno+1)
               sn1 = sls(r-bp-2*sno);
               sn2 = sls(c-sno-bp);
               ps(sn1+bp+2*sno,sn2+bp+sno) = pk(r,c);
           else
               sn1 = sls(r-bp-2*sno);
               sn2 = sls(c-2*sno-bp);
               ps(sn1+bp+2*sno,sn2+bp+2*sno) = pk(r,c);
           end
       end
    end
    pks(:,:,i) = ps;
end
xs(1:3,:) = xs(1:3,:);
end
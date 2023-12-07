clc,clear
tic
load('D:\data\testData\results\zwd_true_value\zwd.mat');
for week=2223:2242
    % 定义一个包含各个站点名的cell数组
    stations = {'hkcl', 'hkfn', 'hkks', 'hkkt', 'hklm', 'hklt', 'hkmw', 'hknp', 'hkoh', 'hkpc', 'hkqt', 'hksc', 'hksl', 'hkss', 'hkst', 'hktk', 'hkws', 'kyc1', 't430'};
    % 选择当前站点
    station_name = stations{7};
    % 调用loadfiles函数加载指定周数和站点的相关文件，返回RINEX观测文件、轨道文件、天线文件和钟差文件
    [rinex,orbitb,anten,clock] = loadfiles(week, station_name);
    % 定义每个RINEX文件的数据点数量
    rn = 2880;
    % 定义需要处理的RINEX文件数量
    num_files = 8;
    % 调用frequencies函数获取频率列表
    [freq,~] = frequencies;
    % 设置滤波器的带宽
    bp = 6;
    % 定义站点编号（或站点类型）
    sit = 3;
    % 初始化存储每个站点处理结果的数组，包括额外的105行作为缓冲
    xs = zeros(bp+105, rn*num_files);
    % 初始化相关联的协方差矩阵
    pks = zeros(bp+105, bp+105, rn*num_files);
    % 初始化每个站点偏移的协方差矩阵
    kofs = zeros(5, 5, rn*num_files);
    % 初始化索引n，通常用于后续循环中 
    n = 1;
    % 初始化存储所有卫星列表的cell数组
    total_satlist = cell(rn, num_files);
    % 初始化存储每个时间点的卫星编号数量的矩阵
    total_satno = double(zeros(rn, num_files));
    % 初始化存储每个文件模型的cell数组
    total_model = cell(1, num_files);
   % 初始化一个全为1的列向量，用于追踪每个文件在处理中的位置（可能用于定位或索引）
    key = ones(num_files,1);
    % 循环处理每个文件
    for fil=1:num_files
        % 调用calculate_zwd函数处理数据，获取卫星列表、卫星编号、模型和位置信息
        fil
        [satlist, satno, model, x_position] = calculate_zwd(rinex(fil,:), orbitb(fil,:), anten, clock(fil,:));
        % 将每个文件的卫星列表存入对应的列
        total_satlist(:,fil) = satlist;
        % 将每个文件的卫星编号数量存入对应的列
        total_satno(:,fil) = satno;
        % 获取模型数组的长度
        [long,~] = size(model);
        % 存储每个文件的模型数据
        total_model{fil} = model;
        % 如果当前文件是第31个文件，则提前终止循环
        if fil == 31
            break;
        else
            % 否则，更新下一个文件的key值，这个值可能是为了跟踪处理的进度或是用于后续计算的索引
            key(fil+1) = key(fil) + long;
        end

    end
    [options]=get_options;
    for i=1:rn*num_files

        day = fix(i / rn) + 1;   % Calculate the day index based on the iteration count i.
        set_i=mod(i,rn);
        if set_i==0&&i~=0
            day=day-1;
            set_i=rn;
        elseif set_i==1
            n=1;
        end
        satlist = total_satlist(:,day);
        satno = total_satno(:,day);
        sls = satlist{set_i};  % Get satellite list at time i.
        sno = satno(set_i);    % Get the number of satellites at time i.
        pno = bp + sno;          % Calculate the total number of parameters (base parameters + satellites).
        nk = n + (4*sno) - 1;
        model=total_model{day};
        meas = model(n:nk,:);
        n  = nk + 1;   
        Q  = zeros(pno);
        Q(4,4) = (1*10^(5))*(30);
        Q(5,5) = (1*10^(-9))*(30);
        Q(6,6) = (1*10^(-7))*(30);
        F = eye(pno);   
        if i == 1 
            %初始化状态向量和协方差矩阵
            x1k = zeros(pno,1);
            p1k = zeros(pno);
            %为协方差矩阵设置初始值
            p1k(1,1) = (1*10^(2))^2;
            p1k(2,2) = (1*10^(2))^2;
            p1k(3,3) = (1*10^(2))^2;
            p1k(4,4) = (1*10^(5))^2;
            p1k(5,5) = (0.5*10^(0))^2;
            p1k(6,6) = (1*10^(2))^2;    
            %为每个卫星参数设置初始协方差
            for u=(bp+1):pno
                p1k(u,u) = (2*10^(1))^2;
            end
        elseif i ~= 1
            x1k = zeros(pno,1);
            x1k(1:bp,1) = xs(1:bp,i-1);
            %根据卫星列表更新状态向量
            for k=1:sno
                snm = sls(k);
                x1k(bp+k,1) = xs(bp+snm,i-1);
            end
            %通过状态转移矩阵更新状态向量
            x1k = F*x1k;    
            %初始化协方差矩阵，并根据先前协方差和卫星列表更新
            p1k = zeros(pno);
            for r=1:size(p1k,1)
                for c=1:size(p1k,2)
                    if r<(bp+1) && c<(bp+1)
                        p1k(r,c) = ps(r,c);
                    elseif r<(bp+1) && c>bp
                        sn = sls(c-bp);
                        p1k(r,c) = ps(r,sn+bp);
                    elseif r>bp && c<(bp+1)
                        sn = sls(r-bp);
                        p1k(r,c) = ps(sn+bp,c);
                    else
                        f1 = sls(r-bp);
                        f2 = sls(c-bp);
                        p1k(r,c) = ps(f1+bp,f2+bp);
                    end
                end
            end
            for k=1:sno
                snm = sls(k)+bp;
                if ps(snm,snm)==0
                    p1k(k+bp,k+bp) = (2*10^(1))^2;
                end
            end      
            p1k = F*p1k*F' + Q;
        end  
        if i == 1
          x1k(1:3,1) = x_position;
         [ xk,pk,kof,res ] = kalman_filtering(sno,sls,pno,meas,x1k,p1k,sit,bp,freq);
        else
         [ xk,pk,kof,res ] = kalman_filtering(sno,sls,pno,meas,x1k,p1k,sit,bp,freq);
        end 
        kofs(:,:,i) = kof(:,:,1);

        xs(1:bp,i) = xk(1:bp,1);
        for k = 1:sno
           snm = sls(k) + bp;
           xs(snm,i) = xk(bp+k,1);
        end
        ps = zeros(bp + 105);
        for r=1:size(pk,1)
           for c=1:size(pk,2)
               if r<(bp+1) && c<(bp+1)
                   ps(r,c) = pk(r,c);
               elseif r<(bp+1) && c>bp
                   sn = sls(c-bp);
                   ps(r,sn+bp) = pk(r,c);
               elseif r>bp && c<(bp+1)
                   sn = sls(r-bp);
                   ps(sn+bp,c) = pk(r,c);
               else
                   sn1 = sls(r-bp); sn2 = sls(c-bp);
                   ps(sn1+bp,sn2+bp) = pk(r,c);
               end
           end
        end
        pks(:,:,i) = ps;
        tzd=xs(5,:)';
    end
    start_number=2880*7*(week-2191)+2881;
    end_number=start_number+2880*7-1;
    zwd.(station_name)(start_number:end_number,1)=tzd(2881:end,1);
    save('D:\\data\\testData\\results\\zwd_true_value\\zwd','zwd');
    clearvars -except week zwd;
    clc;
end
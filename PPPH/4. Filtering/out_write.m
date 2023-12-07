function [] = out_write(xs,Option,model,filename1,filename2,data)


year = data.inf.time.first(1);
doy  = data.inf.time.doy;

dat = xs';

td = model(1,31);
tw = dat(:,5);
tt = td + tw;
% el = 0.25*size(model,1) ;
% angle = zeros(el,1);
% for i = 1 : el
%     satele = 4*i;
%     angle(i,1) = model(satele,26);
% end

name = {'Year','DOY','SOD','(X)','(Y)','(Z)','(DT)','(TH)','(TW)','(TT)','(G01)','(G02)','(G03)','(G04)','(G05)','(G06)','(G07)','(G08)','(G09)','(G10)','(G11)','(G12)','(G13)','(G14)','(G15)','(G16)','(G17)','(G18)','(G19)','(G20)','(G21)','(G22)','(G23)','(G24)','(G25)','(G26)','(G27)','(G28)','(G29)','(G30)','(G31)','(G32)'};
name1 = {'Year','DOY','SOD','(G01)','(G02)','(G03)','(G04)','(G05)','(G06)','(G07)','(G08)','(G09)','(G10)','(G11)','(G12)','(G13)','(G14)','(G15)','(G16)','(G17)','(G18)','(G19)','(G20)','(G21)','(G22)','(G23)','(G24)','(G25)','(G26)','(G27)','(G28)','(G29)','(G30)','(G31)','(G32)'};

% name = {'Year','DOY','SOD','(X)','(Y)','(Z)','(DT)','(TH)','(TW)','(TT)','(G01)','(G02)','(G03)'};
% for i = 1:32
%     G_num = sprintf('%02d', i);  % 将循环变量格式化为两位数的字符串，例如 '01'、'02'、'03'
%     name{end+1} = ['(G', G_num, ')'];  % 将 '(G01)'、'(G02)'、'(G03)' 等字符串添加到 name 列表末尾
% end

format1 = '%4s %3s %5s %13s %13s %13s %10s %7s %7s %7s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s';
format2 = '%4d %3d %5d %13.3f %13.3f %13.3f %10.3f %7.3f %7.3f %7.3f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f';
format3 = '%4s %3s %5s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s %13s';
format4 = '%4d %3d %5d %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f %13.5f';

% format1 = '%4s %3s %5s %13s %13s %13s %10s %7s %7s %7s %13s %13s %13s';
% format2 = '%4d %3d %5d %13.3f %13.3f %13.3f %10.3f %7.3f %7.3f %7.3f %13.5f %13.5f %13.5f';

if Option.TroGrad==1
    nn = size(name,2);
    name{nn+1} = '(TGN)';
    name{nn+2} = '(TGE)';
    
    format1 = [format1,' %7s %7s'];
    format2 = [format2,' %7.4f %7.4f'];
end

if Option.system.glo==1
    nn = size(name,2);
    name{nn+1} = '(SDR)';
    
    format1 = [format1,' %7s'];
    format2 = [format2,' %7.3f'];
end

if Option.system.gal==1
    nn = size(name,2);
    name{nn+1} = '(SDE)';
    
    format1 = [format1,' %7s'];
    format2 = [format2,' %7.3f'];
end

if Option.system.bds==1
    nn = size(name,2);
    name{nn+1} = '(SDC)';
    
    format1 = [format1,' %7s'];
    format2 = [format2,' %7.3f'];
end

fid1 = fopen(filename1,'wt');
fid2 = fopen(filename2,'wt');

fprintf(fid1,[format1,'\n'],name{:});
fprintf(fid2,[format3,'\n'],name1{:});

sm = Option.system.glo + Option.system.gal + Option.system.bds;

for i=1:size(dat,1)
    if Option.TroGrad==1
        if sm==0
            fprintf(fid1,[format2,'\n'],year,doy,data.obs.ep(i,1),...
                dat(i,1:4),td(1),tw(i,1),tt(i,1),dat(i,6:7),dat(i,8:39));
        elseif sm==1
            fprintf(fid1,[format2,'\n'],year,doy,data.obs.ep(i,1),...
                dat(i,1:4),td(1),tw(i,1),tt(i,1),dat(i,6:7),dat(i,9:40),dat(i,8));
        elseif sm==2
            fprintf(fid1,[format2,'\n'],year,doy,data.obs.ep(i,1),...
                dat(i,1:4),td(1),tw(i,1),tt(i,1),dat(i,6:7),dat(i,10:41),dat(i,8:9));
        elseif sm==3
            fprintf(fid1,[format2,'\n'],year,doy,data.obs.ep(i,1),...
                dat(i,1:4),td(1),tw(i,1),tt(i,1),dat(i,6:7),dat(i,11:42),dat(i,8:10));
        end     
    else
        if sm==0
            fprintf(fid1,[format2,'\n'],year,doy,data.obs.ep(i,1),...
                dat(i,1:4),td(1),tw(i,1),tt(i,1),dat(i,6:37));
        elseif sm==1
            fprintf(fid1,[format2,'\n'],year,doy,data.obs.ep(i,1),...
                dat(i,1:4),td(1),tw(i,1),tt(i,1),dat(i,7:38),dat(i,6));
        elseif sm==2
            fprintf(fid1,[format2,'\n'],year,doy,data.obs.ep(i,1),...
                dat(i,1:4),td(1),tw(i,1),tt(i,1),dat(i,8:39),dat(i,6:7));
        elseif sm==3
            fprintf(fid1,[format2,'\n'],year,doy,data.obs.ep(i,1),...
                dat(i,1:4),td(1),tw(i,1),tt(i,1),dat(i,9:40),dat(i,6:8));
        end
    end
    fprintf(fid2,[format4,'\n'],year,doy,data.obs.ep(i,1),data.obs.elv(i,1:32));

end
    fclose(fid1);
    fclose(fid2);
end


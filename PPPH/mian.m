 clc,clear

rinex01 = 'D:\data\testData\ObservationFiles\hkfn3650.22o';
rinex02 = 'D:\data\testData\ObservationFiles\hkfn0010.22o';
rinex03 = 'D:\data\testData\ObservationFiles\hkfn0020.22o';

clock01 = 'D:\data\testData\clc\igs21905.clk_30s';
clock02 = 'D:\data\testData\clc\igs21906.clk_30s';
clock03 = 'D:\data\testData\clc\igs21910.clk_30s';

orbitb01 = 'D:\data\testData\clc\igs21905.sp3';
orbitb02 = 'D:\data\testData\clc\igs21906.sp3';
orbitb03 = 'D:\data\testData\clc\igs21910.sp3';
anten='D:\data\atx\igs20.atx';

files.rinex = rinex01;
files.orbitb = '';
files.orbit = orbitb01;
files.orbita = '';
files.anten = anten;
files.clock = clock01;
files.dcb = '';

options.clock = 'Clock File';
options.system.gps = 1;
options.system.glo = 0;
options.system.gal = 0;
options.system.bds = 0;
options.dcb = 0;

options.clck_int = 30;
[data] = data_hand(files,options);
options.elvangle =8;
options.CSMw = 1;
options.CSGf = 1;
options.clkjump = 1;
options.codsmth = 0;
options.Static = 1;
if options.Static == 1
    options.ProMod = 1; %static
else
    options.ProMod = 0; %kinematic
end

% model
options.SatClk = 1;
options.SatAPC = 1;
options.SatWind = 1;
options.RecAPC = 1;
options.RecARP = 1;
options.AtmTrop = 1;
options.TroGrad = 0;
options.RelClk = 1;
options.RelPath = 1;
options.Solid = 1;
% filter
options.InMethod = 0;
if options.InMethod == 0
    options.IntPos = 1;
    options.IntPos2 = 2;
    options.IntClk = 1;
    options.IntClk2 = 5;
    options.IntTrop = 0.5;
    options.IntTrop2 = 0;
    options.IntSTD = 1;
    options.IntSTD2 = 2;
end
options.IntAmb = 2;
options.IntAmb2 = 1;
options.NosPos = 0;
options.NosPos2 = 0;
options.NosClk = 1;
options.NosClk2 = 5;
options.NosTrop = 1;
options.NosTrop2 = -9;
options.NosSTD = 1;
options.NosSTD2 = -7;
options.ApMethod = 'RINEX';
if strcmp(options.ApMethod,'Specify')
    options.AprioriX = 0; 
    options.AprioriY = 0;
    options.AprioriZ = 0;
end
options.WeMethod = 'Elevation Dependent';
options.CodeStd = 3;
options.PhaseStd = 0.003;
% general
options.from = 0;
options.to = 86370;
[data] = preprocess(data,options);
[model] = nmodel(data,options);
[xs,kofs,pks] = MGNSS_filter(model,data,options);
[ppps] = data_hand(files,options);
app.Data = ppps;
tzd=xs(5,:)';
% boundaries
st = ((0 - app.Data.obs.ep(1,1))/(app.Data.inf.time.int)) + 1;
fn = ((86370 - app.Data.obs.ep(1,1))/(app.Data.inf.time.int)) + 1;
t = app.Data.obs.ep(st:fn,1);
t = t./3600; % hour
% 
% fig = figure('Name','Tropospheric Zenith Total Delay','NumberTitle','off','Color',[0.75 0.75 0.75]);
% ax  = axes(fig);
% plot(t,tzd)
% ax.Title.String = 'Tropospheric Zenith Total Delay';
% ax.XLabel.String = 'Time (Hour)';
% ax.YLabel.String = 'Zenith Total Delay (Meter)';
% ax.XGrid = 'on';
% ax.YGrid = 'on';
% legend(ax,'Zenith Total Delay')
% HKLMref=[-2414045.945;5391602.352;2396878.889];南丫岛
% HKMWref=[-2402484.109;5395262.438;2400726.956];梅窝
% HKOHref=[-2423816.900;5386057.093;2399883.371];石碑山
% HKSLref=[-2393382.416;5393861.175;2412592.411];小冷水
% HKQTREF=[-2421567.892;5384910.563;2404264.394];zei鱼涌
% ref = [-2393382.416;5393861.175;2412592.411];
% [n,e,u,~,~,~] = evaluate(xs,ref);
% boundaries
% st = ((app.PFrom.Value - app.Data.obs.ep(1,1))/(app.Data.inf.time.int)) + 1;
% fn = ((app.PTo.Value - app.Data.obs.ep(1,1))/(app.Data.inf.time.int)) + 1;
% t = app.Data.obs.ep(st:fn,1);
% t = t./3600; % hour
% fig = figure('Name','Positioning Errors and ZWD','NumberTitle','off','Color',[0.75 0.75 0.75]);
% ax = axes(fig);
% plot(ax,t,n(:,1),t,e(:,1),t,u(:,1),t,tzd)
% ax.Title.String = '20210626Positioning Errors and ZWD';
% ax.XLabel.String = 'Time (Hour)';
% ax.YLabel.String = 'Error (Meter)';
% min1 = t(1) - ((t(end)-t(1))*0.1);
% max1 = t(end);
% ax.XLim = [min1 max1];
% ax.XGrid = 'on';
% ax.YGrid = 'on';
% legend(ax,'North','East','Up','ZWD')
plot(tzd);
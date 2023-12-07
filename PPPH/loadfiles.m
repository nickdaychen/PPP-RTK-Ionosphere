function [rinex,orbitb,anten,clock]=loadfiles(week,reference_station)
base_path_observation = 'D:\data\testData\ObservationFiles\';
num_days =8;
start_day01 = (week - 2191) * 7 + 2;
rinex= char(zeros(num_days, 46) + ' ');
for i = 1:num_days
    day_index = start_day01 + i - 2;
    rinex(i, :) = sprintf('%s%s%03d0.22o', base_path_observation, reference_station, day_index);
end
base_path_clock = 'D:\data\productsCLK\';
% Calculate the starting number based on week value
if week>2237
    clock= char(zeros(num_days, 32) + ' ');
    for i = 1:num_days   
        if i==1
            clock(i, :) = [base_path_clock, sprintf('jpl%d6.clk', week-1)];
        else
            clock(i, :) = [base_path_clock, sprintf('jpl%d%d.clk', week,i-2)]; 
        end
    end
else
    clock= char(zeros(num_days, 36) + ' ');
    for i=1:num_days
        if i==1
            clock(i, :) = [base_path_clock, sprintf('igs%d6.clk_30s', week-1)];
        else
            clock(i, :) = [base_path_clock, sprintf('igs%d%d.clk_30s', week,i-2)]; 
        end
    end
end

orbitb= char(zeros(num_days, 32) + ' ');
for i = 1:num_days
    if week>2237
        if i==1
           orbitb(i, :) = [base_path_clock, sprintf('jpl%d6.sp3', week-1)]; 
        else
           orbitb(i, :) = [base_path_clock, sprintf('jpl%d%d.sp3', week,i-2)]; 
        end
    else
        if i==1
           orbitb(i, :) = [base_path_clock, sprintf('igs%d6.sp3', week-1)]; 
        else
           orbitb(i, :) = [base_path_clock, sprintf('igs%d%d.sp3', week,i-2)]; 
        end
    end
end

anten ='D:\data\atx\igs20.atx';

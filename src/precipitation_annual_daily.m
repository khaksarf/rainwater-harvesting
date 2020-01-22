clc
clear
filename = 'Hourly Rainfall - San Antonio - Annual and Seasonal.xlsx';
dir = '../data/rainfall';


excel_file = xlsread(fullfile(dir, filename ));
year_daily = reshape(excel_file(:,4), 24, []);
year_month = reshape(excel_file(:,1), 24, []);

for i = 1 : 4
    season{i} = year_daily(:, ismember(year_month(1,:), 3*(i-1) + 1 : 3 * i ));
    season_logical{i} = logical(season{i});
    season_pattern{i} = sum(season_logical{i}, 2) / sum(sum(season_logical{i}));
end
season_pattern_all = horzcat(season_pattern{:});
season_pattern_all(isnan(season_pattern_all) ) = 0;

xlswrite(fullfile( dir, ['DSRP_' filename '.xlsx']), season_pattern_all);

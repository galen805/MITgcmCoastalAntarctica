% 400 days of input for wind data, 6-hour period, 4*365, need 1461 in total
% first we use the unit of seconds, then interpolate it onto 6-hour change

%close all; clear; clc;

wind_value_scale_6hr_final_mean = 1;

while abs(wind_value_scale_6hr_final_mean - 0.5) >= 0.005
    
    time_in_sec = 0:1:(365*24*3600);
    wind_omega = 2*pi./((48)*3600);
    wind_value_scale_sec = zeros(size(time_in_sec));
    
    
    for i = 1:length(wind_omega)
        wind_value_scale_sec = wind_value_scale_sec + 1*sin(wind_omega(i).*time_in_sec+...
            -1/4*2*pi);
    end
    
    % figure
    % plot(time_in_sec, wind_value_scale_sec)
    wind_value_scale_sec_mean = mean(wind_value_scale_sec);
    
    
    time_in_6hr = time_in_sec(1:21600:end);
    wind_value_scale_6hr = interp1(time_in_sec, wind_value_scale_sec, time_in_6hr);
    
    % mean(wind_value_scale_6hr)
    
    wind_value_scale_6hr = wind_value_scale_6hr - min(wind_value_scale_6hr);
    wind_value_scale_6hr = wind_value_scale_6hr/(max(wind_value_scale_6hr)-min(wind_value_scale_6hr));
    
    
    % figure
    % plot(time_in_6hr, wind_value_scale_6hr)
    
    wind_value_scale_6hr_final_mean = mean(wind_value_scale_6hr);
    
end

%figure
%plot(time_in_6hr, wind_value_scale_6hr)

save wind_scale_6hr_365days_sporadic_data.mat wind_value_scale_6hr
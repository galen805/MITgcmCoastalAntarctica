load wind_scale_6hr_200days_sporadic_data.mat


time_in_sec = 0:1:(400*24*3600);
time_in_6hr = time_in_sec(1:21600:end);

%figure_pos = [100 100 800 300];

figure
plot(time_in_6hr/24/3600, 40*wind_value_scale_6hr);

xlim([0 100])
xlabel('Time (days)','fontsize',14);
ylabel('40\timesScale (0-1)','fontsize',14);
saveas(gcf,['wind_fluc_100days.png']);


figure
plot(time_in_6hr/24/3600, 40*wind_value_scale_6hr);

xlim([0 20])
xlabel('Time (days)','fontsize',14);
ylabel('40\timesScale (0-1) (m/s)','fontsize',14);
title('V-wind speed variations','fontsize',18);
set(gca, 'fontsize',14);
saveas(gcf,['wind_fluc_20days.png']);

% plot actual maximum model wind

for i=1:length(time_in_6hr)
    Uwindmax(i)=min(uwind(:,:,i),[],'all');
    Vwindmax(i)=max(vwind(:,:,i),[],'all');
end

figure
hold on
plot(time_in_6hr/24/3600,Uwindmax,'linewidth',1.5)
plot(time_in_6hr/24/3600,Vwindmax,'linewidth',1.5)
xlabel('Time (days)','fontsize',14);
xlim([0,200])
ylabel('max wind (m/s)','fontsize',14);
title('max wind speed')
legend('U','V')
saveas(gcf,'max_wind_model.png');

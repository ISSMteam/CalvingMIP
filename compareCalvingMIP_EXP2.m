clear 
close all

nameList = {'AWI (HO)', 'HO'};
fileList = {'EXP2_AWI.mat', 'EXP2_HO.mat'};

Nf = numel(fileList);
% load solutions
for i = 1: Nf
	temp{i} = load(fileList{i});
end

% plot front position
figure('Position', [0, 800, 500, 400])
for i = 1: Nf
	plot(temp{i}.time, mean(temp{i}.distance))
	hold on
end
xlim([0,1000])
ylim([640, 770]*1e3)
set(gcf,'Color','w');
xlabel('Time (a)')
ylabel('Mean Front position (m)')
legend(nameList,'location', 'best')

figure('Position', [0, 800, 500, 400])
for i = 1: Nf
   plot(temp{i}.time, std(temp{i}.distance))
   hold on
end
xlim([0,1000])
xlabel('Time (a)')
set(gcf,'Color','w');
ylabel('Std Front position (m)')
legend(nameList,'location', 'best')

% plot front vel
figure('Position', [0, 800, 500, 400])
for i = 1: Nf
   plot(temp{i}.time, mean(temp{i}.vel))
   hold on
end
xlim([0,1000])
ylim([0,800])
set(gcf,'Color','w');
xlabel('Time (a)')
ylabel('Mean Frontal Velocity (m/a)')
legend(nameList,'location', 'best')

figure('Position', [0, 800, 500, 400])
for i = 1: Nf
   plot(temp{i}.time, std(temp{i}.vel))
   hold on
end
xlim([0,1000])
xlabel('Time (a)')
set(gcf,'Color','w');
ylabel('Std Frontal Velocity (m/a)')
legend(nameList,'location', 'best')

% plot front thickness
figure('Position', [0, 800, 500, 400])
for i = 1: Nf
   plot(temp{i}.time, mean(temp{i}.thickness))
   hold on
end
xlim([0,1000])
ylim([0,350])
set(gcf,'Color','w');
xlabel('Time (a)')
ylabel('Mean Frontal ice thickness (m)')
legend(nameList,'location', 'best')

figure('Position', [0, 800, 500, 400])
for i = 1: Nf
   plot(temp{i}.time, std(temp{i}.thickness))
   hold on
end
xlim([0,1000])
xlabel('Time (a)')
set(gcf,'Color','w');
ylabel('Std Frontal ice thickness (m)')
legend(nameList,'location', 'best')

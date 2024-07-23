clear 
close all

nameList = {'AWI (HO)', 'Dart (SSA, 5k)'};
fileList = {'EXP4_AWI.mat', 'EXP4_Dartmouth.mat'};

Nf = numel(fileList);
% load solutions
for i = 1: Nf
	temp{i} = load(['./Results/', fileList{i}]);
end

% plot front position
figure('Position', [0, 800, 1000, 400])
for i = 1: Nf
	subplot(1,2,1)
	plot(temp{i}.time, mean(temp{i}.distance(1:4,:)))
	hold on
	subplot(1,2,2)
	plot(temp{i}.time, mean(temp{i}.distance(5:8,:)))
	hold on
end
subplot(1,2,1)
xlim([0,1000])
ylim([200, 500]*1e3)
set(gcf,'Color','w');
xlabel('Time (a)')
ylabel('Mean Front position (m)')
title('Caprona')
%legend(nameList,'location', 'best')

subplot(1,2,2)
xlim([0,1000])
ylim([450, 750]*1e3)
set(gcf,'Color','w');
xlabel('Time (a)')
ylabel('Mean Front position (m)')
legend(nameList,'location', 'best')
title('Halbrane')

% plot front position
figure('Position', [0, 800, 1000, 400])
for i = 1: Nf
   subplot(1,2,1)
   plot(temp{i}.time, std(temp{i}.distance(1:4,:)))
   hold on
   subplot(1,2,2)
   plot(temp{i}.time, std(temp{i}.distance(5:8,:)))
   hold on
end
subplot(1,2,1)
xlim([0,1000])
ylim([0,7000])
set(gcf,'Color','w');
xlabel('Time (a)')
ylabel('Std Front position (m)')
title('Caprona')
%legend(nameList,'location', 'best')

subplot(1,2,2)
xlim([0,1000])
ylim([0,7000])
set(gcf,'Color','w');
xlabel('Time (a)')
ylabel('Std Front position (m)')
legend(nameList,'location', 'best')
title('Halbrane')

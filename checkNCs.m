clear
close all

EXP = 4;
folder = './Results/2_5kmResults_extrap/';
ExpName = 'CalvingMIP_EXP4_ISSM_SSA_Dartmouth.nc';
% load ice front results from NCs
if EXP >= 3
	pnames = {'Cap', 'Hal'};
	suffix = {'A', 'B', 'C', 'D'};
else
	error('EXP 1 & 2 not supported yet')
end

ncfilename = [folder, ExpName];
% time
if EXP == 4 
	time = ncread(ncfilename, 'Time1');
end
for i = 1:numel(pnames)
	for j = 1:numel(suffix)
		xc(i,j, :) = ncread(ncfilename, ['xcf', pnames{i} suffix{j}]);
		yc(i,j, :) = ncread(ncfilename, ['ycf', pnames{i} suffix{j}]);
		vx(i,j, :) = ncread(ncfilename, ['xvelmeancf', pnames{i} suffix{j}]);
		vy(i,j, :) = ncread(ncfilename, ['yvelmeancf', pnames{i} suffix{j}]);
		thickness(i,j, :) = ncread(ncfilename, ['lithkcf', pnames{i} suffix{j}]);
	end
end
% plot
for i = 1:numel(pnames)
	figure('position',[0,1000,800,600])
	for j = 1:numel(suffix)
		subplot(2,2,1)
		plot(time, squeeze(abs(xc(i,j,:))))
		hold on
		subplot(2,2,2)
		plot(time, squeeze(abs(yc(i,j,:))))
		hold on
		subplot(2,2,3)
		plot(time, squeeze(sqrt(vx(i,j,:).^2+vy(i,j,:).^2)))
		hold on
		subplot(2,2,4)
		plot(time, squeeze(abs(thickness(i,j,:))))
		hold on
	end
end

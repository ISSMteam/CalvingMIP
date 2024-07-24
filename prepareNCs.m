clear
close all


% Setting {{{
glacier = 'Thule';
projPath = './';
folderList = {'../Thule/Models/20230705_EXP3_res_5000/', '../Thule/Models/20240722_EXP4_res_5000/'};
steps = [1];
authorList = 'Gong Cheng(gong.cheng@dartmouth.edu), Helene Seroussi, Mathieu Morlighem';
outputFolder = './Results/5kmResults/';
%}}}
% Loading models {{{
stepName = 'Transient';
Ndata = length(folderList);
parfor i = 1:Ndata
	org{i}=organizer('repository', [projPath, folderList{i}], 'prefix', ['Model_' glacier '_'], 'steps', steps);
	if perform(org{i}, 'Transient')
		disp(['---- Loading the model from ', folderList{i}]);
		mdList{i} = loadmodel(org{i}, [stepName]);
	end
end
%}}}
% copy EXP3 final results to EXP4 model if needed {{{
if ~isfield(mdList{2}.results, 'InitialSolution')
	mdList{2}.results.InitialSolution.GroundedArea = mdList{1}.results.TransientSolution(end).GroundedArea;
	mdList{2}.results.InitialSolution.FloatingArea = mdList{1}.results.TransientSolution(end).FloatingArea;
	mdList{2}.results.InitialSolution.IceVolume = mdList{1}.results.TransientSolution(end).IceVolume;
	mdList{2}.results.InitialSolution.IceVolumeAboveFloatation = mdList{1}.results.TransientSolution(end).IceVolumeAboveFloatation;
	mdList{2}.results.InitialSolution.IcefrontMassFluxLevelset = mdList{1}.results.TransientSolution(end).IcefrontMassFluxLevelset;
	mdList{2}.results.InitialSolution.GroundinglineMassFlux = mdList{1}.results.TransientSolution(end).GroundinglineMassFlux;

	savemodel(org{2}, mdList{2});
end
%}}}
% Create NetCDF {{{ 
results3 = ModelToNetCDF(mdList{1}, 'directoryname', outputFolder, 'EXP', 3, 'author', authorList, 'institution', 'Dartmouth', 'flowequation', 'SSA');
results4 = ModelToNetCDF(mdList{2}, 'directoryname', outputFolder, 'EXP', 4, 'author', authorList, 'institution', 'Dartmouth', 'flowequation', 'SSA');
% }}}

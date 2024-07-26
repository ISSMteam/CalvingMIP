clear
close all


% Setting {{{
glacier = 'Thule';
projPath = './';
%folderList = {'../Thule/Models/20230705_EXP3_res_5000/', '../Thule/Models/20240722_EXP4_res_5000/'};
%outputFolder = './Results/5kmResults/';
folderList = {'../Thule/Models/20230714_EXP3_res_2500//', '../Thule/Models/20240722_EXP4_res_2500/'};
outputFolder = './Results/2_5kmResults_extrap/';
steps = [1];
authorList = 'Gong Cheng(gong.cheng@dartmouth.edu), Helene Seroussi, Mathieu Morlighem';
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
% Generate profiles if not exist {{{
profilefolder = './Profiles/';
cprofilename = [profilefolder, '/Caprona_Profiles.csv'];
hprofilename = [profilefolder, '/Halbrane_Profiles.csv'];
if (~exist(cprofilename, 'file') | (~exist(hprofilename, 'file')) )
	createProfiles('EXP4', profilefolder)
end
profilefolder = './Profiles/finer/';
cprofilename = [profilefolder, '/Caprona_Profiles.csv'];
hprofilename = [profilefolder, '/Halbrane_Profiles.csv'];
if ~exist(profilefolder, 'dir')
	mkdir(profilefolder)
end
if (~exist(cprofilename, 'file') | (~exist(hprofilename, 'file')) )
	createProfiles('EXP4', profilefolder, 1000)
end
%}}}
% Create NetCDF {{{ 
results3 = ModelToNetCDF(mdList{1}, 'directoryname', outputFolder, 'EXP', 3, 'author', authorList, 'institution', 'Dartmouth', 'flowequation', 'SSA');
results4 = ModelToNetCDF(mdList{2}, 'directoryname', outputFolder, 'EXP', 4, 'author', authorList, 'institution', 'Dartmouth', 'flowequation', 'SSA');
% }}}

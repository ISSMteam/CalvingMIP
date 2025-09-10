clear
close all

% Setting {{{
glacier = 'Thule';
projPath = './';
%folderList = {'2025_HO_Exp3/'}; 
folderList = {'2025_HO_Exp3/', '2025_HO_Exp4/'};
%}}}
% Loading data {{{
for i=1:length(folderList)
	org = organizer('repository', [projPath, 'Models/', folderList{i}], 'prefix', ['Model_' glacier '_'], 'steps', 1);
	disp(['Load the model from ', folderList{i}])
	if perform(org, 'Transient') % {{{
		md = loadmodel(org, 'HO');
		% project solutions onto the 2D mesh before extrusion  
		md2 = collapse(md);

		% manually remove the unused data
		md2.results.Transient1 = [];
		md2.results.Transient2 = [];
		% be careful
		md2.results.TransientSolution = [];
		%project solutions to md2
		for i=1:numel(md.results.TransientSolution)
			% scalars
			md2.results.TransientSolution(i).time = md.results.TransientSolution(i).time;
			md2.results.TransientSolution(i).IceVolume = md.results.TransientSolution(i).IceVolume;
			md2.results.TransientSolution(i).IceVolumeAboveFloatation = md.results.TransientSolution(i).IceVolumeAboveFloatation;
			md2.results.TransientSolution(i).GroundedArea = md.results.TransientSolution(i).GroundedArea;
			md2.results.TransientSolution(i).FloatingArea = md.results.TransientSolution(i).FloatingArea;
			md2.results.TransientSolution(i).IcefrontMassFluxLevelset = md.results.TransientSolution(i).IcefrontMassFluxLevelset;
			md2.results.TransientSolution(i).GroundinglineMassFlux = md.results.TransientSolution(i).GroundinglineMassFlux;

			% vectors
			md2.results.TransientSolution(i).Vx = project2d(md,md.results.TransientSolution(i).Vx,md.mesh.numberoflayers);
			md2.results.TransientSolution(i).Vy = project2d(md,md.results.TransientSolution(i).Vy,md.mesh.numberoflayers);
			md2.results.TransientSolution(i).Thickness = project2d(md,md.results.TransientSolution(i).Thickness,md.mesh.numberoflayers);
			md2.results.TransientSolution(i).Base = project2d(md,md.results.TransientSolution(i).Base,md.mesh.numberoflayers);
			md2.results.TransientSolution(i).MaskIceLevelset = project2d(md,md.results.TransientSolution(i).MaskIceLevelset,md.mesh.numberoflayers);
			md2.results.TransientSolution(i).MaskOceanLevelset = project2d(md,md.results.TransientSolution(i).MaskOceanLevelset,md.mesh.numberoflayers);
		end
		% save md2
		savemodel(org, md2);
	end %}}}
end
%}}}

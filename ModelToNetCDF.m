function results = ModelToNetCDF(md, varargin)
	%Create netcdf files for an experiment
	%Check inputs {{{
	if nargout>1
		help runme
		error('runme error message: bad usage');
	end
	%recover options
	options=pairoptions(varargin{:});
	% }}}
	%GET directoryname : './/', directory where the outputs are placed {{{
	directoryname = getfieldvalue(options,'directoryname','.//');
	if exist(directoryname)~=7,
		disp(['directory ' directoryname ' does not exist']);
		mkdir(directoryname);
	end
	% }}}
	%GET EXP : 1, EXP1-EXP4 are currently supported{{{
	expnum = getfieldvalue(options,'EXP', 1);
	expstr = ['EXP', num2str(expnum)];
	% }}}
	%GET modelname : 'ISSM' {{{
	modelname = getfieldvalue(options,'modelname', 'ISSM');
	% }}}
	%GET flowequation : 'SSA' {{{
	flowequation = getfieldvalue(options, 'flowequation', 'SSA');
	% }}}
	%GET institution : 'Dartmouth' {{{
	institution = getfieldvalue(options, 'institution', 'Dartmouth');
	% }}}
	%GET author : '' {{{
	author = getfieldvalue(options, 'author', '');
	% }}}

	%Field variables {{{
	disp(['loading field variables..']);

	if strcmp(expstr,'EXP1') | strcmp(expstr,'EXP3')
		results.Time1 = [0]; % Time1
		results.timegrid=[0];
	elseif strcmp(expstr,'EXP2') | strcmp(expstr,'EXP4')
		results.Time1 = [0:1:1000]; % Time1
		results.timegrid = [0:100:1000]; % Time100
	end

	xmin=-800000; xmax=800000; dx=5000;
	ymin=-800000; ymax=800000; dy=5000;

	results.gridx = xmin:dx:xmax;
	results.gridy = ymin:dy:ymax;
	sizex = numel(results.gridx);
	sizey = numel(results.gridy);

	% allocate matrices
	results.vxmean = nan(sizex,sizey,length(results.timegrid)); %y,x,time
	results.vymean = nan(sizex,sizey,length(results.timegrid)); %y,x,time
	results.thickness= nan(sizex,sizey,length(results.timegrid)); %y,x,time
	results.mask = nan(sizex,sizey,length(results.timegrid)); %y,x,time
	results.base = nan(sizex,sizey,length(results.timegrid)); %y,x,time
	results.icemask = nan(sizex,sizey,length(results.timegrid)); %y,x,time
	results.oceanmask = nan(sizex,sizey,length(results.timegrid)); %y,x,time

	% get the model solution into a matrix
	x = md.mesh.x;
	y = md.mesh.y;
	index = md.mesh.elements;

	if strcmp(expstr,'EXP1') | strcmp(expstr,'EXP3')
		time = 0; 
		vx = md.results.TransientSolution(end).Vx;
		vy = md.results.TransientSolution(end).Vy;
		thickness = md.results.TransientSolution(end).Thickness;
		base = md.results.TransientSolution(end).Base;
		icemask = md.results.TransientSolution(end).MaskIceLevelset;
		oceanmask = md.results.TransientSolution(end).MaskOceanLevelset;
	elseif strcmp(expstr,'EXP2') | strcmp(expstr,'EXP4')
		mdT = [md.results.TransientSolution(:).time];
		% find the indices in model results correspond to results.Time1
		[~, indd] = find(mdT==results.Time1');

		% take all the request time points
		time = [0, md.results.TransientSolution(indd).time]; % time start at 0, add initial conditions at the beginning
		vx = [[md.initialization.vx], md.results.TransientSolution(indd).Vx];
		vy = [[md.initialization.vy], md.results.TransientSolution(indd).Vy];
		thickness = [[md.geometry.thickness], md.results.TransientSolution(indd).Thickness];
		base = [[md.geometry.base], md.results.TransientSolution(indd).Base];
		icemask = [[md.mask.ice_levelset], md.results.TransientSolution(indd).MaskIceLevelset];
		oceanmask = [[md.mask.ocean_levelset], md.results.TransientSolution(indd).MaskOceanLevelset];
	end
	% save the no mask results, for calving front data process
	vx_nomask = vx;
	vy_nomask = vy;
	thickness_nomask = thickness;

	% set vel/thk in the open ocean to NaN
	vx(icemask>0) = NaN;
	vy(icemask>0) = NaN;
	thickness(icemask>0) = NaN;
	base(icemask>0) = NaN;

	% get these data on the time100 grid
	for i = 1:length(results.timegrid)
		if strcmp(expstr,'EXP1') | strcmp(expstr,'EXP3')
			tid = 1;
		elseif strcmp(expstr,'EXP2') | strcmp(expstr,'EXP4')
			[~,tid] = min(abs( results.Time1((i-1)*100+1)-time));
		end
		results.vxmean(:, :, i) = transpose(InterpFromMeshToGrid(index, x, y, vx(:, tid), results.gridx, results.gridy, NaN));
		results.vymean(:, :, i) = transpose(InterpFromMeshToGrid(index, x, y, vy(:, tid), results.gridx, results.gridy, NaN));
		results.thickness(:, :, i) = transpose(InterpFromMeshToGrid(index, x, y, thickness(:, tid), results.gridx, results.gridy, NaN));
		results.base(:, :, i) = transpose(InterpFromMeshToGrid(index, x, y, base(:, tid), results.gridx, results.gridy, NaN));
		results.icemask(:, :, i) = transpose(InterpFromMeshToGrid(index, x, y, icemask(:, tid), results.gridx, results.gridy, NaN));
		results.oceanmask(:, :, i) = transpose(InterpFromMeshToGrid(index, x, y, oceanmask(:, tid), results.gridx, results.gridy, NaN));
	end

	% mask: grounded=1, floating=2, open ocean=3
	mask = results.icemask;
	mask(results.icemask<0) = 1;		% grounded and floating
	mask(results.icemask>0) = 3;		% open ocean
	mask((results.oceanmask<0) & (results.icemask<0)) = 2;
	results.mask = convertLevelsetsToCalvingMIPMask(results.icemask, results.oceanmask);

	% bed is static
	results.bed = transpose(InterpFromMeshToGrid(index, x, y, md.geometry.bed, results.gridx, results.gridy, NaN));

	% Total area of grounded and floating ice in the four quadrants
	results.areafourquad = computeFourQuadrantIceArea(md, icemask);
	%}}}
	%Scalar variables %{{{
	disp(['computing scalar variables..']);

	if strcmp(expstr,'EXP1') | strcmp(expstr,'EXP3')
		results.groundedarea = md.results.TransientSolution(end).GroundedArea; % m^2
		results.floatingarea = md.results.TransientSolution(end).FloatingArea; % m^2
		results.mass = md.results.TransientSolution(end).IceVolume*md.materials.rho_ice; % kg
		results.massaf = md.results.TransientSolution(end).IceVolumeAboveFloatation*md.materials.rho_ice; % kg
		results.totalflux_calving = md.results.TransientSolution(end).IcefrontMassFluxLevelset*10^12; % Gt/yr -> kg/yr
		results.totalflux_groundingline = md.results.TransientSolution(end).GroundinglineMassFlux*10^12;  %Gt/yr -> kg/yr
	elseif strcmp(expstr,'EXP2') | strcmp(expstr,'EXP4')
		% add initial values of these quantities from EXP3
		results.groundedarea = [md.results.InitialSolution.GroundedArea, md.results.TransientSolution(indd).GroundedArea]; % m^2
		results.floatingarea = [md.results.InitialSolution.FloatingArea, md.results.TransientSolution(indd).FloatingArea]; % m^2
		results.mass = [md.results.InitialSolution.IceVolume, md.results.TransientSolution(indd).IceVolume]*md.materials.rho_ice; % kg
		results.massaf = [md.results.InitialSolution.IceVolumeAboveFloatation, md.results.TransientSolution(indd).IceVolumeAboveFloatation]*md.materials.rho_ice; % kg
		results.totalflux_calving = [md.results.InitialSolution.IcefrontMassFluxLevelset, md.results.TransientSolution(indd).IcefrontMassFluxLevelset]*10^12; % Gt/yr -> kg/yr
		results.totalflux_groundingline = [md.results.InitialSolution.GroundinglineMassFlux, md.results.TransientSolution(indd).GroundinglineMassFlux]*10^12;  %Gt/yr -> kg/yr
	else
		error('not implemented yet');	
	end
	%}}}
	%Profile variables {{{
	disp(['loading profile variables..']);
	if strcmp(expstr,'EXP3') | strcmp(expstr,'EXP4')
		nameList = {'A', 'B', 'C', 'D'};
		prfList = {'Caprona', 'Halbrane'};

		for ip = 1:numel(prfList)
			% loof through Caprona and Halbrane
			P = readtable(['./Profiles/', prfList{ip}, '_Profiles.csv']);
			Pc = readtable(['./Profiles/finer/', prfList{ip}, '_Profiles.csv']);
			suffixname = [prfList{ip}, '_'];
			disp(['  Projecting solutions onto ' suffixname, ' profiles'])
			for i = 1:numel(nameList)
				pfx = P.([suffixname, 'Profile_', nameList{i}, '_X']);
				pfy = P.([suffixname, 'Profile_', nameList{i}, '_Y']);

				% use a temporary profile variable to hold all the info
				pf = [];
				% distance from the start
				pf.distance = P.([suffixname, 'Profile_', nameList{i}, '_S']);

				% project ice thickness, icemask, vx and vy
				pf.thickness = InterpFromMeshToMesh2d(index, x, y, thickness, pfx, pfy);
				pf.vx = InterpFromMeshToMesh2d(index, x, y, vx, pfx, pfy);
				pf.vy = InterpFromMeshToMesh2d(index, x, y, vy, pfx, pfy);
				pf.icemask = InterpFromMeshToMesh2d(index, x, y, icemask, pfx, pfy);
				pf.oceanmask = InterpFromMeshToMesh2d(index, x, y, oceanmask, pfx, pfy);
				pf.mask = convertLevelsetsToCalvingMIPMask(pf.icemask, pf.oceanmask);

				% get front 
				pfc = [];
				pfc.distance = Pc.([suffixname, 'Profile_', nameList{i}, '_S']);
				pfc.x = Pc.([suffixname, 'Profile_', nameList{i}, '_X']);
				pfc.y = Pc.([suffixname, 'Profile_', nameList{i}, '_Y']);
				pfc.thickness = InterpFromMeshToMesh2d(index, x, y, thickness, pfc.x, pfc.y);
				pfc.vx = InterpFromMeshToMesh2d(index, x, y, vx, pfc.x, pfc.y);
				pfc.vy = InterpFromMeshToMesh2d(index, x, y, vy, pfc.x, pfc.y);
				pfc.icemask = InterpFromMeshToMesh2d(index, x, y, icemask, pfc.x, pfc.y);

				pf.front = getFrontFromProfiles_extrap(pfc);
				results.profiles.([suffixname, nameList{i}]) = pf;
			end
		end
	else
		error('not implemented yet');	
	end
	%}}}
	%Create netCDF{{{
	ExpName=[directoryname '/CalvingMIP_' expstr '_' modelname '_' flowequation '_' institution '.nc'];
	disp(['Create netCDF files for the ', expstr, ' at ', ExpName])

	% Delete old file
	if exist(ExpName,'file') == 2
		delete(ExpName);
	end

	% Time{{{
	if strcmp(expstr,'EXP1') | strcmp(expstr,'EXP3')
		timeStr = 'Time';
		timeExtraStr = 'Time';
	elseif strcmp(expstr,'EXP2') | strcmp(expstr,'EXP4')
		timeStr = 'Time1';
		timeExtraStr = 'Time100';

		% additional variable
		nccreate(ExpName, timeStr,'Dimensions',{timeStr numel(results.Time1)})
		ncwrite(ExpName,timeStr,results.Time1)
		ncwriteatt(ExpName, timeStr,'units','a');
	end
	nccreate(ExpName, timeExtraStr, 'Dimensions', {timeExtraStr, numel(results.timegrid)})
	ncwrite(ExpName, timeExtraStr, results.timegrid)
	ncwriteatt(ExpName, timeExtraStr,'units','a');
	%}}}
	% Spatial dimensions {{{
	nccreate(ExpName, 'X', 'Dimensions', {'X' 321})
	nccreate(ExpName, 'Y', 'Dimensions', {'Y' 321})
	ncwrite(ExpName, 'X', results.gridx)
	ncwrite(ExpName, 'Y', results.gridy)
	ncwriteatt(ExpName,'X','units','m')
	ncwriteatt(ExpName,'Y','units','m') 
	%}}}
	% Fields: velocity, thickness, masks{{{
	nccreate(ExpName,'xvelmean','Dimensions',{'X' 321 'Y' 321 timeExtraStr numel(results.timegrid)},'FillValue',nan)
	nccreate(ExpName,'yvelmean','Dimensions',{'X' 321 'Y' 321 timeExtraStr numel(results.timegrid)},'FillValue',nan)
	nccreate(ExpName,'lithk','Dimensions',{'X' 321 'Y' 321 timeExtraStr numel(results.timegrid)},'FillValue',nan)
	nccreate(ExpName,'mask','Dimensions',{'X' 321 'Y' 321 timeExtraStr numel(results.timegrid)})
	nccreate(ExpName,'topg','Dimensions',{'X' 321 'Y' 321})

	ncwrite(ExpName, 'xvelmean', results.vxmean)
	ncwrite(ExpName, 'yvelmean', results.vymean)
	ncwrite(ExpName, 'lithk', results.thickness)
	ncwrite(ExpName, 'mask', results.mask)
	ncwrite(ExpName, 'topg', results.bed)

	ncwriteatt(ExpName,'xvelmean','units','m/a')
	ncwriteatt(ExpName,'yvelmean','units','m/a')
	ncwriteatt(ExpName,'lithk','units','m');
	ncwriteatt(ExpName,'mask','flag_values','1, 2, 3');
	ncwriteatt(ExpName,'topg','units','m');

	ncwriteatt(ExpName,'xvelmean','Standard_name','land_ice_vertical_mean_x_velocity')
	ncwriteatt(ExpName,'yvelmean','Standard_name','land_ice_vertical_mean_y_velocity')
	ncwriteatt(ExpName,'lithk','Standard_name','land_ice_thickness');
	ncwriteatt(ExpName,'mask','flag_meanings','1=grounded ice, 2=floating ice, 3=open ocean');
	ncwriteatt(ExpName,'topg','Standard_name','bedrock_altitude');
	%}}}
	% Scalar variables {{{
	nccreate(ExpName,'iareafl',		'Dimensions',{timeStr numel(results.Time1)})
	nccreate(ExpName,'iareagr',		'Dimensions',{timeStr numel(results.Time1)})
	nccreate(ExpName,'lim',				'Dimensions',{timeStr numel(results.Time1)})
	nccreate(ExpName,'limnsw',			'Dimensions',{timeStr numel(results.Time1)})
	nccreate(ExpName,'tendlicalvf',	'Dimensions',{timeStr numel(results.Time1)})
	nccreate(ExpName,'tendligroundf','Dimensions',{timeStr numel(results.Time1)})

	ncwrite(ExpName,'iareafl',results.floatingarea)
	ncwrite(ExpName,'iareagr',results.groundedarea)
	ncwrite(ExpName,'lim',results.mass)
	ncwrite(ExpName,'limnsw',results.massaf)
	ncwrite(ExpName,'tendlicalvf',results.totalflux_calving)
	ncwrite(ExpName,'tendligroundf', results.totalflux_groundingline)

	ncwriteatt(ExpName,'iareafl','units','m^2')
	ncwriteatt(ExpName,'iareagr','units','m^2')
	ncwriteatt(ExpName,'lim','units','kg')
	ncwriteatt(ExpName,'limnsw','units','kg')
	ncwriteatt(ExpName,'tendlicalvf','units','kg/a');
	ncwriteatt(ExpName,'tendligroundf','units','kg/a');

	ncwriteatt(ExpName,'iareafl','Standard_name','grounded_ice_sheet_area')
	ncwriteatt(ExpName,'iareagr','Standard_name','floating_ice_shelf_area')
	ncwriteatt(ExpName,'lim','Standard_name','land_ice_mass')
	ncwriteatt(ExpName,'limnsw','Standard_name','land_ice_mass_not_displacing_sea_water')
	ncwriteatt(ExpName,'tendlicalvf','Standard_name','tendency_of_land_ice_mass_due_to_calving');
	ncwriteatt(ExpName,'tendligroundf','Standard_name','tendency_of_grounded_ice_mass');
	%}}}
	% Four quadrant variables {{{
	nccreate(ExpName,'iareatotalNW',	'Dimensions',{timeStr numel(results.Time1)})
	nccreate(ExpName,'iareatotalNE',	'Dimensions',{timeStr numel(results.Time1)})
	nccreate(ExpName,'iareatotalSW',	'Dimensions',{timeStr numel(results.Time1)})
	nccreate(ExpName,'iareatotalSE',	'Dimensions',{timeStr numel(results.Time1)})

	ncwrite(ExpName,'iareatotalNW',results.areafourquad(2,:))
	ncwrite(ExpName,'iareatotalNE',results.areafourquad(1,:))
	ncwrite(ExpName,'iareatotalSW',results.areafourquad(3,:))
	ncwrite(ExpName,'iareatotalSE',results.areafourquad(4,:))

	ncwriteatt(ExpName,'iareatotalNW','units','m^2')
	ncwriteatt(ExpName,'iareatotalNE','units','m^2')
	ncwriteatt(ExpName,'iareatotalSW','units','m^2')
	ncwriteatt(ExpName,'iareatotalSE','units','m^2')

	ncwriteatt(ExpName,'iareatotalNW','Standard_name','total_ice_area_NorthWest')
	ncwriteatt(ExpName,'iareatotalNE','Standard_name','total_ice_area_NorthEast')
	ncwriteatt(ExpName,'iareatotalSW','Standard_name','total_ice_area_SouthWest')
	ncwriteatt(ExpName,'iareatotalSE','Standard_name','total_ice_area_SouthEast')
	%}}}
	% Profiles {{{
	if strcmp(expstr,'EXP3') | strcmp(expstr,'EXP4')
		fullPre = {'Caprona', 'Halbrane'};		% prefix of the full profile name
		shortPre = {'Cap', 'Hal'};						% prefix of the short profile name
	else
		error('Not implemented')
	end

	for p = 1:numel(fullPre)
		for i = 1:numel(nameList)
			% Profiles
			pf = results.profiles.([fullPre{p}, '_', nameList{i}]);
			% short and long names for the netCDF file
			sN = [shortPre{p}, nameList{i}];
			fN = [fullPre{p}, ' ', nameList{i}];
			% write varibles to nc files
			if strcmp(expstr,'EXP1') | strcmp(expstr,'EXP3')
				% no time dimension for Exp1 and EXP3
				nccreate(ExpName, ['lithk', sN], 'Dimensions', {fN numel(pf.distance)});
				nccreate(ExpName, ['s', sN],'Dimensions',{fN numel(pf.distance)});
				nccreate(ExpName, ['xvelmean', sN], 'Dimensions', {fN numel(pf.distance)}, 'FillValue', nan);
				nccreate(ExpName, ['yvelmean', sN], 'Dimensions', {fN numel(pf.distance)}, 'FillValue', nan);
				nccreate(ExpName, ['mask', sN], 'Dimensions', {fN numel(pf.distance)});
			elseif strcmp(expstr,'EXP2') | strcmp(expstr,'EXP4')
				% add time dimension to Exp2 and EXP4
				nccreate(ExpName, ['lithk', sN], 'Dimensions', {fN numel(pf.distance) timeStr numel(results.Time1)} );
				nccreate(ExpName, ['s', sN],'Dimensions',{fN numel(pf.distance) timeStr numel(results.Time1)});
				nccreate(ExpName, ['xvelmean', sN], 'Dimensions', {fN numel(pf.distance) timeStr numel(results.Time1)}, 'FillValue', nan);
				nccreate(ExpName, ['yvelmean', sN], 'Dimensions',{ fN numel(pf.distance) timeStr numel(results.Time1)}, 'FillValue', nan);
				nccreate(ExpName, ['mask', sN], 'Dimensions', {fN numel(pf.distance) timeStr numel(results.Time1)});
			else
				error('Not implemented')
			end

			ncwrite(ExpName, ['lithk', sN], pf.thickness)
			ncwrite(ExpName, ['s', sN], pf.distance)
			ncwrite(ExpName, ['xvelmean', sN], pf.vx)
			ncwrite(ExpName, ['yvelmean', sN], pf.vy)
			ncwrite(ExpName, ['mask', sN], pf.mask)

			ncwriteatt(ExpName, ['lithk', sN], 'units','m');
			ncwriteatt(ExpName, ['s', sN], 'units', 'm');
			ncwriteatt(ExpName, ['xvelmean' sN], 'units', 'm/a');
			ncwriteatt(ExpName, ['yvelmean' sN], 'units', 'm/a');
			ncwriteatt(ExpName, ['mask', sN], 'flag_values', '1, 2, 3');

			ncwriteatt(ExpName, ['lithk', sN], 'Standard_name',['land_ice_thickness_along_', fN]);
			ncwriteatt(ExpName, ['s', sN], 'Standard_name', ['distance_along_', fN]);
			ncwriteatt(ExpName, ['xvelmean', sN], 'Standard_name', ['land_ice_vertical_mean_x_velocity_along_', fN]);
			ncwriteatt(ExpName, ['yvelmean', sN], 'Standard_name', ['land_ice_vertical_mean_y_velocity_along_', fN]);
			ncwriteatt(ExpName, ['mask', sN], 'flag_meanings', '1=grounded ice, 2=floating ice, 3=open ocean');

			% add Calving front data
			nccreate(ExpName, ['xcf', sN], 'Dimensions', {timeStr numel(results.Time1)})
			nccreate(ExpName, ['ycf', sN], 'Dimensions', {timeStr numel(results.Time1)})
			nccreate(ExpName, ['xvelmeancf', sN], 'Dimensions', {timeStr numel(results.Time1)})
			nccreate(ExpName, ['yvelmeancf', sN], 'Dimensions', {timeStr numel(results.Time1)})
			nccreate(ExpName, ['lithkcf', sN], 'Dimensions', {timeStr numel(results.Time1)})

			ncwrite(ExpName, ['xcf', sN], pf.front.x)
			ncwrite(ExpName, ['ycf', sN], pf.front.y)
			ncwrite(ExpName, ['xvelmeancf', sN], pf.front.vx)
			ncwrite(ExpName, ['yvelmeancf', sN], pf.front.vy)
			ncwrite(ExpName, ['lithkcf', sN], pf.front.thickness)

			ncwriteatt(ExpName, ['xcf', sN], 'units', 'm')
			ncwriteatt(ExpName, ['ycf', sN], 'units', 'm')
			ncwriteatt(ExpName, ['xvelmeancf', sN], 'units', 'm/a')
			ncwriteatt(ExpName, ['yvelmeancf', sN], 'units', 'm/a')
			ncwriteatt(ExpName, ['lithkcf', sN], 'units', 'm')

			ncwriteatt(ExpName, ['xcf', sN], 'Standard_name', ['x_calving_front_on_', fN])
			ncwriteatt(ExpName, ['ycf', sN], 'Standard_name', ['y_calving_front_on_', fN])
			ncwriteatt(ExpName, ['xvelmeancf', sN], 'Standard_name', ['land_ice_vertical_mean_x_velocity_at_calving_front_on_', fN])
			ncwriteatt(ExpName, ['yvelmeancf', sN], 'Standard_name', ['land_ice_vertical_mean_y_velocity_at_calving_front_on_', fN])
			ncwriteatt(ExpName, ['lithkcf', sN], 'Standard_name', ['land_ice_thickness_at_calving_front_on_', fN])
		end
	end
	%}}}
	ncdisp(ExpName)
	%}}}

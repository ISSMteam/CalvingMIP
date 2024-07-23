clear
close all

% Setting 
projPath = './';
EXP = 'EXP4';
filename = '20240722_EXP4_res_5000/Model_Thule_Transient'; 
savename = 'EXP4_Dartmouth_update_calving';
profilefolder = [projPath, '/Profiles/'];

% Loading model
org=organizer('repository', ['../Thule/Models/'], 'prefix', '', 'steps', 0);
%org=organizer('repository', [projPath, 'Models/'], 'prefix', '', 'steps', 0);
disp(['Load the model from ', filename])
md = loadmodel(org, filename);

% create Profiles folder
if ~exist(profilefolder, 'dir')
	mkdir(profilefolder)
end
% create profiles if they don't exist
% get the solutions along the profiles
if strcmp(EXP, 'EXP2')
	suffixname = 'Circle';
	profilename = [profilefolder, '/', suffixname, '_Profiles.csv'];
	if ~exist(profilename, 'file')
		createProfiles(EXP, profilefolder)
	end
	P = readtable(profilename);
	profiles = project2Profiles(md, P, suffixname);
elseif strcmp(EXP, 'EXP4')
	% Caprona
	suffixname = 'Caprona';
	profilename = [profilefolder, '/', suffixname, '_Profiles.csv'];
	if ~exist(profilename, 'file')
		createProfiles(EXP, profilefolder)
	end
	P = readtable(profilename);
	profiles = project2Profiles(md, P, suffixname);
	% Halbrane
	suffixname = 'Halbrane';
	profilename = [profilefolder, '/', suffixname, '_Profiles.csv'];
	if ~exist(profilename, 'file')
		createProfiles(EXP, profilefolder)
	end
	Q = readtable(profilename);
	Q_profiles = project2Profiles(md, Q, suffixname);
	% merge Q_profiles to profiles
	profiles = mergeProfiles(profiles, Q_profiles);
end

% process
names = fieldnames(profiles);
pfnames = names(2:end);
time = profiles.time;
Nt = numel(time);
Np = numel(pfnames);

distance = zeros(Np, Nt);
thickness = zeros(Np, Nt);
vel = zeros(Np, Nt);

for j = 1:Np
	front = getFrontFromProfiles(profiles.(pfnames{j}));
	distance(j,:) = front.distance;
	thickness(j,:) = front.thickness;
	vel(j,:) = front.vel;
end
% save 
if ~exist([projPath, '/Results/'], 'dir')
	mkdir([projPath, '/Results/'])
end
save([projPath, '/Results/', savename, '.mat'], 'time', 'distance', 'thickness', 'vel')

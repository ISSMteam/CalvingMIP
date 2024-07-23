Matlab scripts for analyzing CalvingMIP Experiment results produced by ISSM.

see [CalvingMIP](https://github.com/JRowanJordan/CalvingMIP).

- `processTransientSolutions.m` : the main script to process ISSM results `md`,  project the transient solutions to the profiles, including ice front positions, ice velocity, thickness, etc.
- `convertLevelsetsToCalvingMIPMask.m` : a function to convert ISSM `ice_levelset` and `ocean_levelset` to the mask format in CalvingMIP.
- `createProfiles.m` : generate profiles for CalvingMIP EXP1-4.
- `getFrontFromProfiles.m` : compute the ice front position from the given `profile`, using `profile.distance` and `profile.icemask`.
- `mergeProfiles.m` : merge two profile structs.
- `project2Profiles.m` : project ISSM model solutions to the given profiles.

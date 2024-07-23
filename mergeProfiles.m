function P = mergeProfiles(P,Q)
	% merge all the field from Q, to P, for duplicate names, use P's
	qfnames = fieldnames(Q);
	for i = 1: numel(qfnames)
		if ~isfield(P, qfnames{i})
			P.(qfnames{i}) = Q.(qfnames{i});
		end
	end
end


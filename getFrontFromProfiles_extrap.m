function front=getFrontFromProfiles_extrap(p)
	% find exact ice front position
	x = p.distance;
	x0 = interpZeroPos(x, p.icemask)';
	% interpolate the other variables at x0
	fnames = fieldnames(p);
	for j = 1:numel(fnames)
		if (~strcmp(fnames{j}, 'distance'))
			temp = zeros(size(x0));
			if size(p.(fnames{j}), 2) == numel(x0)
				for t = 1:numel(x0)
					d_temp = p.(fnames{j})(:,t);
					pos = ~isnan(d_temp);
					x_temp = x(pos);
					temp(t) = interp1(x_temp,d_temp(pos), x0(t),'linear','extrap');
				end
			else
				temp = interp1(x,p.(fnames{j}), x0,'linear','extrap');
			end
			front.(fnames{j}) = temp;
		end
	end
	front.distance = x0;
end

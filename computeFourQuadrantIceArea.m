function iceArea = computeFourQuadrantIceArea(md, icemask)
	% compute total ice area for the four quadrant: 1-NE, 2-NW, 3-SW, 4-SE
	% return iceArea, 4 rows, size(icemask,2) columns 
	iceArea = zeros(4, size(icemask,2));
	data = double(icemask<0);

	% quadrant 1
	pos = ((md.mesh.x>0) & (md.mesh.y>=0));
	quaddata = zeros(size(data));
	quaddata(pos,:) = data(pos,:);
	iceArea(1,:) =  integrateOverDomain(md, quaddata);

	% quadrant 2
	pos = ((md.mesh.x<=0) & (md.mesh.y>0));
	quaddata = zeros(size(data));
	quaddata(pos,:) = data(pos,:);
	iceArea(2,:) =  integrateOverDomain(md, quaddata);

	% quadrant 3
	pos = ((md.mesh.x<0) & (md.mesh.y<=0));
	quaddata = zeros(size(data));
	quaddata(pos,:) = data(pos,:);
	iceArea(3,:) =  integrateOverDomain(md, quaddata);

	% quadrant 4
	pos = ((md.mesh.x>=0) & (md.mesh.y<0));
	quaddata = zeros(size(data));
	quaddata(pos,:) = data(pos,:);
	iceArea(4,:) =  integrateOverDomain(md, quaddata);

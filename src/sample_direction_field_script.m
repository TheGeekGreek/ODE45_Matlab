y0 = 10;
tstart = 0;
tend = 10;
interval = 0:1:10;
direction_field(@f, [tstart, tend], y0, interval);

y0 = 10;
tstart = 0;
tend = 10;
interval = -1:0.1:1;
direction_field(@f, [tstart, tend], y0, interval);

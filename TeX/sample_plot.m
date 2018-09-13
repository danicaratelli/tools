F = {[(1:10)' rand(10,1) zeros(10,1)]};
filename = 'sample'
rs = 1;
cs = 1;
specs = struct();
specs.ylabs = {'y',''};
specs.xlabs = {'x','x'};
specs.linestyles = {'solid','dotted'};
specs.szs = [2 1];
specs.legs = {'rand','zero'};
specs.leg_pos = 1;
specs.titles = {'line'};
specs.cols = {'red','black'};
scale = 0.9;

write_texfig(F,filename,rs,cs,specs,scale,1,'\Large','')
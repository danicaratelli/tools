function write_texfig(F,filename,rs,cs,specs,scale,izero,fontsz,cap)
%This file allows to write a .tex file for every individual plot
%xlabs,ylabs,titles,legs,cols,szs,linestyles
num_plots = length(F);
assert(rs*cs>=num_plots,'not enough figure cells for the number of plots')

%for the plot
fileID = fopen([filename '.tex'],'w');
fprintf(fileID,'%s\n', '\begin{figure}[htpb!]');
fprintf(fileID,'%s\n','\begin{center}');
fprintf(fileID,'%s\n',['\begin{tikzpicture}[scale=' num2str(scale) ', transform shape]']);
fprintf(fileID,'%s\n',['\pgfplotsset{every tick label/.append style={font=' fontsz '}};']);
fprintf(fileID,'%s\n',['\pgfplotsset{y tick label style={  font=' fontsz ', /pgf/number format/precision=3,/pgf/number format/fixed}};']);
fprintf(fileID,'%s\n',['\pgfplotsset{y label style={  font=' fontsz '}};']);
fprintf(fileID,'%s\n',['\pgfplotsset{x label style={  font=' fontsz '}};']);
fprintf(fileID,'%s\n',['\begin{groupplot}[group style={group name=allgraphs, group size= ' num2str(cs) ' by ' num2str(cs) ', vertical sep=3cm, horizontal sep=2.5cm}, height = ' num2str(specs.height) 'cm, width = ' num2str(specs.width) 'cm];']);
%individual plots
for i=1:num_plots
    A = F{i}; %ith figure
    xmin = A(1,1); xmax = A(end,1);
    fprintf(fileID,'%s\n',['\pgfmathsetmacro{\xmin}{' num2str(xmin) '};']);
    fprintf(fileID,'%s\n',['\pgfmathsetmacro{\xmax}{' num2str(xmax) '};']);
    fprintf(fileID,'%s\n',['\nextgroupplot[ylabel={' specs.ylabs{i} '}, xlabel={' specs.xlabs{i} '},'...
         'xmin=\xmin, xmax=\xmax, tick label style={/pgf/number format/fixed},' ...
         'legend style={draw=none}, legend style={legend pos=outer north east}, '...         
         'legend columns=' num2str(size(A,2)) ',font=' fontsz '}];']);
     
     for j=2:size(A,2)
         strj = 'coordinates {';
         for k=1:size(A,1)
             strj = [strj '(' num2str(A(k,1)) ',' num2str(A(k,j)) ')'];
         end
         strj = [strj '}'];
         fprintf(fileID,'%s\n',['\addplot[' specs.linestyles{i,j-1} ', line width=' num2str(specs.szs(i,j-1)) 'pt, mark=none, ' specs.cols{i,j-1} '] ' strj ';']);
     end
     if ~isempty(specs.legs) && i == specs.leg_pos
         if izero == 1
             j_end = size(A,2)-1;
         else
             j_end = size(A,2);
         end
         for j=2:j_end
             fprintf(fileID,'%s\n',['\addlegendentry{' specs.legs{j-1} '\,\,};']);
         end
     end
end   
fprintf(fileID,'%s\n','\end{groupplot};');
%% Add titles
for i=1:num_plots
    [I,J] = ind2sub([rs cs],i);
    fprintf(fileID,'%s\n',['\node[above = 0.5cm of allgraphs c' num2str(J) 'r' num2str(I)... 
        ', font=' fontsz '] {' specs.titles{I,J} '};']);
end
fprintf(fileID,'%s\n','\end{tikzpicture}');
fprintf(fileID,'%s\n','\end{center}');
if ~isempty(cap)
    fprintf(fileID,'%s\n',['\caption{' cap '}']);
end
fprintf(fileID,'%s\n','\end{figure}');
fclose(fileID);

end
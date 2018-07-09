function [outputArg1,outputArg2] = write_texfig(F,filename,rs,cs,xlabs,ylabs,titles,legs,cols)
%This file allows to write a .tex file for every individual plot

num_plots = length(F);
assert(rs*cs>=num_plots,'not enough figure cells for the number of plots')

%for the plot
fileID = fopen([filename '.tex'],'w');
fprintf(fileID,'%s\n', '\begin{figure}[t]');
fprintf(fileID,'%s\n','\begin{center}');
fprintf(fileID,'%s\n','\begin{tikzpicture}[scale=0.7, transform shape]');
fprintf(fileID,'%s\n','\pgfplotsset{every tick label/.append style={font=\Large}};');
fprintf(fileID,'%s\n',['\begin{groupplot}[group style={group name=allgraphs, group size= ' num2str(cs) ' by ' num2str(rs) ', vertical sep=3cm, horizontal sep=2.5cm}];']);
%individual plots
for i=1:num_plots
    A = F{i}; %ith figure
    xmin = A(1,1); xmax = A(end,1);
    fprintf(fileID,'%s\n',['\pgfmathsetmacro{\xmin}{' num2str(xmin) '};']);
    fprintf(fileID,'%s\n',['\pgfmathsetmacro{\xmax}{' num2str(xmax) '};']);
    fprintf(fileID,'%s\n',['\nextgroupplot[ylabel={' ylabs{i} '}, xlabel={' xlabs{i} '}, label style={font=\Large},'...
         'xmin=\xmin, xmax=\xmax, tick label style={/pgf/number format/fixed},' ...
         'legend style={draw=none}, legend style={at={(1.15,-0.5)}, anchor=south,'...
         'legend columns=' num2str(size(A,2)) ',font=\Large}];']);
     for j=2:size(A,2)
         strj = 'coordinates {';
         for k=1:size(A,1)
             strj = [strj '(' num2str(A(k,1)) ',' num2str(A(k,j)) ')'];
         end
         strj = [strj '}'];
         fprintf(fileID,'%s\n',['\addplot[line width=2pt, mark=none, ' cols{j-1} '] ' strj ';']);
     end
     for j=2:size(A,2)
         fprintf(fileID,'%s\n',['\addlegendentry{' legs{j-1} '\,\,};']);
     end
end   
fprintf(fileID,'%s\n','\end{groupplot};');
%% Add titles
for i=1:num_plots
    [I,J] = ind2sub([rs cs],i);
    fprintf(fileID,'%s\n',['\node[above = 0.5cm of allgraphs c' num2str(J) 'r' num2str(I)... 
        ', font=\Large] {' titles{i} '};']);
end
fprintf(fileID,'%s\n','\end{tikzpicture}');
fprintf(fileID,'%s\n','\end{center}');
fprintf(fileID,'%s\n','\end{figure}');
fclose(fileID);

end
"""
    write_texbarplot(F,filename,rs,cs,specs,scale,izero,fontsz,cap)
    Convert plot to pretty .tex file.
   Inputs:
   -------
   ``F``        :       Array,  array of array containing each line separately\n
   ``filename`` :       String, name of file\n
   ``rs``       :       Int64,  number of rows in plot (sublots)\n
   ``cs``       :       Int64,  number of columns in plot (sublots)\n
   ``specs``    :       Dict,   plot specifications\n
   ``izero``    :       Bool,   boolean indicating whether to include 0 line\n
   ``fontsz``   :       Int64,  font size\n
   ``cap``      :       String, caption\n
   -------
   ```
"""
function write_texbarplot(F,filename,rs,cs,specs,scale,fontsz,cap)
    #This file allows to write a .tex file for every individual plot
    #xlabs,ylabs,titles,legs,cols,szs,linestyles
    num_plots = length(F);
    @assert(rs*cs>=num_plots,"not enough figure cells for the number of plots")

    #for the plot
    open(filename*".tex","w") do fileID
       write(fileID,"\\begin{figure}[htpb!]\n");
       write(fileID,"\\begin{center}\n");
       write(fileID,"\\textbf{\\fontsize{"*string(specs["title_font"])*"}{0} \\selectfont "*specs["title"]*"}\\par\\medskip\n")
       write(fileID,"\\begin{tikzpicture}[scale="*string(scale)*", transform shape]\n");
       #individual plots
       for i=1:num_plots
           A = F[i]; #ith figure
           write(fileID,"\\begin{axis} [ybar, xtick pos=left, ytick pos=left,"*
                        "height="*string(specs["height"])*
                        "cm, width = "*string(specs["width"])*"cm, xlabel = "*
                        specs["xlabs"][i]*", ylabel = "*specs["ylabs"][i]
                        *", legend style = "*specs["legend style"]*"]\n");
            for j=2:size(A,2)
                strj = "coordinates {";
                for k=1:size(A,1)
                    xcoor = string(A[k,1]);
                    strj = strj*"("*xcoor*","*string(A[k,j])*")";
                end
                strj = strj*"}";
                if isequal(specs["legs"][i,j-1],"")
                    write(fileID,"\\addplot[line width = "*specs["line width"][i,j-1]*", bar width = "*specs["bar width"][i,j-1]*
                                  ", fill = "*specs["fill"][i,j-1]*", bar shift = "*specs["bar shift"][i,j-1]*
                                  ", area legend, fill opacity = "*string(specs["fill opacity"][i,j-1])*", forget plot] "*strj*";\n");
                else
                    write(fileID,"\\addplot[line width = "*specs["line width"][i,j-1]*", bar width = "*specs["bar width"][i,j-1]*
                                  ", fill = "*specs["fill"][i,j-1]*", bar shift = "*specs["bar shift"][i,j-1]*
                                  ", area legend, fill opacity = "*string(specs["fill opacity"][i,j-1])*"] "*strj*";\n");
                end
            end
            lgnd = join(specs["legs"],", ");
            if !isempty(lgnd)
                write(fileID,"\\legend{"*lgnd*"};\n");
            end
       end

       write(fileID,"\\end{axis};\n");
       write(fileID,"\\end{tikzpicture}\n");
       write(fileID,"\\end{center}\n");
       if !isnothing(cap)
           write(fileID,"\\caption{"*cap*"}\n");
       else
           write(fileID,"\\caption[Entry for the List of Figures (LoF)]{}\\label{fig:without}");
       end
       write(fileID,"\\end{figure}");
    end
end

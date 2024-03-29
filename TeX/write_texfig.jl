"""
    write_texfig(F,filename,rs,cs,specs,scale,izero,fontsz,cap)
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
function write_texfig(F,filename,rs,cs,specs,scale,izero,ivert,fontsz,cap)
    #This file allows to write a .tex file for every individual plot
    #xlabs,ylabs,titles,legs,cols,szs,linestyles
    num_plots = length(F);
    @assert(rs*cs>=num_plots,"not enough figure cells for the number of plots")

    #for the plot
    open(filename*".tex","w") do fileID
       write(fileID,"\\begin{figure}[htpb!]\n");
       write(fileID,"\\begin{center}\n");
       write(fileID,"\\textbf{\\fontsize{"*string(specs["title_font"])*"}{0} \\selectfont "*specs["titles"]*"}\\par\\medskip\n")
       write(fileID,"\\begin{tikzpicture}[scale="*string(scale)*", transform shape]\n");
       write(fileID,"\\pgfplotsset{every tick label/.append style={font=\\fontsize{"*string(fontsz)*"}{0}\\selectfont"*"}};\n");
       write(fileID,"\\pgfplotsset{y tick label style={  font=\\fontsize{"*string(fontsz)*"}{0}\\selectfont"*", /pgf/number format/precision=3,/pgf/number format/fixed}};\n");
       write(fileID,"\\pgfplotsset{y label style={  font=\\fontsize{"*string(fontsz)*"}{0}\\selectfont"*"}};\n");
       write(fileID,"\\pgfplotsset{x label style={  font=\\fontsize{"*string(fontsz)*"}{0}\\selectfont"*"}};\n");
       write(fileID,"\\begin{groupplot}[xtick pos=left, ytick pos=left, group style={group name=allgraphs, group size= "*string(cs)*" by "*string(cs)*", vertical sep=3cm, horizontal sep=2.5cm}, height = "*string(specs["height"])*"cm, width = "*string(specs["width"])*"cm];\n");
       #individual plots
       for i=1:num_plots
           A = F[i]; #ith figure
           if isnothing(specs["xmin"]); xmin = A[1,1]; else; xmin = specs["xmin"]; end;
           if isnothing(specs["xmax"]); xmax = A[end,1]; else; xmax = specs["xmax"]; end;
           write(fileID,"\\pgfmathsetmacro{\\xmin}{"*string(xmin)*"};\n");
           write(fileID,"\\pgfmathsetmacro{\\xmax}{"*string(xmax)*"};\n");
           write(fileID,"\\nextgroupplot[")
           if haskey(specs,"xtick")
                write(fileID,"unbounded coords=jump, xtick={" * specs["xtick"] * "}, xticklabels={" * specs["xticklabels"] * "},");
           end
           write(fileID,"ylabel={"*specs["ylabs"][i]*"}, xlabel={"*specs["xlabs"][i]*"},
                xmin=\\xmin, xmax=\\xmax, tick label style={/pgf/number format/fixed},
                legend style={draw=none}, legend style={legend pos=south east},
                legend columns=1,font=\\fontsize{"*string(fontsz)*"}{0}\\selectfont"*"];\n");
            for j=2:size(A,2)
                strj = "coordinates {";
                for k=1:size(A,1)
                    if ivert && j==size(A,2)-1
                        xcoor = "0";
                    else
                        xcoor = string(A[k,1]);
                    end
                    strj = strj*"("*xcoor*","*string(A[k,j])*")";
                end
                strj = strj*"}";
                if  isequal(specs["legs"][j-1],"")
                    write(fileID,"\\addplot["*specs["linestyles"][i,j-1]*", line width="*string(specs["szs"][i,j-1])*"pt, mark=none, "*specs["cols"][i,j-1]*", name path="*specs["path_names"][j-1]*",forget plot] "*strj*";\n");
                else
                    write(fileID,"\\addplot["*specs["linestyles"][i,j-1]*", line width="*string(specs["szs"][i,j-1])*"pt, mark=none, "*specs["cols"][i,j-1]*", name path="*specs["path_names"][j-1]*"] "*strj*";\n");
                end
                #including confidence intervals
                if !isempty(specs["shades"][j-1])
                    write(fileID,"\\tikzfillbetween[\n of="*specs["shades"][j-1]*"\n] {color="*specs["shade_cols"][j-1]*",semitransparent};\n");
                end
            end
            if !isnothing(specs["legs"]) && i == specs["leg_pos"]
                if izero == 1
                    j_end = size(A,2)-1;
                else
                    j_end = size(A,2);
                end
                for j=2:j_end
                    if !isequal(specs["legs"][j-1],"")
                        write(fileID,"\\addlegendentry{"*specs["legs"][j-1]*"\\,\\,};\n");
                    end
                end
            end
       end

       write(fileID,"\\end{groupplot};\n");

       # Add titles
       #=
       for i=1:num_plots
           I,J = Base._ind2sub([rs cs],i);
           write(fileID,"\\node[above = 0.5cm of allgraphs c"*string(J)*"r"*string(I)*",
                         font=\\fontsize{"*string(fontsz)*"}{0}\\selectfont"*"] {"*specs["titles"][I,J]*"};\n");
       end
       =#
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

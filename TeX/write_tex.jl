"""
   ...
   # Description
   writes a Julia DataFrame table to a .tex table in LaTeX

   # Arguments
   - `filename::String`: name of .tex file to which table is written.
   - `tab::DataFrames.DataFrame`: table to be written to file.
   - `caption::String`: (optional) caption for table.
   - `style::String`: (optional) table style, e.g. "c|ccc".
   - `num_prec::Int64`: (optional) numerical precision of numbers in table.
   - `linesep::String`: (optional) add "" if you do not want rows to be separated.

   # Output
   - .tex table will be written to file.
   ...
"""



function write_tex(filename, tab, caption="", style="", num_prec=2,linesep="\\hline")
    nrows = size(tab,1); ncols = size(tab,2);
    col_names = names(tab);
    col_names = map(x->string(col_names[x]),1:ncols);
    if style==""
        style = "@{}"*repeat("l",ncols)*"@{}";
    end
    open(filename*".tex","w") do f
        write(f,"\\begin{table}[htpb!]\r\n");
        write(f,"\\centering\r\n");
        write(f,"\\renewcommand{\\arraystretch}{1.2}\r\n")
        write(f,"\\begin{tabular}{"*style*"}\r\n");

        #header
        header = "";
        for j=1:ncols
            if j==1
                header = col_names[j];
            else
                header = header*" & "*col_names[j];
            end
        end
        header = header*"\\"*"\\"*"\\hline \\hline \r\n";

        tabentry = "";
        for i=1:nrows
            for j=1:ncols
                if typeof(tab[i,j])<:Number
                    tmp_entry = string(round(tab[i,j],num_prec));
                else
                    tmp_entry = tab[i,j];
                end
                if j==1
                    rowi = tmp_entry
                else
                    rowi = rowi*" & "*tmp_entry
                end
            end
            rowi = rowi*"\\"*"\\"*linesep*" \r\n";
            tabentry = tabentry*" "*rowi;
        end
        tablebod = header*tabentry;
        write(f,tablebod);
        write(f,"\\end{tabular}\r\n");
        write(f,"\\caption{"*caption*"} \r\n");
        write(f,"\\label{"*filename*"} \r\n");
        write(f,"\\end{table}\r\n");
    end
    close(filename*".tex")
end

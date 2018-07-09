for write_texfig.m

the .tex file must include a preamble, for example, the following:

\documentclass{article}
\usepackage{tikz}
\usepackage[graphics, active, tightpage]{preview}
\PreviewEnvironment{tikzpicture}

%!tikz preamble begin
\usepackage{tikz}
\usetikzlibrary{shapes.geometric}
\usetikzlibrary{shapes.arrows}
\usetikzlibrary{calc}
\usepackage{array}
\usetikzlibrary{arrows}
\usetikzlibrary{positioning}
\usetikzlibrary{intersections}
\usetikzlibrary{matrix}

\usepackage{pgfplots} 
\usepgfplotslibrary{dateplot}
\usepgfplotslibrary{groupplots}
\pgfplotsset{compat=1.3}

\usepackage{color}
\definecolor{webgreen}{rgb}{0,.5,0}
\definecolor{webbrown}{rgb}{.6,0,0}
\definecolor{webpurple}{rgb}{0.7,0,0.7}
\definecolor{blue}{RGB}{0,114,178}
\definecolor{red2}{RGB}{213,94,0}
\definecolor{green2}{RGB}{0,158,115}

%!tikz preamble end


\begin{document}

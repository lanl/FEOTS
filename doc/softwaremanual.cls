\LoadClass[12pt]{memoir}
% Package for doing custom title-page
\usepackage{titling}
%
\setcounter{tocdepth}{4}
\setcounter{secnumdepth}{4}

\usepackage[margin=1.0in]{geometry}
\renewcommand{\baselinestretch}{1.5}
% AMS and mathematics packages
\usepackage{amsmath,amsthm}
\usepackage{amsfonts}
\usepackage{amssymb}
\usepackage[mathscr]{euscript}

% Figures and controlling packages
\usepackage{float}
\usepackage{wrapfig}
\usepackage{graphicx}
\usepackage{xcolor}
%  ---- Special command for setting up logo for the title page
%\usepackage{adjustbox}
%\newcommand{\logo}{
%  \adjustbox{valign=t}{\includegraphics[width=0.25\textwidth]{images/shallowWater.png}}
%}

\usepackage{pgfplots}
%\usetikzlibrary{pgfplots.groupplots}
\usepackage{helvet}
\usepackage[eulergreek]{sansmath}
%\pgfplotsset{compat=newest,
%  tick label style = {font=\sansmath\sffamily},
%  every axis label = {font=\sansmath\sffamily},
%  legend style = {font=\sansmath\sffamily},
%  label style = {font=\sansmath\sffamily}
%}
\usepackage{verbatim}

%Bibliography and Appendices
\usepackage[round]{natbib}
\usepackage{appendix}

\usepackage{tabularx}
\usepackage{blindtext,pgffor}
%\usepackage{tcolorbox}
%\tcbuselibrary{breakable}
%\usepackage{tcolorbox}
% new tcolorbox environment
% #1: tcolorbox options
% #2: color
% #3: box title
%\newtcolorbox{mybox}[3][]
%{ breakable,
%  colframe = #2!50,
%  colback  = #2!10,
%  coltitle = #2!20!black,  
%  title    = #3,
%  #1,
%}

% Allow in-document and external link referencing with hyperref
\definecolor{webred}{rgb}{0.5, 0, 0} % less intense red
\usepackage{hyperref}
\hypersetup{ ps2pdf,
    bookmarks=true,         % show bookmarks bar?
    bookmarksnumbered=false,
    unicode=false,          % non-Latin characters in Acrobat’s bookmarks
    pdftoolbar=true,        % show Acrobat’s toolbar?
    pdfmenubar=true,        % show Acrobat’s menu?
    pdffitwindow=false,     % window fit to page when opened
    pdfstartview={FitH},    % fits the width of the page to the window
    pdftitle={FEOTS User Guide},    % title
    pdfauthor={Joseph Schoonover},     % author
    pdfnewwindow=true,       % links in new PDF window
    colorlinks=true,         % false: boxed links; true: colored links
    linkcolor=webred,        % color of internal links (change box color with linkbordercolor)
    citecolor=purple,        % color of links to bibliography
    filecolor=magenta,       % color of file links
    urlcolor=blue,           % color of external links
    linktoc=all,
}



% Adjusting the section,chapter, etc. headings
\usepackage{titlesec}
\newcommand*{\justifyheading}{\raggedleft}

\titleformat{\chapter}{\normalfont\sffamily\huge\bfseries\justifyheading\color{black}}{\thechapter }{20pt}{\Huge}
  
  
\titleformat{\section}{\normalfont\sffamily\Large\bfseries\color{black}}{\thesection}{1em}{}

\titleformat{\subsection}{\normalfont\sffamily\large\bfseries\color{black}}{\thesection}{1em}{}

\titleformat{\subsubsection}{\normalfont\sffamily\large\bfseries\color{black}}{\thesection}{1em}{}

% Header and Footer

\makeevenhead{myheadings}{}{}{}
\makeoddhead{myheadings}{}{}{}

\makeevenfoot{myheadings}{\thepage}{}{}

\makeoddfoot{myheadings}{}{}{\thepage}

\makefootrule{myheadings}{\textwidth}{\normalrulethickness}{0ex}


\chapterstyle{mychapter}
%  Adjust the chapter page footer
\copypagestyle{chapter}{plain}
\makeoddfoot{chapter}{}{}{\thepage}
\makeevenfoot{chapter}{\thepage}{}{}

% Add a line above the footer
\makefootrule{chapter}{\textwidth}{\normalrulethickness}{0ex}
% Chapter Style

\raggedbottom



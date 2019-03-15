UQ LaTeX thesis template
Copyright (C) 2018 James Bennett

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

Maintained by james.bennett2 (at) uqconnect.edu.au
Please contact me if there are any bugs.
Latest updated: 04 Jun 2018

**********************************************
Developed by James Bennett.
This contains sections of code provided by Paul Cochrane (2004) [reproduced and modified under the GNU General Public License v2] and Will Petterrsen (2015).
Thanks to Christiaan Bekker, Thomas Bell, Catxere Andrade Casacio, Nicolas Mauranyapin, Varun Prakash, and Alexander Stilgoe for assistance.

**********************************************
Please follow the instructions given below and those given in the template itself. This will produce a thesis with formatting that conforms to UQ's guidelines.

**********************************************
Designed for use with pdfLaTeX. *Requires* packages:
- uqthesis class (included)
- uqtitle (included)
- memoir class
- graphicx
- mathptmx
- booktabs
- longtable
- url
- hyperref
- geometry
- memhfixc
- amsmath, amssymb, amsthm
- setspace
- comment

Will also use packages (style.tex):
- cite
- rotating
- pdfpages
- wrapfig
- bm
- upgreek
- dsfont
- mathtools
- framed
- microtype
- marvosym
- graphic
- color
- transparent
- placeins
- mdframed, mdwlist

**********************************************
This template consists of a .zip file containing the following:
- masterfile.tex is the main file which receives the others as inputs.
	- front.tex contains/produces the front matter (everything prior to the body of the document)
	- back.tex contains the back matter (optional)
	- symbols.tex provides space for the list of abbreviations and (optional) list of symbols
The four directories (Abstract, Chap2, Introduction, and Conclusion) contain source files for each section, as indicated. These make up the body of the thesis. The files have been organised into separate folders so that figure files can be stored separately for each chapter (reduces clutter). Edit and duplicate as necessary.
Note that abstract.tex contains instructions for compiling a stand-alone version of the abstract (for submission purposes).

This template also contains
- uqthesis.cls : loads necessary packages, sets the page and heading styles, defines the linespread, and redefines the \sqrt{} command**.
- uqtitle.cls : sets the title page style and defines storage commands (\currentdegrees, \submittedfor, etc.)
- style.tex : defines 'optional' packages. It also defines some optional text macros. Some elements of this file may be altered if required.
- copyrightstatement.tex : this statement must appear verbatim. Do not alter!
- logo.pdf : the UQ shield 'logo' is provided in pdf and eps formats (default is pdf - you will need to alter the class file if you want to use EPS).
- abbrv-unsrt.bst : numerical citation style by E.Krepska. This provides short numerical citations in order of citation; my personal recommendation.
- GNU_License.txt : copy of GNU General Public License v3

To use:
- compile masterfile.tex using pdflatex and ensure that you have all of the necessary packages. Add or subtract to the packages used as required for your thesis.
- use the class option 'final' to suppress instructional text in the output
- insert the appropriate bibstyle and .bib files into masterfile.tex to allow BibTeX to function correctly
- add each chapter as a separate directory and source .tex file, along with all figures. Add \input{} commands for each in masterfile.tex.
- carefully read the instructional text in front.tex and back.tex and complete these sections. Failure to follow these instructions will likely result in your thesis being rejected.
- complete your list of abbreviations (and symbols) in symbols.tex (there are TeX packages which assist in this; add them if you wish).

Enjoy!

** I've redefined the default length of whitespace following the \sqrt command in order to improve readability (AMS math tends to cram everything very close to the end of the surd). Just comment out the code in the class file if this bothers you.
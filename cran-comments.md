*Test Environments*

-   locally on Apple M1 Pro R version 4.2.1
-   Windows Server 2022, R-devel, 64 bit
-   Ubuntu Linux 20.04.1 LTS, R-devel, GCC

*R CMD check results*

0 errors \| 0 warnings \| 3 notes

*Notes*

-   Windows Server 2022, R-devel, 64 bit
    -   checking Rd cross-references ... NOTE Package unavailable to check Rd xrefs: 'readxl' *[package unavailable on test machine]*
    -   checking for detritus in the temp directory ... NOTE Found the following files/directories: 'lastMiKTeXException' *[all checks run on this rhub environment seem to give this message at the moment, there is no local file or dir with this name]*
-   Ubuntu Linux 20.04.1 LTS, R-devel, GCC
    -   checking Rd cross-references ... NOTE Package unavailable to check Rd xrefs: 'readxl' *[package unavailable on test machine]*
    -   checking HTML version of manual ... NOTE Skipping checking HTML validation: no command 'tidy' found Skipping checking math rendering: package 'V8' unavailable *[Seems to be related to the setup on the test machine rather than mrds]*

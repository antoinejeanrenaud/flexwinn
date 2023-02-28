## Test environments
- R-hub windows-x86_64-devel (r-devel)
- R-hub ubuntu-gcc-release (r-release)
- R-hub fedora-clang-devel (r-devel)

## R CMD check results
❯ On windows-x86_64-devel (r-devel), ubuntu-gcc-release (r-release)
  checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Antoine Jeanrenaud <antoine2104@hotmail.com>'
  
  New submission
  
  Version contains large components (0.0.0.9000)

❯ On windows-x86_64-devel (r-devel)
  checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'

❯ On fedora-clang-devel (r-devel)
  checking CRAN incoming feasibility ... [6s/27s] NOTE
  Maintainer: ‘Antoine Jeanrenaud <antoine2104@hotmail.com>’
  
  New submission
  
  Version contains large components (0.0.0.9000)

❯ On fedora-clang-devel (r-devel)
  checking HTML version of manual ... NOTE
  Skipping checking HTML validation: no command 'tidy' found

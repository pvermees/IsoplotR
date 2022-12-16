# IsoplotR

**IsoplotR** is a free and open-source substitute for Kenneth Ludwig's
popular **Isoplot** add-in to Microsoft Excel.  This Github repository
contains the code for IsoplotR's core commands. A graphical user
interface is provided in a separate repository at
[https://github.com/pvermees/IsoplotRgui](https://github.com/pvermees/IsoplotRgui).

## Installation

You must have **R** installed on your system (see
[https://www.r-project.org/](https://www.r-project.org/)). The most
recent stable version of IsoplotR is available from **CRAN** at
[https://cran.r-project.org/package=IsoplotR](https://cran.r-project.org/package=IsoplotR)
and can be installed on your system as follows:

```
install.packages('IsoplotR')
```

Alternatively, the current development version of IsoplotR can be installed from Github with the **remotes** package:

```
install.packages('remotes')
remotes::install_github('pvermees/IsoplotR')
```

## Example

Once installed, IsoplotR can be loaded into memory by entering the following code at the R command prompt:

```
library(IsoplotR)
```

Now we can issue IsoplotR commands to R. For example:

```
setwd(system.file(package='IsoplotR')) # navigate to the built-in data files
RbSr <- read.data('RbSr1.csv',method='Rb-Sr',format=1)  
isochron(RbSr)
```

## Further information

See [https://isoplotr.london-geochron.com](https://www.ucl.ac.uk/~ucfbpve/isoplotr/)

[Vermeesch, P., 2018, IsoplotR: a free and open toolbox for
geochronology. Geoscience Frontiers, v.9, p.1479-1493, doi:
10.1016/j.gsf.2018.04.001.](https://www.ucl.ac.uk/~ucfbpve/papers/VermeeschGSF2018/)

## Author

[Pieter Vermeesch](https://www.ucl.ac.uk/~ucfbpve/)

## License

This project is licensed under the GPL-3 License

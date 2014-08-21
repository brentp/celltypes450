Adjust for cell composition in 450K data.

Much of the adjust.beta code was originally written by Andres Houseman.
I have only put it into a simple to use function and packaged it.

Installation
============

    > library(devtools)
    > install_github("brentp/celltypes450")

About
=====

See [this thread](https://groups.google.com/forum/#!searchin/epigenomicsforum/brent/epigenomicsforum/IvON6p_hmz4/SZU9hZXsj9wJ) for details.

Compared to minfi, which normalizes the sorted data with the new data:

![comparison](https://17886685092018163118.googlegroups.com/attach/dc8fec95853d9549/ct.png?part=0.1&view=1&vt=ANaJVrFxKvcTQAoYohUVDC5okdDxOR-z5H-Q2p9vSba3kR8pmfNT7NFogl_zkYZraf7iTeDw-nMIgU1CYedBB056qcvApwexhemLEwTmBNyQeFxvwXqm-XY)


To get your beta (0, 1) values adjusted for cell-type composition, use:

```R

library(devtools)
install_github("brentp/celltypes450")
library(celltypes450)
adjusted = adjust.beta(beta)
```

To get the cell-composition estimates such as those from the image above:

```R
cell.types = adjust.beta(beta, est.only=TRUE)
```

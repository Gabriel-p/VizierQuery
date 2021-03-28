
# VizierQuery

Small wrapper around the [`astroquery`](https://astroquery.readthedocs.io) package to download Vizier data for any available catalog.

### Requirements

    Python 3.8, numpy, astropy, astroquery, uncertainties

Can be used in a `conda` environment with:

    $ conda create -n vizierq numpy astropy
    $ conda activate vizierq
    $ conda install -c conda-forge uncertainties
    $ conda install -c astropy astroquery

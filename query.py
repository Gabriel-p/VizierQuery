
import configparser
import numpy as np
import astropy.coordinates as coord
import astropy.units as u
from astropy.io import ascii
from astroquery.vizier import Vizier
from uncertainties import ufloat
from uncertainties import unumpy as unp


__version__ = '1.0'


def main():
    """
    Explore data downloaded via the 'astroquery' package.
    """
    print("\n*******************")
    print(" VizierQuery v{}".format(__version__))
    print("*******************")

    cat, mag_name, mag_max, colors, clusters = readInput()

    # Gaia DR2
    if cat == 'I/345/gaia2':
        DR = '2'
    # Gaia EDR3
    elif cat == 'I/350/gaiaedr3':
        DR = '3'
    else:
        DR = None

    for clust in clusters:
        center, box_s = (clust['cent_ra'], clust['cent_dec']), clust['box_s']

        data = getData(
            cat, mag_name, mag_max, clust['name'], center, box_s)

        N_old = len(data)
        print("{} data read, {} sources".format(clust['name'], N_old))

        if colors:
            print("Obtaining colors and their uncertainties")
            data = uncertMags(DR, data, colors)

        print("Write output file")
        ascii.write(data, clust['name'] + ".dat", overwrite=True)

    print("\nEnd")


def getData(cat, mag_name, mag_max, name, center, box_s):
    """
    Download data using astroquery.
    """
    print("\nDownloading data for {}, from '{}'".format(name, cat))
    print("  {} < {}".format(mag_name, mag_max))
    # Unlimited rows, all columns
    v = Vizier(row_limit=-1, columns=['all'])

    result = v.query_region(coord.SkyCoord(
        ra=center[0], dec=center[1], unit=(u.deg, u.deg), frame='icrs'),
        width=box_s, catalog=[cat],
        column_filters={mag_name: '<{}'.format(mag_max)})
    data = result[cat]

    return data


def uncertMags(DR, data, colors):
    """
    # Gaia DR2 zero points:

    https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/
    chap_cu5pho/sec_cu5pho_calibr/ssec_cu5pho_calibr_extern.html#Ch5.T2

    The G magnitude error:
        unp.std_devs(mag_d['G'])
    is equivalent to:
        np.sqrt((1.0857*(data['e_FG']/data['FG']))**2 + .0018**2)
    as defined in Eq 5.26 of the link above.

    These values are larger than the 'e_Gmag' column given by Vizier by up to
    ~0.002 for the brightest stars, and go to zero for the faintest. I don't
    really know why. I asked Vizier and their answer was:

    > The e_Gmag uncertainties in VizieR were added by the CDS and,
    > apparently, do not take into account the Vegamag corrections; we will
    > continue the investigation.
    >
    > In the meantime, ignoring the e_Gmag in VizieR and using the formula
    > given by Gaia DR2 seems to be the right solution."

    # Gaia EDR3 zero points:

    https://www.cosmos.esa.int/web/gaia/edr3-passbands

    """
    # Zero points for the G,BP,RP magnitudes.
    if DR == '2':
        # Updated October 2017
        Zp_G = ufloat(25.6914396869, 0.0011309370)
        Zp_BP = ufloat(25.3488107670, 0.0004899854)
        Zp_RP = ufloat(24.7626744847, 0.0035071711)
    elif DR == '3':
        Zp_G = ufloat(25.6873668671, 0.0027553202)
        Zp_BP = ufloat(25.3385422158, 0.0027901700)
        Zp_RP = ufloat(24.7478955012, 0.0037793818)

    if DR in ('2', '3'):
        # Fluxes
        I_G = unp.uarray(data['FG'], data['e_FG'])
        I_BP = unp.uarray(data['FBP'], data['e_FBP'])
        I_RP = unp.uarray(data['FRP'], data['e_FRP'])
        # Magnitudes
        mag_d = {
            'Gmag': Zp_G + -2.5 * unp.log10(I_G),
            'BPmag': Zp_BP + -2.5 * unp.log10(I_BP),
            'RPmag': Zp_RP + -2.5 * unp.log10(I_RP)}

    for col in colors:
        if DR in ('2', '3'):
            col_v = mag_d[col[0]] - mag_d[col[1]]
            e_col = unp.std_devs(col_v)
            col_v = np.array(unp.nominal_values(col_v))
        else:
            col_v = data[col[0]] - data[col[1]]
            e_col = np.sqrt(data['e_' + col[0]]**2 + data['e_' + col[1]]**2)

        col_n = col[0] + '-' + col[1]
        # Add color
        try:
            data[col_n]
            print("Column '{}' already exists".format(col_n))
        except KeyError:
            data[col_n] = col_v

        # Add errors for the color
        try:
            data['e_' + col_n]
            print("Column '{}' already exists".format('e_' + col_n))
        except KeyError:
            data['e_' + col_n] = e_col

    return data


def readInput():
    """
    Read 'params.ini' data file.
    """
    in_params = configparser.ConfigParser()
    in_params.read('params.ini')

    pars = in_params['Parameters']
    cat, mag_name, mag_max = pars['cat'], pars['mag_name'],\
        pars.getfloat('mag_max')

    cols_data = in_params.items("Colors")
    colors = []
    for col in cols_data:
        m1, m2 = col[1].split()
        colors.append([m1, m2])

    clust_data = in_params.items("Clusters")
    clusters = []
    for clust in clust_data:
        ra, dec, box_s = clust[1].split()
        clusters.append({
            'name': clust[0], 'cent_ra': float(ra),
            'cent_dec': float(dec), 'box_s': box_s})

    return cat, mag_name, mag_max, columns, colors, clusters


if __name__ == '__main__':
    # To see available catalogs:
    # catalog_list = Vizier.find_catalogs('Pan-STARRS')
    # catalogs = Vizier.get_catalogs(catalog_list.keys())
    # print(catalogs)

    # Generate output dir if it doesn't exist.
    if not exists('out'):
        makedirs('out')
    main()

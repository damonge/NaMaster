#!/bin/bash

if [ ! -f lambda_sfd_ebv.fits ] ; then
    echo "Downloading reddening map"
    wget https://lambda.gsfc.nasa.gov/data/foregrounds/SFD/lambda_sfd_ebv.fits
fi

if [ ! -f mask_lss.fits ] ; then
    echo "Missing mask!"
    exit
fi

if [ ! -f star_template.fits ] ; then
    echo "Missing star template!"
    exit
fi

if [ ! -f nz.txt ] ; then
    echo "Missing N(z)"
    exit
fi

echo "Generating LSS power spectra"
python get_lss_cls.py --plot

echo "Generating LSS contaminant templates"
python get_lss_contaminants.py --plot --nside 512

echo "Generating CMB mask"
python get_cmb_mask.py --plot --nside 256

echo "Generating CMB contaminant templates"
python get_cmb_contaminants.py --plot --nside 256

#!/usr/bin/env python
"""
Create test FITS data for testing the pipeline
"""

import numpy as np
from astropy.io import fits
from pathlib import Path

def create_test_fits(filename, shape=(100, 512), signal_pos=50):
    """
    Create a simple test 2D spectrum FITS file
    
    Parameters
    ----------
    filename : str
        Output filename
    shape : tuple
        (spatial, wavelength) shape
    signal_pos : int
        Position of signal in spatial direction
    """
    np.random.seed(42)
    
    # Create synthetic 2D spectrum
    # Gaussian in spatial direction, with some wavelength features
    spatial_axis = np.arange(shape[0])
    wavelength_axis = np.arange(shape[1])
    
    # Create image with Gaussian profile
    spatial_profile = np.exp(-0.5 * ((spatial_axis - signal_pos) / 5.0) ** 2)
    
    # Add some spectral features
    spectrum_profile = np.ones(shape[1])
    spectrum_profile[100:110] *= 1.5  # Feature 1
    spectrum_profile[250:270] *= 1.8  # Feature 2
    spectrum_profile[400:420] *= 1.3  # Feature 3
    
    # 2D image
    image = np.outer(spatial_profile, spectrum_profile)
    
    # Add noise
    image += np.random.normal(0, 0.1, shape)
    
    # Convert to uint16
    image = (image * 1000).astype(np.uint16)
    
    # Create FITS
    hdu = fits.PrimaryHDU(data=image)
    hdu.header['OBJECT'] = 'TEST'
    hdu.header['INSTRUME'] = 'KOspec'
    hdu.header['NAXIS1'] = shape[1]
    hdu.header['NAXIS2'] = shape[0]
    
    hdu.writeto(filename, overwrite=True)
    print(f"Created {filename}")


if __name__ == '__main__':
    # Create spectra directory
    spectra_dir = Path('spectra')
    spectra_dir.mkdir(exist_ok=True)
    
    # Create test object pairs
    objects = ['star001', 'star002']
    
    for obj in objects:
        # Create A-frame (signal at center)
        create_test_fits(f'spectra/{obj}_A.fits', shape=(100, 512), signal_pos=50)
        
        # Create B-frame (signal shifted or different)
        create_test_fits(f'spectra/{obj}_B.fits', shape=(100, 512), signal_pos=52)
    
    print("Test data created successfully!")
    print("Run: python main.py")

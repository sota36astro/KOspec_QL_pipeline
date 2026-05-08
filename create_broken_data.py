#!/usr/bin/env python
"""
Create broken FITS files for testing error handling
"""

from pathlib import Path

def create_broken_fits():
    """Create various broken FITS files"""
    spectra_dir = Path('spectra')
    
    # Test 1: Truncated file
    with open(spectra_dir / 'broken001_A.fits', 'wb') as f:
        f.write(b'SIMPLE  =                    T / file does conform to FITS standard')
        # Truncate without proper ending
    
    # Test 2: Non-FITS file
    with open(spectra_dir / 'broken001_B.fits', 'w') as f:
        f.write('This is not a FITS file\n')
    
    # Test 3: Empty file
    Path(spectra_dir / 'broken002_A.fits').touch()
    Path(spectra_dir / 'broken002_B.fits').touch()
    
    print("Created broken test files for error handling tests")

if __name__ == '__main__':
    create_broken_fits()

"""
FITS file loader module
"""

import numpy as np
from pathlib import Path
import logging
from astropy.io import fits

logger = logging.getLogger(__name__)


class FITSLoader:
    """Handle FITS file loading with error handling"""
    
    def __init__(self):
        pass
    
    @staticmethod
    def load(filepath):
        """
        Load FITS file safely
        
        Parameters
        ----------
        filepath : str or Path
            Path to FITS file
        
        Returns
        -------
        data : 2D array or None
            Image data, or None if loading failed
        header : dict or None
            FITS header, or None if loading failed
        success : bool
            Whether loading was successful
        error_msg : str
            Error message if unsuccessful
        """
        filepath = Path(filepath)
        
        if not filepath.exists():
            return None, None, False, f"File not found: {filepath}"
        
        try:
            with fits.open(filepath) as hdul:
                if len(hdul) == 0:
                    return None, None, False, f"Empty FITS file: {filepath}"
                
                # Try to find primary HDU with data
                data = None
                header = None
                
                for hdu in hdul:
                    if isinstance(hdu, (fits.PrimaryHDU, fits.ImageHDU)):
                        if hdu.data is not None:
                            if data is None:  # Take first valid HDU
                                data = hdu.data
                                header = dict(hdu.header)
                
                if data is None:
                    return None, None, False, f"No image data in {filepath}"
                
                # Verify data type
                if not isinstance(data, np.ndarray):
                    return None, None, False, f"Data is not array in {filepath}"
                
                if data.ndim != 2:
                    return None, None, False, \
                        f"Data is not 2D (ndim={data.ndim}) in {filepath}"
                
                # Check for NaN/Inf and handle
                if np.all(~np.isfinite(data)):
                    return None, None, False, \
                        f"All data is NaN/Inf in {filepath}"
                
                logger.info(f"Successfully loaded {filepath}: shape={data.shape}")
                return data, header, True, ""
        
        except (EOFError, OSError) as e:
            # Handles truncated files and corrupt FITS
            if 'corrupt' in str(e).lower() or 'eof' in str(e).lower():
                return None, None, False, f"Corrupted or truncated FITS file: {str(e)}"
            else:
                return None, None, False, f"File error: {str(e)}"
        except KeyError as e:
            return None, None, False, f"Missing keyword in header: {str(e)}"
        except Exception as e:
            return None, None, False, f"Error loading {filepath}: {str(e)}"
    
    @staticmethod
    def get_header_value(header, keys, default=None):
        """
        Safely get header value from multiple possible keys
        
        Parameters
        ----------
        header : dict
            FITS header
        keys : str or list
            Header key(s) to try
        default : any
            Default value if key not found
        
        Returns
        -------
        value : any
            Header value or default
        """
        if header is None:
            return default
        
        if isinstance(keys, str):
            keys = [keys]
        
        for key in keys:
            if key in header:
                return header[key]
        
        return default

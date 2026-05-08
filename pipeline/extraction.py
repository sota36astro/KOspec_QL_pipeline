"""
Spectral extraction module
Handles trace estimation and 1D spectrum extraction
"""

import numpy as np
import logging
from .utils import estimate_trace_from_profile, extract_1d_spectrum, \
    estimate_sky_level, sky_subtract

logger = logging.getLogger(__name__)


class SpectralExtraction:
    """Extract 1D spectra from 2D spectrum images"""
    
    def __init__(self, aperture_width=10, sky_offset=20, 
                 spatial_axis=0, trace_method='gaussian_fit'):
        """
        Initialize extraction parameters
        
        Parameters
        ----------
        aperture_width : int
            Width of extraction aperture in pixels
        sky_offset : int
            Offset from trace for sky estimation
        spatial_axis : int
            Spatial axis (0=rows, 1=cols)
        trace_method : str
            Method for trace estimation
        """
        self.aperture_width = aperture_width
        self.sky_offset = sky_offset
        self.spatial_axis = spatial_axis
        self.trace_method = trace_method
        self.trace_positions = {}
        self.extraction_info = {}
    
    def estimate_trace(self, image_2d, positive_peak_only=True, y_range=None):
        """
        Estimate spectral trace position
        
        Parameters
        ----------
        image_2d : 2D array
            2D spectrum image (can contain both positive and negative peaks)
        positive_peak_only : bool
            If True, trace only the positive peak
        y_range : tuple or None
            (y_min, y_max) range for trace estimation. If None, use full range.
        
        Returns
        -------
        trace_y : float or None
            Estimated trace position (in full image coordinates)
        spatial_profile : 1D array
            Spatial profile used for estimation
        """
        if image_2d is None:
            return None, None
        
        try:
            if self.spatial_axis == 0:
                # Spatial axis is rows (y-direction)
                image_height = image_2d.shape[0]
                
                if y_range is not None:
                    y_min, y_max = y_range
                    # Ensure range is within image bounds
                    y_min = max(0, min(y_min, image_height - 1))
                    y_max = max(y_min + 1, min(y_max, image_height))
                    
                    if y_max <= y_min:
                        logger.warning(f"Invalid y_range {y_range} for image height {image_height}, using full range")
                        image_subset = image_2d
                        y_offset = 0
                    else:
                        image_subset = image_2d[y_min:y_max, :]
                        y_offset = y_min
                else:
                    image_subset = image_2d
                    y_offset = 0
                
                # Collapse along cols
                if positive_peak_only:
                    # Use positive part only
                    positive = np.maximum(image_subset, 0)
                    spatial_profile = np.sum(positive, axis=1)
                else:
                    spatial_profile = np.sum(np.abs(image_subset), axis=1)
            else:
                # Spatial axis is cols
                if y_range is not None:
                    logger.warning("y_range parameter is only supported for spatial_axis=0")
                image_subset = image_2d
                y_offset = 0
                
                if positive_peak_only:
                    positive = np.maximum(image_subset, 0)
                    spatial_profile = np.sum(positive, axis=0)
                else:
                    spatial_profile = np.sum(np.abs(image_subset), axis=0)
            
            # Normalize
            if np.max(spatial_profile) > 0:
                spatial_profile = spatial_profile / np.max(spatial_profile)
            
            # Estimate trace
            trace_y_local = estimate_trace_from_profile(spatial_profile, 
                                                      method=self.trace_method)
            
            if trace_y_local is not None:
                # Convert back to full image coordinates
                trace_y = trace_y_local + y_offset
            else:
                trace_y = None
            
            logger.info(f"Estimated trace position: {trace_y:.1f} (y_range: {y_range})")
            return trace_y, spatial_profile
        
        except Exception as e:
            logger.error(f"Trace estimation failed: {str(e)}")
            return None, None
    
    def extract(self, image_2d, object_name=None, positive_peak_only=True):
        """
        Extract 1D spectrum with sky subtraction
        
        Parameters
        ----------
        image_2d : 2D array
            2D spectrum image
        object_name : str
            Name of the object (for logging)
        positive_peak_only : bool
            If True, trace only positive peak
        
        Returns
        -------
        spectrum_1d : 1D array or None
            Extracted 1D spectrum
        sky_spectrum : 1D array or None
            Sky spectrum
        trace_y : float or None
            Estimated trace position
        success : bool
            Whether extraction was successful
        error_msg : str
            Error message if unsuccessful
        """
        if image_2d is None:
            return None, None, None, False, "Input image is None"
        
        try:
            # Estimate trace
            trace_y, spatial_profile = self.estimate_trace(
                image_2d, positive_peak_only=positive_peak_only,
                y_range=(380, 620)  # Limit trace estimation to y=380-620
            )
            
            if trace_y is None:
                return None, None, None, False, "Trace estimation failed"
            
            # Extract spectrum using fixed-width aperture
            spectrum_1d = extract_1d_spectrum(
                image_2d, trace_y, 
                aperture_width=self.aperture_width,
                spatial_axis=self.spatial_axis
            )
            
            if spectrum_1d is None:
                return None, None, None, False, "Aperture extraction failed"
            
            # Estimate and subtract sky
            sky_spectrum = estimate_sky_level(
                image_2d, trace_y,
                aperture_width=self.aperture_width,
                sky_offset=self.sky_offset,
                spatial_axis=self.spatial_axis
            )
            
            spectrum_sky_subtracted = sky_subtract(
                spectrum_1d, sky_spectrum, self.aperture_width
            )
            
            # Store info
            if object_name:
                self.trace_positions[object_name] = trace_y
                self.extraction_info[object_name] = {
                    'trace_y': trace_y,
                    'aperture_width': self.aperture_width,
                    'sky_offset': self.sky_offset,
                    'spatial_profile': spatial_profile
                }
            
            logger.info(f"Successfully extracted spectrum for {object_name}")
            return spectrum_sky_subtracted, sky_spectrum, trace_y, True, ""
        
        except Exception as e:
            error_msg = f"Extraction error: {str(e)}"
            logger.error(error_msg)
            return None, None, None, False, error_msg
    
    def get_trace_position(self, object_name):
        """Get stored trace position for object"""
        return self.trace_positions.get(object_name)
    
    def get_extraction_info(self, object_name):
        """Get stored extraction info for object"""
        return self.extraction_info.get(object_name)

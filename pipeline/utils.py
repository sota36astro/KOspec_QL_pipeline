"""
Utility functions for the KOspec pipeline
"""

import numpy as np
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def setup_logging(verbose=False):
    """Setup logging configuration"""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )


def safe_log(func):
    """Decorator to safely log errors without stopping pipeline"""
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            logger.error(f"Error in {func.__name__}: {str(e)}", exc_info=True)
            return None
    return wrapper


def estimate_trace_from_profile(spatial_profile, 
                                method='gaussian_fit', 
                                smooth_window=5):
    """
    Estimate trace position from spatial profile
    
    Parameters
    ----------
    spatial_profile : 1D array
        Spatial profile (collapsed spectrum)
    method : str
        Method to use: 'center_of_mass', 'peak', 'gaussian_fit'
    smooth_window : int
        Window size for smoothing before fitting
    
    Returns
    -------
    trace_y : float or array
        Estimated y position(s) of the trace
    """
    if spatial_profile is None or len(spatial_profile) == 0:
        return None
    
    # Smooth the profile
    from scipy.ndimage import uniform_filter1d
    smoothed = uniform_filter1d(spatial_profile, size=smooth_window, mode='nearest')
    
    if method == 'peak':
        return np.argmax(smoothed)
    elif method == 'center_of_mass':
        # Use center of mass above threshold
        threshold = np.max(smoothed) * 0.1
        masked_profile = np.where(smoothed > threshold, smoothed, 0)
        if np.sum(masked_profile) > 0:
            y_indices = np.arange(len(smoothed))
            return np.average(y_indices, weights=masked_profile)
        else:
            return np.argmax(smoothed)
    elif method == 'gaussian_fit':
        try:
            from scipy.optimize import curve_fit
            y = np.arange(len(smoothed))
            # Initial guess
            peak_idx = np.argmax(smoothed)
            peak_val = smoothed[peak_idx]
            # Estimate sigma from profile width
            above_half = np.where(smoothed > peak_val / 2)[0]
            if len(above_half) > 1:
                sigma_guess = (above_half[-1] - above_half[0]) / 2.355
            else:
                sigma_guess = 3.0
            
            def gaussian(x, amp, mu, sigma, offset):
                return amp * np.exp(-0.5 * ((x - mu) / sigma) ** 2) + offset
            
            popt, _ = curve_fit(gaussian, y, smoothed,
                              p0=[peak_val, peak_idx, sigma_guess, 0],
                              maxfev=1000)
            return popt[1]  # Return center position
        except Exception as e:
            logger.warning(f"Gaussian fit failed: {e}. Using peak.")
            return np.argmax(smoothed)
    else:
        return np.argmax(smoothed)


def extract_1d_spectrum(image_2d, trace_y, aperture_width=5, 
                       spatial_axis=0):
    """
    Extract 1D spectrum using fixed-width aperture
    
    Parameters
    ----------
    image_2d : 2D array
        2D spectrum image
    trace_y : float or array
        Trace position(s)
    aperture_width : int
        Width of extraction aperture in pixels
    spatial_axis : int
        Spatial axis (0=rows, 1=cols)
    
    Returns
    -------
    spectrum_1d : 1D array
        Extracted 1D spectrum
    """
    if trace_y is None:
        return None
    
    trace_y = int(np.round(trace_y))
    half_width = aperture_width // 2
    
    if spatial_axis == 0:
        # Spatial axis is rows, wavelength axis is cols
        y_start = max(0, trace_y - half_width)
        y_end = min(image_2d.shape[0], trace_y + half_width + 1)
        spectrum = np.sum(image_2d[y_start:y_end, :], axis=0)
    else:
        # Spatial axis is cols
        x_start = max(0, trace_y - half_width)
        x_end = min(image_2d.shape[1], trace_y + half_width + 1)
        spectrum = np.sum(image_2d[:, x_start:x_end], axis=1)
    
    return spectrum


def estimate_sky_level(image_2d, trace_y, aperture_width=5, 
                       sky_offset=20, spatial_axis=0):
    """
    Estimate sky level from regions away from trace
    
    Parameters
    ----------
    image_2d : 2D array
        2D spectrum image
    trace_y : float
        Trace position
    aperture_width : int
        Width of extraction aperture
    sky_offset : int
        Offset from trace to sky region
    spatial_axis : int
        Spatial axis (0=rows, 1=cols)
    
    Returns
    -------
    sky_spectrum : 1D array
        Sky spectrum
    """
    if trace_y is None:
        return None
    
    trace_y = int(np.round(trace_y))
    half_width = aperture_width // 2
    aperture_pixels = 2 * half_width + 1
    
    if spatial_axis == 0:
        # Use all valid sky regions and scale the mean sky level to the
        # extraction aperture. This returns the sky contribution in summed
        # aperture units, matching extract_1d_spectrum().
        y_top_start = max(0, trace_y - half_width - sky_offset - half_width)
        y_top_end = max(0, trace_y - half_width - sky_offset)
        y_bot_start = min(image_2d.shape[0], trace_y + half_width + sky_offset)
        y_bot_end = min(image_2d.shape[0], 
                        trace_y + half_width + sky_offset + half_width)

        sky_regions = []
        if y_top_end - y_top_start > 0:
            sky_regions.append(image_2d[y_top_start:y_top_end, :])
        if y_bot_end - y_bot_start > 0:
            sky_regions.append(image_2d[y_bot_start:y_bot_end, :])

        if sky_regions:
            sky_pixels = np.concatenate(sky_regions, axis=0)
            sky = np.mean(sky_pixels, axis=0) * aperture_pixels
        else:
            sky = np.zeros(image_2d.shape[1])
    else:
        # Similar for column-wise
        x_left_start = max(0, trace_y - half_width - sky_offset - half_width)
        x_left_end = max(0, trace_y - half_width - sky_offset)
        x_right_start = min(image_2d.shape[1], trace_y + half_width + sky_offset)
        x_right_end = min(image_2d.shape[1],
                         trace_y + half_width + sky_offset + half_width)

        sky_regions = []
        if x_left_end - x_left_start > 0:
            sky_regions.append(image_2d[:, x_left_start:x_left_end])
        if x_right_end - x_right_start > 0:
            sky_regions.append(image_2d[:, x_right_start:x_right_end])

        if sky_regions:
            sky_pixels = np.concatenate(sky_regions, axis=1)
            sky = np.mean(sky_pixels, axis=1) * aperture_pixels
        else:
            sky = np.zeros(image_2d.shape[0])
    
    return sky


def sky_subtract(spectrum, sky_spectrum, aperture_width):
    """
    Subtract sky from spectrum
    
    Parameters
    ----------
    spectrum : 1D array
        Extracted spectrum
    sky_spectrum : 1D array
        Sky spectrum
    aperture_width : int
        Aperture width used in extraction
    
    Returns
    -------
    sky_subtracted : 1D array
        Sky-subtracted spectrum
    """
    if sky_spectrum is None:
        return spectrum
    
    return spectrum - sky_spectrum

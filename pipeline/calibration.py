"""
Wavelength calibration module
"""

import numpy as np
import logging

logger = logging.getLogger(__name__)

DEFAULT_POLY_VALID_NPIX = 1056


# Wavelengths are stored in Angstrom throughout the pipeline.
# Primary emission lines
EMISSION_LINES = {
    # Hydrogen Balmer series
    'Halpha': 6563.0,
    'Hbeta': 4861.0,
    'Hgamma': 4340.0,
    'Hdelta': 4102.0,
    'H_Balmer_alpha': 6563.0,
    'H_Balmer_beta': 4861.0,
    'H_Balmer_gamma': 4340.0,
    'H_Balmer_delta': 4102.0,
    'H_Balmer_epsilon': 3970.0,
    'H_Balmer_zeta': 3889.0,
    'H_Balmer_eta': 3835.0,
    
    # Helium
    'He_I_3889': 3889.0,
    'He_I_4026': 4026.0,
    'He_I_4471': 4471.0,
    'He_I_4713': 4713.0,
    'He_I_4922': 4922.0,
    'He_I_5016': 5016.0,
    'He_I_5875': 5876.0,
    'He_I_6678': 6678.0,
    'He_I_7065': 7065.0,
    'He_II_4200': 4200.0,
    'He_II_4542': 4542.0,
    'He_II_4686': 4686.0,
    'He_II_5411': 5411.0,
    'He_II_6560': 6560.0,
    
    # Sodium
    'Na_I_D': 5890.0,  # D lines (average)
    'Na_I_D1': 5889.95,
    'Na_I_D2': 5895.92,
    
    # Oxygen
    '[O_I]_6300': 6300.0,
    '[O_I]_6364': 6364.0,
    '[O_III]_5007': 5007.0,
    '[O_III]_4959': 4959.0,
    
    # Calcium
    'Ca_II_3968': 3968.0,
    'Ca_II_3933': 3933.0,
    'Ca_II_NIR_8498': 8498.0,
    'Ca_II_NIR_8542': 8542.0,
    'Ca_II_NIR_8662': 8662.0,
    
    # Sulfur
    '[S_II]_6717': 6717.0,
    '[S_II]_6731': 6731.0,
    
    # Nitrogen
    'N_III_4097': 4097.0,
    'N_III_4634': 4634.0,
    'N_III_4640': 4640.0,
    'N_III_4642': 4642.0,
    'N_IV_4058': 4058.0,
    'N_IV_7111': 7111.0,
    'N_IV_7123': 7123.0,
    'N_V_4604': 4604.0,
    'N_V_4620': 4620.0,
    '[N_II]_6548': 6548.0,
    '[N_II]_6584': 6584.0,

    # Carbon
    'C_III_4647': 4647.0,
    'C_III_4650': 4650.0,
    'C_III_5696': 5696.0,
    'C_IV_5801': 5801.0,
    'C_IV_5812': 5812.0,
}

# Lines used by default
DEFAULT_LINES = [
    'H_Balmer', 'He_I', 'He_II', 'Na_I_D', '[O_I]_6300', 'Ca_II_NIR'
]


class WavelengthCalibration:
    """Wavelength calibration for spectra"""
    
    def __init__(self, wave_solution=None, ctype='PIXEL'):
        """
        Initialize wavelength calibration
        
        Parameters
        ----------
        wave_solution : dict or None
            Wavelength solution parameters
            For polynomial: {'type': 'poly', 'coeffs': [c0, c1, c2, ...]}
            For linear: {'type': 'linear', 'wave0': w0, 'disp': disp, 'pix0': p0}
        ctype : str
            Calibration type (PIXEL, LINEAR, POLYNOM, etc.)
        """
        if wave_solution is None:
            # Default: polynomial solution converted from nm to Angstrom.
            # f_nm(x) = 1043.947468362323 - 2.491445451495574*x
            #         + 0.006468775026217139*x^2
            #         - 0.00001183580177752003*x^3
            #         + 0.00000001365974323268898*x^4
            #         - 0.00000000000873649357565434*x^5
            #         + 0.000000000000002335710049781674*x^6
            self.wave_solution = {
                'type': 'poly',
                'valid_npix': DEFAULT_POLY_VALID_NPIX,
                'coeffs': [
                    10439.47468362323,       # c0
                    -24.91445451495574,      # c1
                    6.468775026217139e-2,    # c2
                    -1.183580177752003e-4,   # c3
                    1.365974323268898e-7,    # c4
                    -8.73649357565434e-11,   # c5
                    2.335710049781674e-14    # c6
                ]
            }
        else:
            self.wave_solution = wave_solution
        
        self.ctype = ctype
        self._warned_outside_valid_range = False
    
    def pixel_to_wavelength(self, pixels):
        """
        Convert pixel coordinates to wavelength
        
        Parameters
        ----------
        pixels : scalar or array
            Pixel coordinate(s)
        
        Returns
        -------
        wavelength : scalar or array
            Wavelength in Angstrom
        """
        sol = self.wave_solution
        
        if sol['type'] == 'linear':
            # Linear: w = w0 + disp * (pix - pix0)
            return sol['wave0'] + sol['disp'] * (pixels - sol['pix0'])
        
        elif sol['type'] == 'poly':
            # Polynomial: w = sum(c_i * pix^i)
            valid_npix = sol.get('valid_npix')
            if valid_npix is not None and not self._warned_outside_valid_range:
                pixel_values = np.asarray(pixels)
                if np.any((pixel_values < 0) | (pixel_values >= valid_npix)):
                    logger.warning(
                        "Polynomial wavelength solution is calibrated for "
                        f"pixels 0-{valid_npix - 1}; requested pixels extend "
                        f"to {np.nanmin(pixel_values):.0f}-{np.nanmax(pixel_values):.0f}"
                    )
                    self._warned_outside_valid_range = True
            coeffs = sol['coeffs']
            return np.polyval(coeffs[::-1], pixels)
        
        else:
            logger.warning(f"Unknown solution type: {sol['type']}")
            return pixels
    
    def set_linear_solution(self, wave0, pix0, disp):
        """
        Set linear wavelength solution
        
        Parameters
        ----------
        wave0 : float
            Reference wavelength (Angstrom)
        pix0 : float
            Reference pixel
        disp : float
            Dispersion (Angstrom/pixel)
        """
        self.wave_solution = {
            'type': 'linear',
            'wave0': wave0,
            'pix0': pix0,
            'disp': disp
        }
        logger.info(f"Set linear solution: w0={wave0}, pix0={pix0}, disp={disp}")
    
    def set_polynomial_solution(self, coeffs):
        """
        Set polynomial wavelength solution
        
        Parameters
        ----------
        coeffs : list
            Polynomial coefficients [c0, c1, c2, ...] where w = sum(c_i * pix^i)
        """
        self.wave_solution = {
            'type': 'poly',
            'coeffs': coeffs
        }
        logger.info(f"Set polynomial solution with {len(coeffs)} coefficients")
    
    def apply_redshift(self, wavelength, redshift):
        """
        Apply redshift to wavelengths
        
        Parameters
        ----------
        wavelength : scalar or array
            Wavelength(s) in Angstrom
        redshift : float
            Redshift z
        
        Returns
        -------
        redshifted : scalar or array
            Redshifted wavelength(s)
        """
        return wavelength * (1 + redshift)
    
    def get_line_wavelengths(self, redshift=0.0, line_list=None):
        """
        Get emission line wavelengths
        
        Parameters
        ----------
        redshift : float
            Redshift z for observed frame
        line_list : list or None
            List of line names. If None, use default list.
            Prefixes such as 'Ca_II_NIR' expand to matching lines.
        
        Returns
        -------
        lines : dict
            Dictionary mapping line name to observed wavelength
        """
        if line_list is None:
            line_list = DEFAULT_LINES
        
        lines = {}
        for line_name in line_list:
            if line_name in EMISSION_LINES:
                rest_wave = EMISSION_LINES[line_name]
                obs_wave = self.apply_redshift(rest_wave, redshift)
                lines[line_name] = obs_wave
                continue

            # Handle prefix groups such as 'He_I', 'Ca_II_NIR', or '[O_I]'.
            # Require the next character to be '_' so He_I does not match He_II.
            group_prefix = f"{line_name}_"
            matching = {k: v for k, v in EMISSION_LINES.items()
                        if k.startswith(group_prefix)}
            for k, v in matching.items():
                obs_wave = self.apply_redshift(v, redshift)
                lines[k] = obs_wave
        
        return lines
    
    @staticmethod
    def create_default_solution(npix=1024, wave_start=4000, wave_end=9000):
        """
        Create a default wavelength solution
        
        Parameters
        ----------
        npix : int
            Number of pixels
        wave_start : float
            Starting wavelength (Angstrom)
        wave_end : float
            Ending wavelength (Angstrom)
        
        Returns
        -------
        calib : WavelengthCalibration
            Calibration object
        """
        # Use the user-specified polynomial solution as default
        calib = WavelengthCalibration()
        calib.wave_solution['npix'] = npix
        return calib

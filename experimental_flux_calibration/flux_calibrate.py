#!/usr/bin/env python3
"""
Experimental flux calibration from quicklook 1D spectra.

This module is kept outside the stable quicklook pipeline on purpose.
It consumes quicklook text spectra and a reference standard-star flux table.
"""

import argparse
import logging
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from pipeline.runtime import configure_matplotlib_cache
configure_matplotlib_cache()

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.ndimage import median_filter

from pipeline.calibration import WavelengthCalibration


logger = logging.getLogger(__name__)


TELLURIC_BANDS = [
    (6860.0, 6950.0),
    (7160.0, 7340.0),
    (7590.0, 7700.0),
    (8120.0, 8400.0),
    (8900.0, 9800.0),
]
TELLURIC_BAND_LABELS = ['O2 B', 'H2O', 'O2 A', 'H2O', 'H2O']

STANDARD_TEMPLATES = {
    'HR4468': {
        'teff': 10500.0,
        'v_mag': 4.70,
        'wave_min': 3500.0,
        'wave_max': 10500.0,
        'wave_step': 5.0,
        'description': 'Approximate B9.5V blackbody continuum template',
    },
}

STANDARD_FEATURE_MASKS = [
    (4080.0, 4130.0),  # Hdelta
    (4310.0, 4370.0),  # Hgamma
    (4820.0, 4900.0),  # Hbeta
    (6510.0, 6620.0),  # Halpha
]

V_ZEROPOINT_FLAMBDA = 3.631e-9  # erg s^-1 cm^-2 A^-1 near 5500 A
V_REFERENCE_WAVELENGTH = 5500.0
DEFAULT_VALID_WAVE_MIN = 4500.0
DEFAULT_VALID_WAVE_MAX = 8500.0
DEFAULT_SPECTRA_DIR = 'spectra'
AIRMASS_KEYS = ('AIRMASS', 'SECZ', 'AM', 'AMSTART', 'AIRM-STR')
LINE_ELEMENT_COLORS = {
    'H': 'crimson',
    'He': 'darkorange',
    'Na': 'goldenrod',
    'O': 'forestgreen',
    'Ca': 'purple',
    'S': 'teal',
    'N': 'royalblue',
    'C': 'darkcyan',
}
DEFAULT_EXTINCTION_TABLE = np.asarray([
    [3500.0, 0.55],
    [4000.0, 0.32],
    [4500.0, 0.22],
    [5000.0, 0.16],
    [5500.0, 0.13],
    [6000.0, 0.11],
    [6500.0, 0.09],
    [7000.0, 0.08],
    [7500.0, 0.07],
    [8000.0, 0.06],
    [8500.0, 0.055],
    [9000.0, 0.05],
    [10000.0, 0.045],
], dtype=float)


def load_quicklook_spectrum(path):
    """Load quicklook 1D text spectrum."""
    path = Path(path)
    data = np.loadtxt(path)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Expected at least two columns in {path}")

    wavelength = data[:, 0]
    flux = data[:, 1]
    sky = data[:, 2] if data.shape[1] >= 3 else None
    return wavelength, flux, sky


def load_reference_flux(path):
    """Load reference standard-star flux table."""
    path = Path(path)
    data = np.loadtxt(path)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Expected at least two columns in {path}")
    return data[:, 0], data[:, 1]


def object_name_from_quicklook_path(path):
    """Infer object name from a quicklook 1D text filename."""
    name = Path(path).stem
    if name.endswith('_1d'):
        name = name[:-3]
    return name


def normalize_object_name(name):
    """Normalize FITS OBJECT names by dropping A/B suffixes."""
    if name is None:
        return ''
    name = str(name).strip()
    for suffix in ('_A', '_B'):
        if name.endswith(suffix):
            return name[:-2]
    return name


def read_exptime_from_fits(path):
    """Read EXPTIME from one FITS file."""
    path = Path(path)
    with fits.open(path) as hdul:
        header = hdul[0].header
        exptime = header.get('EXPTIME')
    if exptime is None:
        raise ValueError(f"EXPTIME not found in {path}")
    exptime = float(exptime)
    if exptime <= 0:
        raise ValueError(f"EXPTIME must be positive in {path}: {exptime}")
    return exptime


def read_airmass_from_fits(path):
    """Read airmass from one FITS file."""
    path = Path(path)
    with fits.open(path) as hdul:
        header = hdul[0].header
        for key in AIRMASS_KEYS:
            if key in header:
                airmass = float(header[key])
                if airmass <= 0:
                    raise ValueError(
                        f"{key} must be positive in {path}: {airmass}"
                    )
                return airmass

        altitude = header.get('ALTITUDE')
        if altitude is not None:
            altitude = float(altitude)
            if 0 < altitude <= 90:
                zenith_distance = np.deg2rad(90.0 - altitude)
                airmass = 1.0 / np.cos(zenith_distance)
                logger.debug(
                    "Computed airmass %.6f from ALTITUDE in %s",
                    airmass, path
                )
                return float(airmass)

    raise ValueError(f"AIRMASS not found in {path}")


def fits_position(path):
    """Infer A/B position from FITS header or filename."""
    path = Path(path)
    try:
        with fits.open(path) as hdul:
            header = hdul[0].header
            object_name = str(header.get('OBJECT', '')).strip()
            position = str(header.get('POSITION', '')).strip().upper()
    except Exception:
        object_name = ''
        position = ''

    if object_name.endswith('_A') or position == 'A':
        return 'A'
    if object_name.endswith('_B') or position == 'B':
        return 'B'
    if '_A' in path.name:
        return 'A'
    if '_B' in path.name:
        return 'B'
    return None


def find_object_fits(object_name, spectra_dir=DEFAULT_SPECTRA_DIR):
    """Find FITS files in spectra_dir whose OBJECT matches object_name."""
    spectra_dir = Path(spectra_dir)
    matches = []
    if not spectra_dir.exists():
        return matches

    for path in sorted(spectra_dir.glob('*.fits')):
        try:
            with fits.open(path) as hdul:
                header = hdul[0].header
                header_object = normalize_object_name(header.get('OBJECT'))
                if header_object == object_name:
                    matches.append(path)
        except Exception as exc:
            logger.debug("Skipping %s while searching FITS headers: %s", path, exc)
    return matches


def resolve_fits_paths(label, quicklook_path, fits_paths,
                       spectra_dir=DEFAULT_SPECTRA_DIR):
    """Resolve FITS paths and prefer the positive A frame used by quicklook."""
    paths = [Path(path) for path in fits_paths]
    object_name = object_name_from_quicklook_path(quicklook_path)
    if not paths:
        paths = find_object_fits(object_name, spectra_dir=spectra_dir)

    if not paths:
        raise ValueError(
            f"Could not find {label} FITS files. Provide --{label}-fits."
        )

    positive_paths = [path for path in paths if fits_position(path) == 'A']
    if positive_paths:
        return positive_paths, 'positive A-frame'
    return paths, 'all matched FITS files'


def resolve_exptime(label, quicklook_path, explicit_exptime, fits_paths,
                    spectra_dir=DEFAULT_SPECTRA_DIR):
    """Resolve exposure time from CLI value, FITS paths, or spectra_dir search."""
    if explicit_exptime is not None:
        if explicit_exptime <= 0:
            raise ValueError(f"{label} exposure time must be positive")
        logger.info("Using %s exposure time from CLI: %.3f s", label, explicit_exptime)
        return explicit_exptime

    try:
        paths_for_exptime, exposure_source = resolve_fits_paths(
            label, quicklook_path, fits_paths, spectra_dir=spectra_dir
        )
    except ValueError as exc:
        raise ValueError(
            f"Could not determine {label} exposure time. "
            f"Provide --{label}-fits or --{label}-exptime."
        ) from exc

    exptimes = np.asarray(
        [read_exptime_from_fits(path) for path in paths_for_exptime],
        dtype=float
    )
    if np.nanmax(exptimes) - np.nanmin(exptimes) > 1e-6:
        logger.warning(
            "%s FITS files used for EXPTIME have different values %s; using their sum",
            label, ', '.join(f"{value:.3f}" for value in exptimes)
        )

    total_exptime = float(np.nansum(exptimes))
    logger.info(
        "Using %s exposure time from %s header(s): %.3f s (%s)",
        label, exposure_source, total_exptime,
        ', '.join(path.name for path in paths_for_exptime)
    )
    return total_exptime


def resolve_airmass(label, quicklook_path, explicit_airmass, fits_paths,
                    spectra_dir=DEFAULT_SPECTRA_DIR):
    """Resolve airmass from CLI value, FITS paths, or spectra_dir search."""
    if explicit_airmass is not None:
        if explicit_airmass <= 0:
            raise ValueError(f"{label} airmass must be positive")
        logger.info("Using %s airmass from CLI: %.6f", label, explicit_airmass)
        return explicit_airmass

    try:
        paths_for_airmass, airmass_source = resolve_fits_paths(
            label, quicklook_path, fits_paths, spectra_dir=spectra_dir
        )
    except ValueError as exc:
        raise ValueError(
            f"Could not determine {label} airmass. "
            f"Provide --{label}-fits or --{label}-airmass."
        ) from exc

    airmasses = np.asarray(
        [read_airmass_from_fits(path) for path in paths_for_airmass],
        dtype=float
    )
    if np.nanmax(airmasses) - np.nanmin(airmasses) > 1e-3:
        logger.warning(
            "%s FITS files used for AIRMASS have different values %s; "
            "using exposure-time-unweighted mean",
            label, ', '.join(f"{value:.6f}" for value in airmasses)
        )

    airmass = float(np.nanmean(airmasses))
    logger.info(
        "Using %s airmass from %s header(s): %.6f (%s)",
        label, airmass_source, airmass,
        ', '.join(path.name for path in paths_for_airmass)
    )
    return airmass


def load_extinction_curve(path=None):
    """Load atmospheric extinction curve in mag/airmass."""
    if path is None:
        return DEFAULT_EXTINCTION_TABLE[:, 0], DEFAULT_EXTINCTION_TABLE[:, 1]

    data = np.loadtxt(path)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Expected at least two columns in {path}")
    return data[:, 0], data[:, 1]


def extinction_at_wavelength(wavelength, extinction_wave, extinction_mag):
    """Interpolate atmospheric extinction in mag/airmass."""
    extinction_wave, extinction_mag = sort_by_wavelength(
        extinction_wave, extinction_mag
    )
    return np.interp(
        wavelength, extinction_wave, extinction_mag,
        left=extinction_mag[0], right=extinction_mag[-1]
    )


def differential_airmass_correction(wavelength, standard_airmass,
                                    target_airmass, extinction_wave,
                                    extinction_mag):
    """Return factor correcting target flux for target-standard airmass difference."""
    delta_airmass = target_airmass - standard_airmass
    extinction = extinction_at_wavelength(
        wavelength, extinction_wave, extinction_mag
    )
    return 10.0 ** (0.4 * extinction * delta_airmass)


def planck_lambda_angstrom(wavelength, temperature):
    """Return Planck function shape per Angstrom for wavelength in Angstrom."""
    wavelength_cm = wavelength * 1e-8
    h = 6.62607015e-27
    c = 2.99792458e10
    k = 1.380649e-16
    exponent = h * c / (wavelength_cm * k * temperature)
    intensity_per_cm = (2.0 * h * c**2) / (
        wavelength_cm**5 * np.expm1(exponent)
    )
    return intensity_per_cm * 1e-8


def make_blackbody_template(wavelength, temperature, v_mag):
    """Build a blackbody continuum normalized to an approximate V magnitude."""
    flux_shape = planck_lambda_angstrom(wavelength, temperature)
    reference_shape = planck_lambda_angstrom(
        np.asarray([V_REFERENCE_WAVELENGTH]), temperature
    )[0]
    reference_flux = V_ZEROPOINT_FLAMBDA * 10 ** (-0.4 * v_mag)
    return flux_shape / reference_shape * reference_flux


def load_template_reference(template_name, wave_min=None, wave_max=None,
                            wave_step=None, teff=None, v_mag=None):
    """Load a built-in approximate standard-star reference template."""
    template_key = template_name.upper()
    if template_key not in STANDARD_TEMPLATES:
        choices = ', '.join(sorted(STANDARD_TEMPLATES))
        raise ValueError(f"Unknown template {template_name!r}. Choices: {choices}")

    template = STANDARD_TEMPLATES[template_key]
    wave_min = template['wave_min'] if wave_min is None else wave_min
    wave_max = template['wave_max'] if wave_max is None else wave_max
    wave_step = template['wave_step'] if wave_step is None else wave_step
    teff = template['teff'] if teff is None else teff
    v_mag = template['v_mag'] if v_mag is None else v_mag

    wavelength = np.arange(wave_min, wave_max + wave_step, wave_step)
    flux = make_blackbody_template(wavelength, teff, v_mag)
    return wavelength, flux


def sort_by_wavelength(wavelength, *arrays):
    """Return wavelength and arrays sorted by increasing wavelength."""
    order = np.argsort(wavelength)
    sorted_arrays = [arr[order] if arr is not None else None for arr in arrays]
    return (wavelength[order], *sorted_arrays)


def telluric_mask(wavelength):
    """Return True for wavelengths outside common telluric bands."""
    mask = np.ones_like(wavelength, dtype=bool)
    for wave_min, wave_max in TELLURIC_BANDS:
        mask &= ~((wavelength >= wave_min) & (wavelength <= wave_max))
    return mask


def standard_feature_mask(wavelength):
    """Return True outside strong standard-star absorption regions."""
    mask = np.ones_like(wavelength, dtype=bool)
    for wave_min, wave_max in STANDARD_FEATURE_MASKS:
        mask &= ~((wavelength >= wave_min) & (wavelength <= wave_max))
    return mask


def wavelength_range_mask(wavelength, wave_min, wave_max):
    """Return True inside the trusted wavelength range."""
    return (wavelength >= wave_min) & (wavelength <= wave_max)


def interpolate_reference(reference_wave, reference_flux, wavelength):
    """Interpolate reference flux onto spectrum wavelength grid."""
    ref_wave, ref_flux = sort_by_wavelength(reference_wave, reference_flux)
    interp_flux = np.interp(wavelength, ref_wave, ref_flux, left=np.nan, right=np.nan)
    return interp_flux


def smooth_sensitivity(sensitivity, valid_mask, smooth_window):
    """Median-smooth sensitivity using only valid points."""
    smoothed = np.full_like(sensitivity, np.nan, dtype=float)
    if np.count_nonzero(valid_mask) < 2:
        return smoothed

    indices = np.arange(len(sensitivity))
    filled = np.interp(indices, indices[valid_mask], sensitivity[valid_mask])

    if smooth_window > 1:
        if smooth_window % 2 == 0:
            smooth_window += 1
        filled = median_filter(filled, size=smooth_window, mode='nearest')

    smoothed[valid_mask] = filled[valid_mask]
    smoothed[~valid_mask] = filled[~valid_mask]
    return smoothed


def build_sensitivity(standard_wave, standard_counts, reference_wave,
                      reference_flux, standard_exptime=1.0,
                      smooth_window=51, mask_telluric=True,
                      mask_standard_features=True,
                      min_standard_count_frac=0.02,
                      valid_wave_min=DEFAULT_VALID_WAVE_MIN,
                      valid_wave_max=DEFAULT_VALID_WAVE_MAX):
    """Build a smoothed sensitivity curve from a standard star."""
    if standard_exptime <= 0:
        raise ValueError("standard_exptime must be positive")

    standard_wave, standard_counts = sort_by_wavelength(
        standard_wave, standard_counts
    )
    reference_interp = interpolate_reference(
        reference_wave, reference_flux, standard_wave
    )
    counts_rate = standard_counts / standard_exptime

    raw_sensitivity = np.full_like(standard_wave, np.nan, dtype=float)
    valid = (
        np.isfinite(standard_wave)
        & np.isfinite(reference_interp)
        & np.isfinite(counts_rate)
        & (reference_interp > 0)
        & (counts_rate > 0)
    )
    if mask_telluric:
        valid &= telluric_mask(standard_wave)
    if mask_standard_features:
        valid &= standard_feature_mask(standard_wave)
    valid &= wavelength_range_mask(standard_wave, valid_wave_min, valid_wave_max)
    if min_standard_count_frac > 0:
        positive_counts = counts_rate[np.isfinite(counts_rate) & (counts_rate > 0)]
        if positive_counts.size:
            count_floor = (
                np.percentile(positive_counts, 95) * min_standard_count_frac
            )
            valid &= counts_rate > count_floor

    if np.count_nonzero(valid) < 2:
        raise ValueError(
            "Not enough valid standard-star points to build sensitivity. "
            "Check reference wavelength coverage and standard counts."
        )

    raw_sensitivity[valid] = reference_interp[valid] / counts_rate[valid]
    sensitivity = smooth_sensitivity(raw_sensitivity, valid, smooth_window)
    return standard_wave, raw_sensitivity, sensitivity, valid


def apply_sensitivity(target_wave, target_counts, sensitivity_wave,
                      sensitivity, target_exptime=1.0,
                      airmass_correction=None):
    """Apply a sensitivity curve to a target spectrum."""
    if target_exptime <= 0:
        raise ValueError("target_exptime must be positive")

    target_counts_rate = target_counts / target_exptime
    sens_interp = np.interp(
        target_wave, sensitivity_wave, sensitivity,
        left=np.nan, right=np.nan
    )
    calibrated_flux = target_counts_rate * sens_interp
    if airmass_correction is not None:
        calibrated_flux = calibrated_flux * airmass_correction
    return calibrated_flux, sens_interp


def save_calibrated_spectrum(output_path, wavelength, calibrated_flux,
                             target_counts, sky, sensitivity,
                             airmass_correction=None):
    """Save calibrated target spectrum."""
    columns = [wavelength, calibrated_flux, target_counts]
    names = ["Wavelength(Angstrom)", "Flux_calibrated", "Target_counts"]
    if sky is not None:
        columns.append(sky)
        names.append("BG_counts")
    columns.append(sensitivity)
    names.append("Sensitivity")
    if airmass_correction is not None:
        columns.append(airmass_correction)
        names.append("Airmass_correction")

    data = np.column_stack(columns)
    header = "  ".join(names)
    np.savetxt(output_path, data, header=header, fmt='%.8e', comments='# ')


def save_airmass_correction(output_path, wavelength, correction,
                            extinction_mag, standard_airmass, target_airmass):
    """Save the differential airmass correction applied to target flux."""
    data = np.column_stack([wavelength, correction, extinction_mag])
    header = (
        "Wavelength(Angstrom)  Airmass_correction  "
        "Extinction_mag_per_airmass\n"
        f"standard_airmass={standard_airmass:.6f} "
        f"target_airmass={target_airmass:.6f} "
        f"delta_airmass={target_airmass - standard_airmass:.6f}"
    )
    np.savetxt(output_path, data, header=header, fmt='%.8e', comments='# ')


def plot_airmass_correction(output_path, wavelength, correction,
                            extinction_mag, standard_airmass, target_airmass,
                            valid_wave_min, valid_wave_max):
    """Plot the differential airmass correction."""
    fig, (ax_corr, ax_ext) = plt.subplots(
        2, 1, figsize=(12, 7), sharex=True,
        gridspec_kw={'height_ratios': [1.1, 0.9]}
    )
    for ax in (ax_corr, ax_ext):
        shade_untrusted_wavelengths(
            ax, np.nanmin(wavelength), np.nanmax(wavelength),
            valid_wave_min, valid_wave_max
        )

    delta_airmass = target_airmass - standard_airmass
    ax_corr.plot(wavelength, correction, 'k-', linewidth=1.5)
    ax_corr.axhline(1.0, color='0.5', linestyle='--', linewidth=1)
    ax_corr.set_ylabel('Flux correction')
    ax_corr.set_title(
        f'Differential airmass correction: target - standard = '
        f'{delta_airmass:.3f}'
    )
    ax_corr.grid(True, alpha=0.3)

    ax_ext.plot(wavelength, extinction_mag, color='tab:blue', linewidth=1.2)
    ax_ext.set_xlabel('Wavelength (Angstrom)')
    ax_ext.set_ylabel('k(lambda)\nmag/airmass')
    ax_ext.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(output_path, dpi=120, bbox_inches='tight')
    plt.close(fig)


def save_sensitivity(output_path, wavelength, raw_sensitivity, sensitivity,
                     valid, reference_flux, counts_rate):
    """Save sensitivity curve and the quantities used to build it."""
    data = np.column_stack([
        wavelength, raw_sensitivity, sensitivity, valid.astype(int),
        reference_flux, counts_rate
    ])
    header = (
        "Wavelength(Angstrom)  Raw_sensitivity  Smoothed_sensitivity  "
        "Valid_for_fit  Template_or_reference_flux  Standard_counts_rate"
    )
    np.savetxt(output_path, data, header=header, fmt='%.8e', comments='# ')


def set_robust_ylim(ax, values, lower=1, upper=99):
    """Set y limits from robust percentiles for display only."""
    finite = np.asarray(values)
    finite = finite[np.isfinite(finite)]
    if finite.size < 2:
        return

    y_min, y_max = np.percentile(finite, [lower, upper])
    if y_min == y_max:
        return

    pad = 0.08 * (y_max - y_min)
    ax.set_ylim(y_min - pad, y_max + pad)


def line_element(line_name):
    """Return the element key used to color a line ID."""
    clean_name = line_name.strip('[]')
    if clean_name.startswith('He'):
        return 'He'
    if clean_name.startswith('H'):
        return 'H'
    return clean_name.split('_', 1)[0]


def line_color(line_name):
    """Return the plot color for a line ID."""
    return LINE_ELEMENT_COLORS.get(line_element(line_name), '0.35')


def resolve_line_ids(wavelength, redshift=None, line_list=None):
    """Return selected line IDs that fall inside the wavelength range."""
    if redshift is None:
        return None

    calib = WavelengthCalibration()
    line_ids = calib.get_line_wavelengths(redshift=redshift, line_list=line_list)
    wave_min = np.nanmin(wavelength)
    wave_max = np.nanmax(wavelength)
    return {
        name: wave for name, wave in line_ids.items()
        if wave_min <= wave <= wave_max
    }


def parse_line_list(line_list):
    """Normalize comma-separated and space-separated line list options."""
    if not line_list:
        return None

    parsed = []
    for item in line_list:
        parsed.extend(part.strip() for part in item.split(',') if part.strip())
    return parsed or None


def draw_line_ids(ax, wavelength, values, line_ids):
    """Draw vertical line IDs on a wavelength plot."""
    if not line_ids:
        return

    wave_min = np.nanmin(wavelength)
    wave_max = np.nanmax(wavelength)
    for line_name, wave in line_ids.items():
        if wave_min <= wave <= wave_max:
            color = line_color(line_name)
            ax.axvline(
                wave, color=color, linestyle=':',
                alpha=0.7, linewidth=1.2
            )
            ax.text(
                wave, 0.98, line_name, rotation=90,
                transform=ax.get_xaxis_transform(),
                verticalalignment='top', horizontalalignment='right',
                fontsize=7, color=color, clip_on=True
            )


def draw_telluric_bands(ax, wavelength):
    """Shade common atmospheric absorption bands in wavelength plots."""
    finite = np.asarray(wavelength)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return

    wave_min = np.nanmin(finite)
    wave_max = np.nanmax(finite)
    label_added = False
    for band_name, (band_min, band_max) in zip(TELLURIC_BAND_LABELS, TELLURIC_BANDS):
        if band_max < wave_min or band_min > wave_max:
            continue

        label = 'Telluric absorption' if not label_added else None
        ax.axvspan(
            band_min, band_max, color='0.5', alpha=0.12,
            linewidth=0, label=label, zorder=0
        )
        ax.text(
            0.5 * (band_min + band_max), 0.98, band_name,
            transform=ax.get_xaxis_transform(),
            ha='center', va='top', fontsize=7,
            color='0.35', alpha=0.8
        )
        label_added = True


def shade_untrusted_wavelengths(ax, wave_min, wave_max, trusted_min,
                                trusted_max):
    """Shade wavelengths outside the trusted calibration range."""
    spans = []
    if wave_min < trusted_min:
        spans.append((wave_min, min(trusted_min, wave_max)))
    if wave_max > trusted_max:
        spans.append((max(trusted_max, wave_min), wave_max))

    label_added = False
    for span_min, span_max in spans:
        if span_max <= span_min:
            continue
        label = 'Outside nominal range' if not label_added else None
        ax.axvspan(
            span_min, span_max, facecolor='0.85', alpha=0.25,
            hatch='///', edgecolor='0.65', linewidth=0.0,
            label=label, zorder=0
        )
        label_added = True


def plot_sensitivity(output_path, wavelength, raw_sensitivity, sensitivity,
                     valid, valid_wave_min, valid_wave_max):
    """Plot raw and smoothed sensitivity."""
    fig, ax = plt.subplots(figsize=(12, 5))
    shade_untrusted_wavelengths(
        ax, np.nanmin(wavelength), np.nanmax(wavelength),
        valid_wave_min, valid_wave_max
    )
    ax.plot(wavelength, raw_sensitivity, '.', color='0.65', markersize=3,
            label='Raw sensitivity', zorder=2)
    ax.plot(wavelength, sensitivity, 'k-', linewidth=1.5,
            label='Smoothed sensitivity', zorder=3)
    if np.any(~valid):
        ax.plot(wavelength[~valid], sensitivity[~valid], '.', color='tab:red',
                markersize=2, alpha=0.4, label='Masked/interpolated', zorder=2)
    ax.set_xlabel('Wavelength (Angstrom)')
    ax.set_ylabel('Sensitivity')
    set_robust_ylim(ax, sensitivity)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best')
    fig.tight_layout()
    fig.savefig(output_path, dpi=120, bbox_inches='tight')
    plt.close(fig)


def plot_standard_comparison(output_path, wavelength, counts_rate,
                             reference_flux, sensitivity, valid,
                             valid_wave_min, valid_wave_max):
    """Plot observed standard counts against the reference/template shape."""
    reference_scaled = reference_flux / sensitivity

    finite = (
        np.isfinite(counts_rate)
        & np.isfinite(reference_scaled)
        & (counts_rate > 0)
        & (reference_scaled > 0)
        & valid
    )
    if np.count_nonzero(finite) >= 2:
        scale = np.nanmedian(counts_rate[finite] / reference_scaled[finite])
    else:
        scale = 1.0
    reference_scaled *= scale

    fig, (ax_compare, ax_ratio) = plt.subplots(
        2, 1, figsize=(12, 7), sharex=True,
        gridspec_kw={'height_ratios': [1.2, 0.8]}
    )

    for ax in (ax_compare, ax_ratio):
        shade_untrusted_wavelengths(
            ax, np.nanmin(wavelength), np.nanmax(wavelength),
            valid_wave_min, valid_wave_max
        )

    ax_compare.plot(wavelength, counts_rate, color='0.25', linewidth=1,
                    label='Observed standard counts/sec', zorder=2)
    ax_compare.plot(wavelength, reference_scaled, color='tab:blue',
                    linewidth=1.5, label='Blackbody/reference scaled to counts',
                    zorder=3)
    if np.any(~valid):
        ax_compare.plot(wavelength[~valid], counts_rate[~valid], '.',
                        color='tab:red', markersize=2, alpha=0.35,
                        label='Masked')
    ax_compare.set_ylabel('Counts/sec')
    set_robust_ylim(ax_compare, counts_rate)
    ax_compare.grid(True, alpha=0.3)
    ax_compare.legend(loc='best')

    ratio = np.full_like(wavelength, np.nan, dtype=float)
    ratio[finite] = counts_rate[finite] / reference_scaled[finite]
    ax_ratio.plot(wavelength, ratio, 'k-', linewidth=1)
    ax_ratio.axhline(1.0, color='0.5', linestyle='--', linewidth=1)
    ax_ratio.set_xlabel('Wavelength (Angstrom)')
    ax_ratio.set_ylabel('Obs / template')
    set_robust_ylim(ax_ratio, ratio, lower=2, upper=98)
    ax_ratio.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(output_path, dpi=120, bbox_inches='tight')
    plt.close(fig)


def plot_calibrated_spectrum(output_path, wavelength, calibrated_flux,
                             target_counts, valid_wave_min, valid_wave_max,
                             line_ids=None):
    """Plot calibrated flux and original counts."""
    fig, (ax_flux, ax_counts) = plt.subplots(
        2, 1, figsize=(12, 7), sharex=True,
        gridspec_kw={'height_ratios': [1.2, 0.8]}
    )
    for ax in (ax_flux, ax_counts):
        draw_telluric_bands(ax, wavelength)
        shade_untrusted_wavelengths(
            ax, np.nanmin(wavelength), np.nanmax(wavelength),
            valid_wave_min, valid_wave_max
        )
    ax_flux.plot(wavelength, calibrated_flux, 'b-', linewidth=1)
    ax_flux.set_ylabel('Calibrated flux')
    set_robust_ylim(ax_flux, calibrated_flux)
    draw_line_ids(ax_flux, wavelength, calibrated_flux, line_ids)
    ax_flux.grid(True, alpha=0.3)

    ax_counts.plot(wavelength, target_counts, color='0.25', linewidth=1)
    ax_counts.set_xlabel('Wavelength (Angstrom)')
    ax_counts.set_ylabel('Original counts')
    set_robust_ylim(ax_counts, target_counts)
    draw_line_ids(ax_counts, wavelength, target_counts, line_ids)
    ax_counts.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(output_path, dpi=120, bbox_inches='tight')
    plt.close(fig)


def calibrate(args):
    """Run flux calibration."""
    if args.valid_wave_min >= args.valid_wave_max:
        raise ValueError("--valid-wave-min must be smaller than --valid-wave-max")

    standard_wave, standard_counts, _ = load_quicklook_spectrum(args.standard)
    standard_exptime = resolve_exptime(
        'standard', args.standard, args.standard_exptime,
        args.standard_fits, spectra_dir=args.spectra_dir
    )
    standard_airmass = None

    if args.reference:
        reference_wave, reference_flux = load_reference_flux(args.reference)
    else:
        reference_wave, reference_flux = load_template_reference(
            args.template_standard,
            wave_min=args.template_wave_min,
            wave_max=args.template_wave_max,
            wave_step=args.template_wave_step,
            teff=args.template_teff,
            v_mag=args.template_vmag,
        )

    sensitivity_wave, raw_sensitivity, sensitivity, valid = build_sensitivity(
        standard_wave, standard_counts, reference_wave, reference_flux,
        standard_exptime=standard_exptime,
        smooth_window=args.smooth_window,
        mask_telluric=not args.no_telluric_mask,
        mask_standard_features=not args.no_standard_feature_mask,
        min_standard_count_frac=args.min_standard_count_frac,
        valid_wave_min=args.valid_wave_min,
        valid_wave_max=args.valid_wave_max
    )

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    standard_name = Path(args.standard).stem.replace('_1d', '')

    reference_interp = interpolate_reference(
        reference_wave, reference_flux, sensitivity_wave
    )
    standard_counts_sorted = sort_by_wavelength(
        standard_wave, standard_counts
    )[1]
    counts_rate = standard_counts_sorted / standard_exptime

    sensitivity_txt_path = output_dir / f"{standard_name}_sensitivity.txt"
    sensitivity_path = output_dir / f"{standard_name}_sensitivity.png"
    comparison_path = output_dir / f"{standard_name}_standard_vs_template.png"
    template_path = output_dir / f"{args.template_standard}_template_reference.txt"

    save_sensitivity(
        sensitivity_txt_path, sensitivity_wave, raw_sensitivity, sensitivity,
        valid, reference_interp, counts_rate
    )
    plot_sensitivity(
        sensitivity_path, sensitivity_wave, raw_sensitivity, sensitivity, valid,
        args.valid_wave_min, args.valid_wave_max
    )
    plot_standard_comparison(
        comparison_path, sensitivity_wave, counts_rate,
        reference_interp, sensitivity, valid,
        args.valid_wave_min, args.valid_wave_max
    )

    if not args.reference and args.save_template:
        np.savetxt(
            template_path, np.column_stack([reference_wave, reference_flux]),
            header="Wavelength(Angstrom)  Template_flux",
            fmt='%.8e', comments='# '
        )
        logger.info("Saved template reference: %s", template_path)

    logger.info("Saved sensitivity table: %s", sensitivity_txt_path)
    logger.info("Saved sensitivity plot: %s", sensitivity_path)
    logger.info("Saved standard comparison plot: %s", comparison_path)

    if args.target:
        target_wave, target_counts, target_sky = load_quicklook_spectrum(args.target)
        target_exptime = resolve_exptime(
            'target', args.target, args.target_exptime,
            args.target_fits, spectra_dir=args.spectra_dir
        )
        target_wave, target_counts, target_sky = sort_by_wavelength(
            target_wave, target_counts, target_sky
        )

        airmass_correction = None
        if not args.no_airmass_correction:
            standard_airmass = resolve_airmass(
                'standard', args.standard, args.standard_airmass,
                args.standard_fits, spectra_dir=args.spectra_dir
            )
            target_airmass = resolve_airmass(
                'target', args.target, args.target_airmass,
                args.target_fits, spectra_dir=args.spectra_dir
            )
            extinction_wave, extinction_mag_table = load_extinction_curve(
                args.extinction_curve
            )
            extinction_mag = extinction_at_wavelength(
                target_wave, extinction_wave, extinction_mag_table
            )
            airmass_correction = differential_airmass_correction(
                target_wave, standard_airmass, target_airmass,
                extinction_wave, extinction_mag_table
            )
            logger.info(
                "Applying differential airmass correction: "
                "standard X=%.6f, target X=%.6f, delta X=%.6f",
                standard_airmass, target_airmass,
                target_airmass - standard_airmass
            )

        calibrated_flux, target_sensitivity = apply_sensitivity(
            target_wave, target_counts, sensitivity_wave, sensitivity,
            target_exptime=target_exptime,
            airmass_correction=airmass_correction
        )
        line_ids = resolve_line_ids(
            target_wave,
            redshift=getattr(args, 'line_redshift', None),
            line_list=parse_line_list(getattr(args, 'line_list', None))
        )

        target_name = Path(args.target).stem.replace('_1d', '')
        spectrum_path = output_dir / f"{target_name}_flux_calibrated.txt"
        plot_path = output_dir / f"{target_name}_flux_calibrated.png"
        airmass_txt_path = output_dir / f"{target_name}_airmass_correction.txt"
        airmass_plot_path = output_dir / f"{target_name}_airmass_correction.png"

        save_calibrated_spectrum(
            spectrum_path, target_wave, calibrated_flux,
            target_counts, target_sky, target_sensitivity,
            airmass_correction=airmass_correction
        )
        plot_calibrated_spectrum(
            plot_path, target_wave, calibrated_flux, target_counts,
            args.valid_wave_min, args.valid_wave_max,
            line_ids=line_ids
        )
        if airmass_correction is not None:
            save_airmass_correction(
                airmass_txt_path, target_wave, airmass_correction,
                extinction_mag, standard_airmass, target_airmass
            )
            plot_airmass_correction(
                airmass_plot_path, target_wave, airmass_correction,
                extinction_mag, standard_airmass, target_airmass,
                args.valid_wave_min, args.valid_wave_max
            )
            logger.info("Saved airmass correction table: %s", airmass_txt_path)
            logger.info("Saved airmass correction plot: %s", airmass_plot_path)
        logger.info("Saved calibrated spectrum: %s", spectrum_path)
        logger.info("Saved calibrated plot: %s", plot_path)


def parse_args():
    parser = argparse.ArgumentParser(
        description='Experimental flux calibration from quicklook 1D spectra'
    )
    parser.add_argument('--standard', required=True,
                        help='Quicklook 1D spectrum of the standard star')
    parser.add_argument('--reference',
                        help='Reference standard-star flux table')
    parser.add_argument('--template-standard', default='HR4468',
                        choices=sorted(STANDARD_TEMPLATES),
                        help='Built-in approximate standard template if --reference is omitted')
    parser.add_argument('--template-teff', type=float,
                        help='Override template effective temperature in K')
    parser.add_argument('--template-vmag', type=float,
                        help='Override template V magnitude')
    parser.add_argument('--template-wave-min', type=float,
                        help='Template wavelength minimum in Angstrom')
    parser.add_argument('--template-wave-max', type=float,
                        help='Template wavelength maximum in Angstrom')
    parser.add_argument('--template-wave-step', type=float,
                        help='Template wavelength step in Angstrom')
    parser.add_argument('--target',
                        help='Quicklook 1D spectrum of the science target')
    parser.add_argument('--output-dir', default='flux_calibrated',
                        help='Output directory')
    parser.add_argument('--standard-fits', nargs='*', default=[],
                        help='FITS file(s) for the standard star; EXPTIME is summed')
    parser.add_argument('--target-fits', nargs='*', default=[],
                        help='FITS file(s) for the target; EXPTIME is summed')
    parser.add_argument('--spectra-dir', default=DEFAULT_SPECTRA_DIR,
                        help='Directory searched for FITS files when *_fits is omitted')
    parser.add_argument('--standard-exptime', type=float,
                        help='Override standard-star exposure time in seconds')
    parser.add_argument('--target-exptime', type=float,
                        help='Override target exposure time in seconds')
    parser.add_argument('--standard-airmass', type=float,
                        help='Override standard-star airmass')
    parser.add_argument('--target-airmass', type=float,
                        help='Override target airmass')
    parser.add_argument('--extinction-curve',
                        help='Two-column wavelength/extinction-mag-per-airmass table')
    parser.add_argument('--no-airmass-correction', action='store_true',
                        help='Disable differential airmass correction')
    parser.add_argument('--smooth-window', type=int, default=51,
                        help='Median smoothing window in pixels')
    parser.add_argument('--no-telluric-mask', action='store_true',
                        help='Do not mask common telluric bands in sensitivity')
    parser.add_argument('--no-standard-feature-mask', action='store_true',
                        help='Do not mask Balmer absorption regions in the standard')
    parser.add_argument('--min-standard-count-frac', type=float, default=0.02,
                        help='Mask standard points below this fraction of the 95th percentile count rate')
    parser.add_argument('--valid-wave-min', type=float,
                        default=DEFAULT_VALID_WAVE_MIN,
                        help='Lower bound of nominally trusted wavelength range in Angstrom')
    parser.add_argument('--valid-wave-max', type=float,
                        default=DEFAULT_VALID_WAVE_MAX,
                        help='Upper bound of nominally trusted wavelength range in Angstrom')
    parser.add_argument('--line-z', dest='line_redshift', type=float,
                        help='Redshift for line ID markers on calibrated target plots')
    parser.add_argument('--line-list', nargs='*',
                        help='Line names or prefix groups for calibrated target plots')
    parser.add_argument('--save-template', action='store_true',
                        help='Save the built-in template reference used')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Verbose logging')
    return parser.parse_args()


def main():
    args = parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    calibrate(args)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

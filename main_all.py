#!/usr/bin/env python3
"""
Run quicklook reduction and experimental flux calibration in one command.

This is a wrapper. The stable quicklook implementation stays in main.py, and
the experimental flux calibration stays in experimental_flux_calibration/.
"""

import argparse
import logging
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from main import (
    create_default_wavelength_solution,
    parse_line_list,
    parse_name_list,
    print_file_mapping,
    print_summary,
    process_all_pairs,
)
from pipeline.preprocessing import PreprocessingPipeline
from pipeline.extraction import SpectralExtraction
from pipeline.visualization import SpectrumVisualizer
from pipeline.utils import setup_logging
from experimental_flux_calibration.flux_calibrate import (
    calibrate,
    set_robust_ylim,
    shade_untrusted_wavelengths,
)


logger = logging.getLogger(__name__)


def find_result_by_object(results, object_name):
    """Return the processing result for object_name."""
    for result in results:
        if result.get('object_name') == object_name:
            return result
    return None


def successful_1d_results(results):
    """Return results with a saved quicklook 1D text spectrum."""
    return [
        result for result in results
        if result.get('stage_1d_plot') and result.get('files', {}).get('1d_txt')
    ]


def build_flux_args(args, standard_txt, target_txt=None):
    """Build an argparse-like namespace for flux_calibrate.calibrate()."""
    return SimpleNamespace(
        standard=str(standard_txt),
        reference=args.flux_reference,
        template_standard=args.template_standard,
        template_teff=args.template_teff,
        template_vmag=args.template_vmag,
        template_wave_min=args.template_wave_min,
        template_wave_max=args.template_wave_max,
        template_wave_step=args.template_wave_step,
        target=str(target_txt) if target_txt is not None else None,
        output_dir=args.flux_output_dir,
        standard_fits=[],
        target_fits=[],
        spectra_dir=args.spectra_dir,
        standard_exptime=args.standard_exptime,
        target_exptime=args.target_exptime,
        standard_airmass=args.standard_airmass,
        target_airmass=args.target_airmass,
        extinction_curve=args.extinction_curve,
        no_airmass_correction=args.no_airmass_correction,
        smooth_window=args.smooth_window,
        no_telluric_mask=args.no_telluric_mask,
        no_standard_feature_mask=args.no_standard_feature_mask,
        min_standard_count_frac=args.min_standard_count_frac,
        valid_wave_min=args.valid_wave_min,
        valid_wave_max=args.valid_wave_max,
        line_redshift=args.z if args.z != 0 or args.line_list else None,
        line_list=parse_line_list(args.line_list),
        save_template=args.save_template,
    )


def target_name_from_txt(path):
    """Return object name from quicklook 1D text path."""
    name = Path(path).stem
    if name.endswith('_1d'):
        name = name[:-3]
    return name


def draw_emission_lines(ax, wavelength, spectrum, emission_lines):
    """Draw emission-line markers using the quicklook visual style."""
    if not emission_lines:
        return

    wave_min = np.nanmin(wavelength)
    wave_max = np.nanmax(wavelength)
    for line_name, wave in emission_lines.items():
        if wave_min <= wave <= wave_max:
            color = SpectrumVisualizer._line_color(line_name)
            ax.axvline(
                x=wave, color=color, linestyle=':',
                alpha=0.7, linewidth=1.4
            )
            SpectrumVisualizer._add_line_label(ax, wave, line_name, color)


def plot_summary_all(target_result, calibrated_txt, output_png,
                     valid_wave_min=4500.0, valid_wave_max=8500.0,
                     dpi=120):
    """Create a summary figure with 2D zoom, quicklook 1D, and calibrated flux."""
    calibrated_txt = Path(calibrated_txt)
    output_png = Path(output_png)
    plot_data = target_result.get('plot_data') or {}

    required = ['image_2d', 'wavelength', 'spectrum_1d']
    missing = [key for key in required if plot_data.get(key) is None]
    if missing:
        raise ValueError(f"Missing quicklook plot data: {', '.join(missing)}")

    image_2d = plot_data['image_2d']
    wavelength = plot_data['wavelength']
    spectrum_1d = plot_data['spectrum_1d']
    sky_spectrum = plot_data.get('sky_spectrum')
    emission_lines = plot_data.get('emission_lines')
    aperture_width = plot_data.get('aperture_width')
    trace_y = target_result.get('trace_y')
    object_name = target_result.get('object_name', target_name_from_txt(calibrated_txt))

    calibrated_data = np.loadtxt(calibrated_txt)
    if calibrated_data.ndim != 2 or calibrated_data.shape[1] < 2:
        raise ValueError(f"Expected at least two columns in {calibrated_txt}")

    calibrated_wavelength = calibrated_data[:, 0]
    calibrated_flux = calibrated_data[:, 1]

    zoom_half_height = 120
    y_center = int(np.round(trace_y))
    y_min = max(0, y_center - zoom_half_height)
    y_max = min(image_2d.shape[0], y_center + zoom_half_height + 1)
    zoom_image = image_2d[y_min:y_max, :]

    fig, (ax_2d, ax_1d, ax_flux) = plt.subplots(
        3, 1, figsize=(14, 12), dpi=dpi,
        gridspec_kw={'height_ratios': [1.05, 1.0, 0.82]},
        constrained_layout=True
    )

    vmin, vmax = SpectrumVisualizer._get_display_limits(zoom_image)
    y_extent_max = y_min + zoom_image.shape[0] - 1
    im = ax_2d.imshow(
        zoom_image, origin='lower', cmap='RdBu_r', aspect='auto',
        vmin=vmin, vmax=vmax, interpolation='nearest',
        extent=(-0.5, zoom_image.shape[1] - 0.5,
                y_min - 0.5, y_extent_max + 0.5)
    )
    SpectrumVisualizer._draw_trace_overlay(
        ax_2d, zoom_image.shape[1], trace_y,
        aperture_width=aperture_width
    )
    ax_2d.set_title(f'Quicklook + Flux Summary - {object_name}')
    ax_2d.set_xlabel('Wavelength direction (pixels)')
    ax_2d.set_ylabel('Spatial direction (pixels)')
    SpectrumVisualizer._add_wavelength_axis(ax_2d, wavelength)
    ax_2d.legend(loc='upper right')
    fig.colorbar(im, ax=ax_2d, label='Intensity', pad=0.01)

    ax_1d.plot(wavelength, spectrum_1d, 'b-', label='Spectrum',
               linewidth=1.4, zorder=2)
    SpectrumVisualizer._draw_atmospheric_bands(ax_1d, wavelength)
    if sky_spectrum is not None:
        ax_1d.plot(wavelength, sky_spectrum, 'r--', alpha=0.5,
                   label='BG', linewidth=1, zorder=2)
    draw_emission_lines(ax_1d, wavelength, spectrum_1d, emission_lines)
    ax_1d.set_xlabel('Wavelength (Angstrom)')
    ax_1d.set_ylabel('Flux (ADU)')
    ax_1d.set_xlim(np.nanmin(wavelength), np.nanmax(wavelength))
    ax_1d.grid(True, alpha=0.3)
    ax_1d.legend(loc='upper right')

    ax_flux.plot(calibrated_wavelength, calibrated_flux,
                 color='tab:blue', linewidth=1.1)
    SpectrumVisualizer._draw_atmospheric_bands(ax_flux, calibrated_wavelength)
    draw_emission_lines(
        ax_flux, calibrated_wavelength, calibrated_flux, emission_lines
    )
    ax_flux.set_xlim(np.nanmin(wavelength), np.nanmax(wavelength))
    shade_untrusted_wavelengths(
        ax_flux, np.nanmin(wavelength), np.nanmax(wavelength),
        valid_wave_min, valid_wave_max
    )
    ax_flux.set_title('Flux Calibrated Spectrum', fontsize=12)
    ax_flux.set_xlabel('Wavelength (Angstrom)')
    ax_flux.set_ylabel('Flux')
    ax_flux.grid(True, alpha=0.3)
    set_robust_ylim(ax_flux, calibrated_flux)

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    logger.info("Saved combined summary: %s", output_png)
    return output_png


def run_quicklook(args):
    """Run the stable quicklook pipeline and return processing results."""
    prep = PreprocessingPipeline(
        spectra_dir=args.spectra_dir,
        pattern_a=args.pattern_a,
        pattern_b=args.pattern_b,
    )
    selected_objects = parse_name_list(args.objects)
    print_file_mapping(prep, object_names=selected_objects)

    extraction = SpectralExtraction(aperture_width=args.aperture)
    visualizer = SpectrumVisualizer(output_dir=args.output_dir)
    calib = create_default_wavelength_solution(npix=1024)

    results = process_all_pairs(
        prep, extraction, visualizer, calib,
        redshift=args.z, verbose=args.verbose,
        line_list=parse_line_list(args.line_list),
        object_names=selected_objects
    )
    print_summary(results)
    return results


def run_flux_calibration(args, results):
    """Run experimental flux calibration using quicklook 1D outputs."""
    if args.no_flux_calibration:
        logger.info("Skipping flux calibration (--no-flux-calibration)")
        return True

    standard_result = find_result_by_object(results, args.standard_object)
    if standard_result is None:
        logger.error("Standard object %s was not processed", args.standard_object)
        return False

    standard_txt = standard_result.get('files', {}).get('1d_txt')
    if not standard_txt:
        logger.error("Standard object %s has no quicklook 1D text output",
                     args.standard_object)
        return False

    candidates = successful_1d_results(results)
    targets = [
        result for result in candidates
        if result.get('object_name') != args.standard_object
    ]
    if args.flux_targets:
        target_names = set(args.flux_targets)
        targets = [
            result for result in targets
            if result.get('object_name') in target_names
        ]

    print("\n" + "=" * 70)
    print("EXPERIMENTAL FLUX CALIBRATION")
    print("=" * 70)
    print(f"Standard: {args.standard_object} ({standard_txt})")

    if not targets:
        print("No science targets selected. Building standard sensitivity only.")
        calibrate(build_flux_args(args, standard_txt, target_txt=None))
        print("=" * 70)
        return True

    success = True
    for target in targets:
        target_name = target.get('object_name')
        target_txt = target.get('files', {}).get('1d_txt')
        print(f"\nCalibrating target: {target_name} ({target_txt})")
        try:
            calibrate(build_flux_args(args, standard_txt, target_txt=target_txt))
            flux_name = target_name_from_txt(target_txt)
            calibrated_txt = (
                Path(args.flux_output_dir) / f"{flux_name}_flux_calibrated.txt"
            )
            if target.get('plot_data') and calibrated_txt.exists():
                output_png = (
                    Path(args.flux_output_dir) / f"{flux_name}_summary_all.png"
                )
                plot_summary_all(
                    target, calibrated_txt, output_png,
                    valid_wave_min=args.valid_wave_min,
                    valid_wave_max=args.valid_wave_max
                )
            else:
                logger.warning(
                    "Skipping summary_all for %s: plot_data=%s calibrated_txt=%s",
                    target_name, bool(target.get('plot_data')), calibrated_txt
                )
        except Exception as exc:
            logger.error("Flux calibration failed for %s: %s", target_name, exc)
            success = False

    print("=" * 70)
    return success


def parse_args():
    parser = argparse.ArgumentParser(
        description='Run KOspec quicklook and experimental flux calibration'
    )

    parser.add_argument('--spectra-dir', default='./spectra',
                        help='Input FITS directory (default: ./spectra)')
    parser.add_argument('--output-dir', default='./quicklook',
                        help='Quicklook output directory (default: ./quicklook)')
    parser.add_argument('--pattern-a', default='_A.fits',
                        help='Pattern for A-position frames')
    parser.add_argument('--pattern-b', default='_B.fits',
                        help='Pattern for B-position frames')
    parser.add_argument('--aperture', type=int, default=10,
                        help='Aperture width in pixels')
    parser.add_argument('--z', type=float, default=0.0,
                        help='Redshift for line marking')
    parser.add_argument('--line-list', nargs='*',
                        help='Line names or prefix groups to mark')
    parser.add_argument('--objects', nargs='*',
                        help='Only process these object names (default: all pairs)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Verbose logging')

    parser.add_argument('--no-flux-calibration', action='store_true',
                        help='Only run the quicklook pipeline')
    parser.add_argument('--standard-object', default='HR4468',
                        help='Object name used as flux standard')
    parser.add_argument('--flux-targets', nargs='*',
                        help='Specific object names to flux-calibrate')
    parser.add_argument('--flux-output-dir', default='flux_calibrated',
                        help='Flux calibration output directory')
    parser.add_argument('--flux-reference',
                        help='External reference standard flux table')
    parser.add_argument('--template-standard', default='HR4468',
                        help='Built-in template standard name')
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
    parser.add_argument('--standard-exptime', type=float,
                        help='Override standard exposure time in seconds')
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
                        help='Median smoothing window for sensitivity')
    parser.add_argument('--no-telluric-mask', action='store_true',
                        help='Do not mask telluric bands for sensitivity')
    parser.add_argument('--no-standard-feature-mask', action='store_true',
                        help='Do not mask Balmer regions in the standard')
    parser.add_argument('--min-standard-count-frac', type=float, default=0.02,
                        help='Mask standard points below this fraction of the 95th percentile')
    parser.add_argument('--valid-wave-min', type=float, default=4500.0,
                        help='Lower bound of trusted wavelength range')
    parser.add_argument('--valid-wave-max', type=float, default=8500.0,
                        help='Upper bound of trusted wavelength range')
    parser.add_argument('--save-template', action='store_true',
                        help='Save the built-in template reference')

    return parser.parse_args()


def main():
    args = parse_args()
    setup_logging(verbose=args.verbose)

    logger.info("Starting combined quicklook + flux calibration")
    logger.info("Input directory: %s", args.spectra_dir)
    logger.info("Quicklook output directory: %s", args.output_dir)
    logger.info("Flux output directory: %s", args.flux_output_dir)

    results = run_quicklook(args)
    quicklook_success = any(result.get('stage_1d_plot') for result in results)
    if not quicklook_success:
        return 1

    flux_success = run_flux_calibration(args, results)
    return 0 if flux_success else 1


if __name__ == '__main__':
    sys.exit(main())

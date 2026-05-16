#!/usr/bin/python3
"""
Main quicklook pipeline script for spectroscopic data
"""

import argparse
import sys
import logging
import subprocess
import shutil
from pathlib import Path
import numpy as np

from pipeline.preprocessing import PreprocessingPipeline
from pipeline.loader import FITSLoader
from pipeline.extraction import SpectralExtraction
from pipeline.calibration import WavelengthCalibration
from pipeline.visualization import SpectrumVisualizer
from pipeline.utils import setup_logging

logger = logging.getLogger(__name__)

FITS_SUFFIXES = {'.fits', '.fit', '.fts'}


def _strip_config_value(value):
    """Return a plain scalar value from a small YAML-like config line."""
    value = value.split('#', 1)[0].strip()
    if not value:
        return None
    if (
        (value.startswith('"') and value.endswith('"'))
        or (value.startswith("'") and value.endswith("'"))
    ):
        value = value[1:-1]
    return value


def load_config_value(config_path, keys, default=None):
    """Load a nested scalar from config.yaml without requiring PyYAML."""
    config_path = Path(config_path)
    if not config_path.exists():
        return default

    stack = []
    try:
        lines = config_path.read_text(encoding='utf-8').splitlines()
    except OSError:
        return default

    for line in lines:
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        if ':' not in line:
            continue

        indent = len(line) - len(line.lstrip(' '))
        level = indent // 2
        key, value = line.strip().split(':', 1)
        stack = stack[:level]
        stack.append(key.strip())

        if stack == list(keys):
            parsed = _strip_config_value(value)
            return default if parsed is None else parsed

    return default


def candidate_fits_files(rawdata_dir, pattern):
    """Return FITS candidates sorted newest first."""
    rawdata_dir = Path(rawdata_dir)
    if not rawdata_dir.exists():
        raise FileNotFoundError(f"Rawdata directory not found: {rawdata_dir}")

    if pattern and pattern != "*.fits":
        candidates = [
            path for path in rawdata_dir.glob(pattern)
            if path.is_file() and path.suffix.lower() in FITS_SUFFIXES
        ]
    else:
        candidates = [
            path for path in rawdata_dir.iterdir()
            if path.is_file() and path.suffix.lower() in FITS_SUFFIXES
        ]

    return sorted(candidates, key=lambda path: path.stat().st_mtime, reverse=True)


def frame_metadata(filepath, loader, pattern_a, pattern_b):
    """Return object name and A/B position metadata for a FITS frame."""
    _, header, success, msg = loader.load(filepath)
    if not success:
        raise ValueError(f"Cannot read FITS header for {filepath}: {msg}")

    original_object = loader.get_header_value(header, 'OBJECT', 'UNKNOWN')
    if original_object == 'UNKNOWN':
        raise ValueError(f"No OBJECT keyword in {filepath}")

    object_name = str(original_object)
    position = None
    if object_name.endswith('_A'):
        object_name = object_name[:-2]
        position = 'A'
    elif object_name.endswith('_B'):
        object_name = object_name[:-2]
        position = 'B'

    if position is None:
        position = str(loader.get_header_value(header, 'POSITION', '')).upper()

    if position not in ['A', 'B']:
        name = filepath.name
        if pattern_a in name or '_A' in name:
            position = 'A'
        elif pattern_b in name or '_B' in name:
            position = 'B'

    if position not in ['A', 'B']:
        raise ValueError(
            f"Cannot determine A/B position for {filepath.name} "
            f"(OBJECT={original_object})"
        )

    return object_name, position


def import_latest_pair(args):
    """Copy the newest raw A/B FITS pair into the import staging directory."""
    rawdata_dir = args.rawdata_dir or load_config_value(
        args.config, ['rawdata', 'dir'], default=None
    )
    if not rawdata_dir:
        raise ValueError(
            "Rawdata directory is not set. Add rawdata.dir to config.yaml "
            "or pass --rawdata-dir."
        )

    pattern = args.latest_pattern or load_config_value(
        args.config, ['rawdata', 'pattern'], default='*.fits'
    )
    import_dir = Path(args.import_dir or load_config_value(
        args.config, ['rawdata', 'import_dir'], default='./data/spectra/latest'
    ))
    import_dir.mkdir(parents=True, exist_ok=True)

    loader = FITSLoader()
    candidates = candidate_fits_files(rawdata_dir, pattern)
    if len(candidates) < 2:
        raise ValueError(f"Need at least two FITS files in {rawdata_dir}")

    latest_two = candidates[:2]
    metadata = [
        frame_metadata(path, loader, args.pattern_a, args.pattern_b)
        for path in latest_two
    ]

    object_names = {item[0] for item in metadata}
    positions = {item[1] for item in metadata}
    if len(object_names) != 1 or positions != {'A', 'B'}:
        details = ', '.join(
            f"{path.name}: OBJECT={object_name}, POS={position}"
            for path, (object_name, position) in zip(latest_two, metadata)
        )
        raise ValueError(f"Latest two FITS files are not an A/B pair: {details}")

    for old_file in import_dir.iterdir():
        if old_file.is_file() and old_file.suffix.lower() in FITS_SUFFIXES:
            old_file.unlink()

    copied = []
    for source in sorted(latest_two):
        destination = import_dir / source.name
        shutil.copy2(source, destination)
        copied.append(destination)

    object_name = next(iter(object_names))
    print("\nImported latest raw A/B pair:", flush=True)
    for source, destination in zip(sorted(latest_two), copied):
        print(f"  {source} -> {destination}", flush=True)

    args.spectra_dir = str(import_dir)
    args.objects = [object_name]
    return copied


def create_default_wavelength_solution(npix):
    """
    Create default wavelength solution
    
    Parameters
    ----------
    npix : int
        Number of wavelength pixels
    
    Returns
    -------
    calib : WavelengthCalibration
        Calibration object
    """
    # Default: 4000 to 9000 Angstrom
    return WavelengthCalibration.create_default_solution(
        npix=npix, wave_start=4000, wave_end=9000
    )


def print_file_mapping(prep_pipeline, object_names=None):
    """Print mapping of files to objects"""
    print("\n" + "="*70)
    print("FILE TO OBJECT MAPPING")
    print("="*70)
    
    pairs = [
        pair for pair in prep_pipeline.find_ab_pairs()
        if prep_pipeline._object_is_selected(pair[2], object_names)
    ]
    if not pairs:
        print("No A-B pairs found.")
        return
    
    print(f"Found {len(pairs)} object(s):")
    for filepath_a, filepath_b, object_name in pairs:
        print(f"\n  Object: {object_name}")
        print(f"    A-frame: {filepath_a.name}")
        print(f"    B-frame: {filepath_b.name}")
    
    print("="*70)


def print_processing_progress(object_name, stage, status, details=""):
    """Print processing progress with nice formatting"""
    status_icons = {
        'start': '🔄',
        'ab_diff': '📊',
        '2d_plot': '🖼️',
        'extract': '📈',
        '1d_plot': '📊',
        'complete': '✅',
        'error': '❌',
        'warning': '⚠️'
    }
    
    icon = status_icons.get(status, '•')
    
    if status == 'start':
        print(f"\n{icon} Processing {object_name}...")
    elif status == 'complete':
        print(f"{icon} {object_name}: Complete")
    elif status == 'error':
        print(f"{icon} {object_name}: ERROR - {details}")
    elif status == 'warning':
        print(f"{icon} {object_name}: WARNING - {details}")
    else:
        stage_names = {
            'ab_diff': 'A-B difference',
            '2d_plot': '2D spectrum plot',
            'extract': '1D spectrum extraction',
            '1d_plot': '1D spectrum plot'
        }
        stage_name = stage_names.get(stage, stage)
        if status == 'success':
            print(f"  {icon} {stage_name}: OK")
        elif status == 'failed':
            print(f"  {icon} {stage_name}: FAILED - {details}")
        else:
            print(f"  {icon} {stage_name}: {status}")


def parse_line_list(line_list):
    """Normalize comma-separated and space-separated line list options."""
    if not line_list:
        return None

    parsed = []
    for item in line_list:
        parsed.extend(part.strip() for part in item.split(',') if part.strip())
    return parsed or None


def parse_name_list(name_list):
    """Normalize comma-separated and space-separated object name options."""
    if not name_list:
        return None

    parsed = []
    for item in name_list:
        parsed.extend(part.strip() for part in item.split(',') if part.strip())
    return parsed or None


def process_all_pairs(prep_pipeline, extraction, visualizer, calib,
                      redshift=0.0, verbose=False, line_list=None,
                      object_names=None):
    """
    Process all A-B pairs through the full pipeline
    
    Parameters
    ----------
    prep_pipeline : PreprocessingPipeline
        Preprocessing pipeline
    extraction : SpectralExtraction
        Extraction pipeline
    visualizer : SpectrumVisualizer
        Visualization tool
    calib : WavelengthCalibration
        Wavelength calibration
    redshift : float
        Redshift for emission lines
    verbose : bool
        Verbose logging
    line_list : list or None
        Line names or prefix groups to display
    object_names : list or None
        Object names to process. If None, process all A-B pairs.
    
    Returns
    -------
    results : list of dict
        List of processing results
    """
    results = []
    
    for pair_data in prep_pipeline.process_all_pairs(object_names=object_names):
        object_name = pair_data['object_name']
        data_ab = pair_data['data']
        headers = pair_data['headers']
        success_ab = pair_data['success']
        error_ab = pair_data['error']
        
        # Show processing start
        print_processing_progress(object_name, 'start', 'start')
        
        result = {
            'object_name': object_name,
            'stage_ab_diff': success_ab,
            'error_ab_diff': error_ab,
            'stage_2d_plot': False,
            'stage_1d_extract': False,
            'stage_1d_plot': False,
            'files': {}
        }
        
        # If A-B failed, skip but record
        if not success_ab:
            print_processing_progress(object_name, 'ab_diff', 'error', error_ab)
            results.append(result)
            continue
        
        print_processing_progress(object_name, 'ab_diff', 'success')
        
        # Step 1: Plot 2D spectrum (even if 1D extraction fails later)
        try:
            png_2d = visualizer.plot_2d_spectrum(
                data_ab, object_name, trace_y=None, aperture_width=None
            )
            if png_2d:
                result['stage_2d_plot'] = True
                result['files']['2d_png'] = str(png_2d)
                print_processing_progress(object_name, '2d_plot', 'success')
            else:
                print_processing_progress(object_name, '2d_plot', 'failed', 'Plot creation failed')
        except Exception as e:
            print_processing_progress(object_name, '2d_plot', 'error', str(e))
        
        # Step 2: Extract 1D spectrum
        spectrum_1d, sky_spectrum, trace_y, success_extract, error_extract = \
            extraction.extract(data_ab, object_name=object_name, 
                             positive_peak_only=True)
        
        if not success_extract:
            print_processing_progress(object_name, 'extract', 'error', error_extract)
            results.append(result)
            continue
        
        result['stage_1d_extract'] = True
        result['trace_y'] = trace_y
        print_processing_progress(object_name, 'extract', 'success', f"Trace at y={trace_y:.1f}")
        
        # Step 3: Create wavelength array
        try:
            # Get spectrum length
            spec_len = len(spectrum_1d)
            
            # Update calibration if needed
            if 'npix' not in calib.wave_solution or calib.wave_solution['npix'] != spec_len:
                calib_local = create_default_wavelength_solution(spec_len)
            else:
                calib_local = calib
            
            wavelength = calib_local.pixel_to_wavelength(
                np.arange(spec_len)
            )
            
            # Get emission lines if redshift specified
            emission_lines = None
            if redshift is not None and (redshift != 0 or line_list):
                emission_lines = calib_local.get_line_wavelengths(
                    redshift=redshift, line_list=line_list
                )
            
            # Step 4: Save 1D spectrum as text
            txt_path = visualizer.save_1d_spectrum_txt(
                wavelength, spectrum_1d, object_name, 
                sky_spectrum=sky_spectrum,
                comments=f"# Object: {object_name}, Redshift: {redshift}, Trace: {trace_y:.1f}"
            )
            if txt_path:
                result['files']['1d_txt'] = str(txt_path)
            
            # Step 5: Plot 1D spectrum
            png_1d = visualizer.plot_1d_spectrum(
                wavelength, spectrum_1d, object_name,
                sky_spectrum=sky_spectrum,
                emission_lines=emission_lines
            )
            if png_1d:
                result['stage_1d_plot'] = True
                result['files']['1d_png'] = str(png_1d)
                print_processing_progress(object_name, '1d_plot', 'success')

            # Step 6: Create combined 2D zoom + 1D summary
            summary_png = visualizer.plot_summary(
                data_ab, wavelength, spectrum_1d, object_name,
                trace_y=trace_y,
                aperture_width=extraction.aperture_width,
                sky_spectrum=sky_spectrum,
                emission_lines=emission_lines
            )
            if summary_png:
                result['files']['summary_png'] = str(summary_png)

            result['plot_data'] = {
                'image_2d': data_ab,
                'wavelength': wavelength,
                'spectrum_1d': spectrum_1d,
                'sky_spectrum': sky_spectrum,
                'emission_lines': emission_lines,
                'aperture_width': extraction.aperture_width,
            }
            
            # Update 2D plot with trace overlay
            png_2d_updated = visualizer.plot_2d_spectrum(
                data_ab, object_name, trace_y=trace_y,
                aperture_width=extraction.aperture_width,
                wavelength=wavelength
            )
            if png_2d_updated:
                result['files']['2d_png'] = str(png_2d_updated)
                png_2d_zoom = visualizer.output_dir / f"{object_name}_2d_zoom.png"
                if png_2d_zoom.exists():
                    result['files']['2d_zoom_png'] = str(png_2d_zoom)
            
            print_processing_progress(object_name, 'complete', 'complete')
        
        except Exception as e:
            print_processing_progress(object_name, 'complete', 'error', str(e))
            result['error_pipeline'] = str(e)
        
        results.append(result)
    
    return results


def print_summary(results):
    """Print processing summary"""
    print("\n" + "="*70)
    print("QUICKLOOK PIPELINE SUMMARY")
    print("="*70)
    
    total = len(results)
    successful = sum(1 for r in results if r['stage_1d_plot'])
    partial = sum(1 for r in results if r['stage_2d_plot'] and not r['stage_1d_plot'])
    failed = sum(1 for r in results if not r['stage_2d_plot'])
    
    print(f"📊 Total objects: {total}")
    print(f"✅ Fully processed (1D extracted): {successful}")
    print(f"🖼️  Partial (2D only): {partial}")
    print(f"❌ Failed: {failed}")
    
    if successful > 0:
        success_rate = successful / total * 100
        print(f"📈 Success rate: {success_rate:.1f}%")
    
    print("\n📋 Detailed results:")
    for r in results:
        object_name = r['object_name']
        if r['stage_1d_plot']:
            status = "✅"
            details = "Complete (1D + 2D)"
        elif r['stage_2d_plot']:
            status = "🖼️"
            details = "Partial (2D only)"
        else:
            status = "❌"
            details = "Failed"
        
        print(f"  {status} {object_name}: {details}")
        
        # Show file information
        files = r.get('files', {})
        if files:
            file_list = []
            if '2d_png' in files:
                file_list.append("2D plot")
            if '2d_zoom_png' in files:
                file_list.append("2D zoom")
            if 'summary_png' in files:
                file_list.append("Summary")
            if '1d_png' in files:
                file_list.append("1D plot")
            if '1d_txt' in files:
                file_list.append("1D data")
            if file_list:
                print(f"      Files: {', '.join(file_list)}")
        
        # Show errors
        if r.get('error_ab_diff'):
            print(f"      ⚠️  A-B Error: {r['error_ab_diff']}")
        if r.get('error_pipeline'):
            print(f"      ⚠️  Pipeline Error: {r['error_pipeline']}")
    
    print("="*70)


def open_summary_pngs(results):
    """Open generated quicklook summary PNG files with the OS open command."""
    summary_paths = [
        Path(r['files']['summary_png'])
        for r in results
        if r.get('files', {}).get('summary_png')
    ]

    if not summary_paths:
        print("\nNo summary PNG files to open.")
        return

    print("\nOpening summary PNG file(s):")
    for summary_path in summary_paths:
        print(f"  {summary_path}")
        try:
            subprocess.Popen(
                ['open', str(summary_path)],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        except FileNotFoundError:
            print("  WARNING: 'open' command was not found on this PC.")
            break


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description='Quicklook pipeline for spectroscopic data'
    )
    parser.add_argument('--spectra-dir', default='./data/spectra',
                       help='Input FITS directory (default: ./data/spectra)')
    parser.add_argument('--output-dir', default='./outputs/quicklook',
                       help='Output directory (default: ./outputs/quicklook)')
    parser.add_argument('--pattern-a', default='_A.fits',
                       help='Pattern for A-position frames (default: _A.fits)')
    parser.add_argument('--pattern-b', default='_B.fits',
                       help='Pattern for B-position frames (default: _B.fits)')
    parser.add_argument('--aperture', type=int, default=10,
                       help='Aperture width (pixels, default: 10)')
    parser.add_argument('--z', type=float, default=0.0,
                       help='Redshift for emission line marking (default: 0)')
    parser.add_argument('--line-list', nargs='*',
                       help='Line names or prefix groups to mark')
    parser.add_argument('--objects', nargs='*',
                       help='Only process these object names (default: all pairs)')
    parser.add_argument('--config', default='config.yaml',
                       help='Local config file (default: config.yaml)')
    parser.add_argument('--import-latest', action='store_true',
                       help='Copy newest raw A/B FITS pair before running QL')
    parser.add_argument('--rawdata-dir',
                       help='Raw FITS directory for --import-latest')
    parser.add_argument('--latest-pattern',
                       help='Glob pattern for raw FITS search (default from config or *.fits)')
    parser.add_argument('--import-dir',
                       help='Staging directory for imported latest pair (default from config or ./data/spectra/latest)')
    parser.add_argument('--no-open-summary', action='store_true',
                       help='Do not open generated summary PNG files at the end')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Verbose logging')
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(verbose=args.verbose)
    
    logger.info(f"Starting quicklook pipeline")
    logger.info(f"Redshift: {args.z}")

    if args.import_latest:
        try:
            import_latest_pair(args)
        except Exception as exc:
            logger.error(f"Failed to import latest raw A/B pair: {exc}")
            return 1

    logger.info(f"Input directory: {args.spectra_dir}")
    logger.info(f"Output directory: {args.output_dir}")
    
    # Initialize pipelines
    prep = PreprocessingPipeline(spectra_dir=args.spectra_dir,
                                pattern_a=args.pattern_a,
                                pattern_b=args.pattern_b)
    
    # Show file mapping
    selected_objects = parse_name_list(args.objects)
    print_file_mapping(prep, object_names=selected_objects)
    
    extraction = SpectralExtraction(aperture_width=args.aperture)
    visualizer = SpectrumVisualizer(output_dir=args.output_dir)
    calib = create_default_wavelength_solution(npix=1024)
    
    # Process
    results = process_all_pairs(
        prep, extraction, visualizer, calib,
        redshift=args.z, verbose=args.verbose,
        line_list=parse_line_list(args.line_list),
        object_names=selected_objects
    )
    
    # Summary
    print_summary(results)

    if not args.no_open_summary:
        open_summary_pngs(results)
    
    # Return exit code
    successful = sum(1 for r in results if r['stage_1d_plot'])
    if successful > 0:
        return 0
    else:
        return 1


if __name__ == '__main__':
    sys.exit(main())

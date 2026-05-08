#!/usr/bin/env python3
"""
Interactively choose atomic-line markers for a saved quicklook 1D spectrum.
"""

import argparse
from pathlib import Path

import numpy as np

from pipeline.calibration import DEFAULT_LINES, WavelengthCalibration


LINE_ELEMENT_COLORS = {
    'H': 'crimson',
    'He': 'darkorange',
    'Na': 'goldenrod',
    'O': 'forestgreen',
    'Ca': 'purple',
    'S': 'teal',
    'N': 'royalblue',
}

ATMOSPHERIC_ABSORPTION_BANDS = [
    ('O2 B', 6860.0, 6950.0),
    ('H2O', 7160.0, 7340.0),
    ('O2 A', 7590.0, 7700.0),
    ('H2O', 8120.0, 8400.0),
    ('H2O', 8900.0, 9800.0),
]


def line_element(line_name):
    clean_name = line_name.strip('[]')
    if clean_name.startswith('He'):
        return 'He'
    if clean_name.startswith('H'):
        return 'H'
    return clean_name.split('_', 1)[0]


def line_color(line_name):
    return LINE_ELEMENT_COLORS.get(line_element(line_name), '0.35')


def draw_atmospheric_bands(ax, wavelength):
    finite = wavelength[np.isfinite(wavelength)]
    if finite.size == 0:
        return

    wave_min = np.nanmin(finite)
    wave_max = np.nanmax(finite)
    label_added = False
    for band_name, band_min, band_max in ATMOSPHERIC_ABSORPTION_BANDS:
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


def parse_line_args(line_args):
    if not line_args:
        return list(DEFAULT_LINES)

    lines = []
    for item in line_args:
        lines.extend(part.strip() for part in item.split(',') if part.strip())
    return lines


def load_spectrum(path):
    data = np.loadtxt(path)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Expected at least two columns in {path}")

    wavelength = data[:, 0]
    flux = data[:, 1]
    sky = data[:, 2] if data.shape[1] >= 3 else None
    return wavelength, flux, sky


def visible_lines(wavelength, redshift, line_list):
    calib = WavelengthCalibration()
    lines = calib.get_line_wavelengths(redshift=redshift, line_list=line_list)
    wave_min = np.nanmin(wavelength)
    wave_max = np.nanmax(wavelength)
    return {
        name: wave for name, wave in lines.items()
        if wave_min <= wave <= wave_max
    }


def set_all(active, labels, check):
    for index, label in enumerate(labels):
        if check.get_status()[index] != active:
            check.set_active(index)


def interactive_plot(args):
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Button, CheckButtons

    spectrum_path = Path(args.spectrum)
    wavelength, flux, sky = load_spectrum(spectrum_path)
    line_list = parse_line_args(args.lines)
    lines = visible_lines(wavelength, args.z, line_list)

    if not lines:
        raise ValueError("No selected lines fall inside the spectrum wavelength range.")

    fig = plt.figure(figsize=(14, 7))
    ax = fig.add_axes([0.08, 0.12, 0.68, 0.78])
    check_ax = fig.add_axes([0.79, 0.20, 0.18, 0.68])
    all_ax = fig.add_axes([0.79, 0.12, 0.08, 0.05])
    none_ax = fig.add_axes([0.89, 0.12, 0.08, 0.05])

    draw_atmospheric_bands(ax, wavelength)
    ax.plot(wavelength, flux, 'b-', label='Spectrum', linewidth=1.3, zorder=2)
    if sky is not None:
        ax.plot(wavelength, sky, 'r--', alpha=0.5, label='Sky',
                linewidth=1, zorder=2)

    artists = {}
    for name, wave in lines.items():
        color = line_color(name)
        marker = ax.axvline(
            wave, color=color, linestyle=':', alpha=0.75,
            linewidth=1.4, visible=True
        )
        label = ax.text(
            wave, 0.98, name, rotation=90,
            transform=ax.get_xaxis_transform(),
            verticalalignment='top', horizontalalignment='right',
            fontsize=8, color=color, visible=True, clip_on=True
        )
        artists[name] = (marker, label)

    labels = list(lines.keys())
    check = CheckButtons(check_ax, labels, [True] * len(labels))

    def toggle(label):
        for artist in artists[label]:
            artist.set_visible(not artist.get_visible())
        fig.canvas.draw_idle()

    check.on_clicked(toggle)

    all_button = Button(all_ax, 'All')
    none_button = Button(none_ax, 'None')
    all_button.on_clicked(lambda event: set_all(True, labels, check))
    none_button.on_clicked(lambda event: set_all(False, labels, check))

    object_name = spectrum_path.stem.removesuffix('_1d')
    ax.set_title(f'Interactive Line Selection - {object_name}')
    ax.set_xlabel('Wavelength (Angstrom)')
    ax.set_ylabel('Flux (ADU)')
    ax.set_xlim(np.nanmin(wavelength), np.nanmax(wavelength))
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right')

    plt.show()


def parse_args():
    parser = argparse.ArgumentParser(
        description='Interactively toggle line markers on a quicklook 1D spectrum'
    )
    parser.add_argument('spectrum', help='Quicklook 1D text spectrum')
    parser.add_argument('--z', type=float, default=0.0,
                        help='Redshift for line marking')
    parser.add_argument(
        '--lines', nargs='*',
        help='Line names or prefix groups, comma-separated or space-separated'
    )
    return parser.parse_args()


def main():
    interactive_plot(parse_args())


if __name__ == '__main__':
    main()

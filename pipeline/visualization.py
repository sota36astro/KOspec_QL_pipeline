"""
Visualization module for spectra
"""

import numpy as np
import logging
from pathlib import Path
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)


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

ATMOSPHERIC_ABSORPTION_BANDS = [
    ('O2 B', 6860.0, 6950.0),
    ('H2O', 7160.0, 7340.0),
    ('O2 A', 7590.0, 7700.0),
    ('H2O', 8120.0, 8400.0),
    ('H2O', 8900.0, 9800.0),
]


class SpectrumVisualizer:
    """Create visualizations of 2D and 1D spectra"""
    
    def __init__(self, output_dir="./quicklook", dpi=100):
        """
        Initialize visualizer
        
        Parameters
        ----------
        output_dir : str
            Output directory for plots
        dpi : int
            DPI for PNG output
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.dpi = dpi

    @staticmethod
    def _get_display_limits(image_2d, vmin=None, vmax=None):
        """Calculate robust display limits for a 2D image."""
        finite_data = image_2d[np.isfinite(image_2d)]
        if finite_data.size == 0:
            return 0, 1

        if vmin is None or vmax is None:
            if np.any(finite_data < 0):
                vmax_auto = np.abs(np.percentile(finite_data, [1, 99])).max()
                if vmax_auto == 0:
                    vmax_auto = 1
                vmin_auto = -vmax_auto
            else:
                vmin_auto = np.percentile(finite_data, 1)
                vmax_auto = np.percentile(finite_data, 99)
                if vmin_auto == vmax_auto:
                    vmax_auto = vmin_auto + 1

            if vmin is None:
                vmin = vmin_auto
            if vmax is None:
                vmax = vmax_auto

        return vmin, vmax

    @staticmethod
    def _draw_trace_overlay(ax, image_width, trace_y, aperture_width=None):
        """Draw trace and aperture boundaries without hiding the spectrum."""
        ax.axhline(
            y=trace_y, color='yellow', linestyle='--',
            linewidth=1.0, alpha=0.9, label=f'Trace: {trace_y:.1f}'
        )

        if aperture_width is not None:
            half_w = aperture_width / 2
            y_low = trace_y - half_w
            y_high = trace_y + half_w
            ax.hlines(
                [y_low, y_high], xmin=-0.5, xmax=image_width - 0.5,
                colors='lime', linestyles='-', linewidth=1.0, alpha=0.9,
                label=f'Aperture: {aperture_width}px'
            )

    @staticmethod
    def _add_wavelength_axis(ax, wavelength):
        """Add a top wavelength axis while keeping the image x-axis in pixels."""
        if wavelength is None:
            return

        wavelength = np.asarray(wavelength)
        finite = np.isfinite(wavelength)
        if wavelength.ndim != 1 or np.count_nonzero(finite) < 2:
            return

        pixels = np.arange(len(wavelength))[finite]
        waves = wavelength[finite]
        order = np.argsort(pixels)
        pixels = pixels[order]
        waves = waves[order]

        def pixel_to_wave(pixel):
            return np.interp(pixel, pixels, waves)

        def wave_to_pixel(wave):
            if waves[0] <= waves[-1]:
                return np.interp(wave, waves, pixels)
            return np.interp(wave, waves[::-1], pixels[::-1])

        secax = ax.secondary_xaxis('top', functions=(pixel_to_wave, wave_to_pixel))
        secax.set_xlabel('Wavelength (Angstrom)')
        tick_min = np.ceil(np.nanmin(waves) / 1000) * 1000
        tick_max = np.floor(np.nanmax(waves) / 1000) * 1000
        ticks = np.arange(tick_min, tick_max + 1, 1000)
        if ticks.size:
            secax.set_xticks(ticks)

    def _save_2d_panel(self, image_2d, object_name, output_path, title,
                       trace_y=None, aperture_width=None, cmap='RdBu_r',
                       vmin=None, vmax=None, y_offset=0, wavelength=None):
        """Save one 2D spectrum panel."""
        fig, ax = plt.subplots(figsize=(12, 6), dpi=self.dpi)

        vmin, vmax = self._get_display_limits(image_2d, vmin=vmin, vmax=vmax)
        y_max = y_offset + image_2d.shape[0] - 1
        im = ax.imshow(
            image_2d, origin='lower', cmap=cmap, aspect='auto',
            vmin=vmin, vmax=vmax, interpolation='nearest',
            extent=(-0.5, image_2d.shape[1] - 0.5, y_offset - 0.5, y_max + 0.5)
        )

        if trace_y is not None:
            self._draw_trace_overlay(
                ax, image_2d.shape[1], trace_y,
                aperture_width=aperture_width
            )
            ax.legend(loc='upper right')

        ax.set_xlabel('Wavelength direction (pixels)')
        ax.set_ylabel('Spatial direction (pixels)')
        ax.set_title(title)
        self._add_wavelength_axis(ax, wavelength)
        plt.colorbar(im, ax=ax, label='Intensity')

        plt.tight_layout()
        fig.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        plt.close(fig)

        logger.info(f"Saved 2D spectrum plot: {output_path}")
        return output_path

    @staticmethod
    def _line_element(line_name):
        """Return the element key used to color an emission/absorption line."""
        clean_name = line_name.strip('[]')
        if clean_name.startswith('He'):
            return 'He'
        if clean_name.startswith('H'):
            return 'H'
        return clean_name.split('_', 1)[0]

    @classmethod
    def _line_color(cls, line_name):
        element = cls._line_element(line_name)
        return LINE_ELEMENT_COLORS.get(element, '0.35')

    @staticmethod
    def _add_line_label(ax, wave, line_name, color, fontsize=8):
        """Place a line ID label at the top of the panel, away from the spectrum."""
        ax.text(
            wave, 0.98, line_name, rotation=90,
            transform=ax.get_xaxis_transform(),
            verticalalignment='top', horizontalalignment='right',
            fontsize=fontsize, color=color, clip_on=True
        )

    @staticmethod
    def _draw_atmospheric_bands(ax, wavelength):
        """Shade common atmospheric absorption bands in wavelength plots."""
        if wavelength is None:
            return

        wavelength = np.asarray(wavelength)
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
            band_center = 0.5 * (band_min + band_max)
            ax.text(
                band_center, 0.98, band_name,
                transform=ax.get_xaxis_transform(),
                ha='center', va='top', fontsize=7,
                color='0.35', alpha=0.8
            )
            label_added = True
    
    def plot_2d_spectrum(self, image_2d, object_name, trace_y=None, 
                        aperture_width=None, cmap='RdBu_r', 
                        vmin=None, vmax=None, wavelength=None):
        """
        Create 2D spectrum visualization
        
        Parameters
        ----------
        image_2d : 2D array
            2D spectrum image
        object_name : str
            Object name for filename/title
        trace_y : float or None
            Trace position to overlay
        aperture_width : int or None
            Aperture width to visualize
        cmap : str
            Colormap name
        vmin, vmax : float or None
            Color scale limits
        wavelength : 1D array or None
            Wavelength array in Angstrom for the top x-axis
        
        Returns
        -------
        filepath : Path or None
            Path to saved PNG, or None if failed
        """
        if image_2d is None:
            return None
        
        try:
            output_path = self.output_dir / f"{object_name}_2d.png"
            self._save_2d_panel(
                image_2d, object_name, output_path,
                title=f'2D Spectrum - {object_name}',
                trace_y=trace_y, aperture_width=aperture_width,
                cmap=cmap, vmin=vmin, vmax=vmax, y_offset=0,
                wavelength=wavelength
            )

            if trace_y is not None:
                zoom_half_height = 120
                y_center = int(np.round(trace_y))
                y_min = max(0, y_center - zoom_half_height)
                y_max = min(image_2d.shape[0], y_center + zoom_half_height + 1)
                zoom_image = image_2d[y_min:y_max, :]
                zoom_path = self.output_dir / f"{object_name}_2d_zoom.png"
                self._save_2d_panel(
                    zoom_image, object_name, zoom_path,
                    title=f'2D Spectrum Zoom - {object_name}',
                    trace_y=trace_y, aperture_width=aperture_width,
                    cmap=cmap, vmin=None, vmax=None, y_offset=y_min,
                    wavelength=wavelength
                )

            return output_path
        
        except Exception as e:
            logger.error(f"Failed to plot 2D spectrum: {str(e)}")
            return None
    
    def plot_1d_spectrum(self, wavelength, spectrum_1d, object_name,
                        sky_spectrum=None, emission_lines=None, 
                        y_log=False):
        """
        Create 1D spectrum visualization
        
        Parameters
        ----------
        wavelength : 1D array
            Wavelength array (Angstrom)
        spectrum_1d : 1D array
            1D spectrum
        object_name : str
            Object name
        sky_spectrum : 1D array or None
            Sky spectrum to overlay
        emission_lines : dict or None
            Emission lines to mark {name: wavelength}
        y_log : bool
            Use log scale for y-axis
        
        Returns
        -------
        filepath : Path or None
            Path to saved PNG
        """
        if spectrum_1d is None:
            return None
        
        try:
            fig, ax = plt.subplots(figsize=(14, 6), dpi=self.dpi)

            self._draw_atmospheric_bands(ax, wavelength)
            
            # Plot spectrum
            ax.plot(wavelength, spectrum_1d, 'b-', label='Spectrum', 
                   linewidth=1.5, zorder=2)
            
            # Plot sky if available
            if sky_spectrum is not None:
                ax.plot(wavelength, sky_spectrum, 'r--', alpha=0.5,
                       label='BG', linewidth=1, zorder=2)
            
            # Mark emission lines
            if emission_lines:
                for line_name, wave in emission_lines.items():
                    # Check if line is in wavelength range
                    if np.min(wavelength) <= wave <= np.max(wavelength):
                        color = self._line_color(line_name)
                        ax.axvline(x=wave, color=color, linestyle=':', 
                                  alpha=0.7, linewidth=1.5)
                        
                        self._add_line_label(ax, wave, line_name, color)
            
            # Labels and formatting
            ax.set_xlabel('Wavelength (Angstrom)')
            ax.set_ylabel('Flux (ADU)')
            ax.set_title(f'1D Spectrum - {object_name}')
            ax.legend(loc='upper right')
            ax.grid(True, alpha=0.3)
            
            if y_log and np.all(spectrum_1d > 0):
                ax.set_yscale('log')
            
            plt.tight_layout()
            
            # Save
            output_path = self.output_dir / f"{object_name}_1d.png"
            fig.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
            plt.close(fig)
            
            logger.info(f"Saved 1D spectrum plot: {output_path}")
            return output_path
        
        except Exception as e:
            logger.error(f"Failed to plot 1D spectrum: {str(e)}")
            return None

    def plot_summary(self, image_2d, wavelength, spectrum_1d, object_name,
                     trace_y, aperture_width=None, sky_spectrum=None,
                     emission_lines=None, cmap='RdBu_r'):
        """
        Create a combined quicklook summary with 2D zoom and 1D spectrum.

        Parameters
        ----------
        image_2d : 2D array
            A-B 2D spectrum image
        wavelength : 1D array
            Wavelength array (Angstrom)
        spectrum_1d : 1D array
            Sky-subtracted 1D spectrum
        object_name : str
            Object name
        trace_y : float
            Trace position
        aperture_width : int or None
            Aperture width to visualize
        sky_spectrum : 1D array or None
            Sky spectrum to overlay
        emission_lines : dict or None
            Emission lines to mark {name: wavelength}
        cmap : str
            Colormap name for 2D panel

        Returns
        -------
        filepath : Path or None
            Path to saved PNG
        """
        if image_2d is None or spectrum_1d is None or trace_y is None:
            return None

        try:
            zoom_half_height = 120
            y_center = int(np.round(trace_y))
            y_min = max(0, y_center - zoom_half_height)
            y_max = min(image_2d.shape[0], y_center + zoom_half_height + 1)
            zoom_image = image_2d[y_min:y_max, :]

            fig, (ax_2d, ax_1d) = plt.subplots(
                2, 1, figsize=(14, 9), dpi=self.dpi,
                gridspec_kw={'height_ratios': [1.05, 1.0]},
                constrained_layout=True
            )

            vmin, vmax = self._get_display_limits(zoom_image)
            y_extent_max = y_min + zoom_image.shape[0] - 1
            im = ax_2d.imshow(
                zoom_image, origin='lower', cmap=cmap, aspect='auto',
                vmin=vmin, vmax=vmax, interpolation='nearest',
                extent=(-0.5, zoom_image.shape[1] - 0.5,
                        y_min - 0.5, y_extent_max + 0.5)
            )
            self._draw_trace_overlay(
                ax_2d, zoom_image.shape[1], trace_y,
                aperture_width=aperture_width
            )
            ax_2d.set_title(f'Quicklook Summary - {object_name}')
            ax_2d.set_xlabel('Wavelength direction (pixels)')
            ax_2d.set_ylabel('Spatial direction (pixels)')
            self._add_wavelength_axis(ax_2d, wavelength)
            ax_2d.legend(loc='upper right')
            fig.colorbar(im, ax=ax_2d, label='Intensity', pad=0.01)

            ax_1d.plot(wavelength, spectrum_1d, 'b-', label='Spectrum',
                       linewidth=1.4, zorder=2)
            self._draw_atmospheric_bands(ax_1d, wavelength)
            if sky_spectrum is not None:
                ax_1d.plot(wavelength, sky_spectrum, 'r--', alpha=0.5,
                           label='BG', linewidth=1, zorder=2)

            if emission_lines:
                wave_min = np.nanmin(wavelength)
                wave_max = np.nanmax(wavelength)
                for line_name, wave in emission_lines.items():
                    if wave_min <= wave <= wave_max:
                        color = self._line_color(line_name)
                        ax_1d.axvline(
                            x=wave, color=color, linestyle=':',
                            alpha=0.7, linewidth=1.4
                        )
                        self._add_line_label(ax_1d, wave, line_name, color)

            ax_1d.set_xlabel('Wavelength (Angstrom)')
            ax_1d.set_ylabel('Flux (ADU)')
            ax_1d.grid(True, alpha=0.3)
            ax_1d.legend(loc='upper right')

            output_path = self.output_dir / f"{object_name}_summary.png"
            fig.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
            plt.close(fig)

            logger.info(f"Saved summary plot: {output_path}")
            return output_path

        except Exception as e:
            logger.error(f"Failed to plot summary: {str(e)}")
            return None
    
    def save_1d_spectrum_txt(self, wavelength, spectrum_1d, object_name,
                            sky_spectrum=None, comments=None):
        """
        Save 1D spectrum as text file
        
        Parameters
        ----------
        wavelength : 1D array
            Wavelength array (Angstrom)
        spectrum_1d : 1D array
            1D spectrum
        object_name : str
            Object name
        sky_spectrum : 1D array or None
            Sky spectrum
        comments : str or None
            Additional comments
        
        Returns
        -------
        filepath : Path or None
            Path to saved file
        """
        if spectrum_1d is None:
            return None
        
        try:
            output_path = self.output_dir / f"{object_name}_1d.txt"
            
            # Prepare data
            if sky_spectrum is not None:
                data = np.column_stack([wavelength, spectrum_1d, sky_spectrum])
                header = "Wavelength(Angstrom)  Flux(ADU)  BG(ADU)"
            else:
                data = np.column_stack([wavelength, spectrum_1d])
                header = "Wavelength(Angstrom)  Flux(ADU)"
            
            # Add comments
            if comments:
                header = comments + "\n" + header
            
            # Save
            np.savetxt(output_path, data, header=header, fmt='%.6e',
                      comments='# ')
            
            logger.info(f"Saved 1D spectrum text: {output_path}")
            return output_path
        
        except Exception as e:
            logger.error(f"Failed to save 1D spectrum text: {str(e)}")
            return None

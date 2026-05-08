"""
Preprocessing module for spectroscopic data
Handles A-B difference, dark subtraction, etc.
"""

import numpy as np
import logging
from pathlib import Path
from .loader import FITSLoader

logger = logging.getLogger(__name__)


class PreprocessingPipeline:
    """Preprocessing pipeline for spectroscopic observations"""
    
    def __init__(self, spectra_dir="./spectra", pattern_a="_A.fits", 
                 pattern_b="_B.fits"):
        """
        Initialize preprocessing pipeline
        
        Parameters
        ----------
        spectra_dir : str
            Directory containing FITS files
        pattern_a : str
            Pattern for A-position frames
        pattern_b : str
            Pattern for B-position frames
        """
        self.spectra_dir = Path(spectra_dir)
        self.pattern_a = pattern_a
        self.pattern_b = pattern_b
        self.loader = FITSLoader()

    @staticmethod
    def _base_object_name(object_name):
        """Drop the automatic pair suffix from object names such as OBJ_01."""
        parts = str(object_name).rsplit('_', 1)
        if len(parts) == 2 and parts[1].isdigit():
            return parts[0]
        return str(object_name)

    @classmethod
    def _object_is_selected(cls, object_name, object_names):
        """Return True when object_name is selected for processing."""
        if not object_names:
            return True

        selected = set(object_names)
        return (
            object_name in selected
            or cls._base_object_name(object_name) in selected
        )
    
    def find_ab_pairs(self):
        """
        Find paired A-B observations based on OBJECT header keyword
        
        Returns
        -------
        pairs : list of tuples
            List of (filename_a, filename_b, object_name) tuples
        """
        if not self.spectra_dir.exists():
            logger.warning(f"Spectra directory not found: {self.spectra_dir}")
            return []
        
        # Find all FITS files
        all_files = sorted(self.spectra_dir.glob("*.fits"))
        
        # Group files by OBJECT name. A target can have multiple A/B exposures,
        # so keep every candidate instead of overwriting earlier files.
        object_groups = {}
        
        for filepath in all_files:
            # Load header to get OBJECT name
            _, header, success, msg = self.loader.load(filepath)
            if not success:
                logger.warning(f"Failed to load header for {filepath.name}: {msg}")
                continue
            
            # Extract object name and position from OBJECT keyword
            original_object = self.loader.get_header_value(header, 'OBJECT', 'UNKNOWN')
            if original_object == 'UNKNOWN':
                logger.warning(f"No OBJECT keyword in {filepath.name}")
                continue
            
            # Extract object name by removing position suffix (_A, _B)
            object_name = original_object
            position_from_object = None
            if original_object.endswith('_A'):
                object_name = original_object[:-2]  # Remove '_A'
                position_from_object = 'A'
            elif original_object.endswith('_B'):
                object_name = original_object[:-2]  # Remove '_B'
                position_from_object = 'B'
            
            # Determine position: OBJECT suffix > POSITION header > filename pattern
            position = position_from_object
            if position is None:
                position = self.loader.get_header_value(header, 'POSITION', '').upper()
            
            if position in ['A', 'B']:
                is_a = (position == 'A')
                is_b = (position == 'B')
            else:
                # Fallback to filename patterns
                is_a = (self.pattern_a in filepath.name or '_A' in filepath.name)
                is_b = (self.pattern_b in filepath.name or '_B' in filepath.name)
            
            if not (is_a or is_b):
                logger.warning(f"Cannot determine position for {filepath.name} (OBJECT: {original_object})")
                continue
            
            if object_name not in object_groups:
                object_groups[object_name] = {'A': [], 'B': []}
            
            if is_a:
                object_groups[object_name]['A'].append(filepath)
            elif is_b:
                object_groups[object_name]['B'].append(filepath)
        
        # Create pairs
        pairs = []
        for object_name, files in object_groups.items():
            a_files = sorted(files['A'])
            b_files = sorted(files['B'])
            n_pairs = min(len(a_files), len(b_files))

            if n_pairs == 0:
                logger.warning(f"Incomplete pair for {object_name}: A={a_files}, B={b_files}")
                continue

            if len(a_files) != len(b_files):
                logger.warning(
                    f"Uneven A-B count for {object_name}: "
                    f"A={len(a_files)}, B={len(b_files)}; processing {n_pairs} pair(s)"
                )

            for idx, (filepath_a, filepath_b) in enumerate(zip(a_files, b_files), start=1):
                pair_object_name = object_name
                if n_pairs > 1:
                    pair_object_name = f"{object_name}_{idx:02d}"

                pairs.append((filepath_a, filepath_b, pair_object_name))
                logger.info(
                    f"Found pair for {pair_object_name}: "
                    f"{filepath_a.name} <-> {filepath_b.name}"
                )
        
        return pairs
    
    def ab_difference(self, data_a, data_b):
        """
        Perform A-B difference to cancel dark and sky
        
        Parameters
        ----------
        data_a : 2D array
            A-position frame
        data_b : 2D array
            B-position frame
        
        Returns
        -------
        diff : 2D array or None
            A-B difference, or None if shapes don't match
        """
        if data_a is None or data_b is None:
            return None
        
        if data_a.shape != data_b.shape:
            logger.error(f"Shape mismatch: {data_a.shape} vs {data_b.shape}")
            return None
        
        # A-B difference
        diff = data_a.astype(float) - data_b.astype(float)
        
        return diff
    
    def process_pair(self, filepath_a, filepath_b, object_name=None):
        """
        Process a single A-B pair
        
        Parameters
        ----------
        filepath_a : Path
            Path to A-position FITS file
        filepath_b : Path
            Path to B-position FITS file
        object_name : str or None
            Object name. If None, extracted from header
        
        Returns
        -------
        processed_data : 2D array or None
            A-B difference image
        headers : tuple of dict or None
            (header_a, header_b)
        object_name : str
            Object name extracted from header
        success : bool
            Whether processing was successful
        error_msg : str
            Error message if unsuccessful
        """
        # Load A frame
        data_a, header_a, success_a, msg_a = self.loader.load(filepath_a)
        if not success_a:
            return None, None, None, False, f"Load A failed: {msg_a}"
        
        # Load B frame
        data_b, header_b, success_b, msg_b = self.loader.load(filepath_b)
        if not success_b:
            return None, None, None, False, f"Load B failed: {msg_b}"
        
        # A-B difference
        diff = self.ab_difference(data_a, data_b)
        if diff is None:
            return None, None, None, False, "A-B difference failed"
        
        # Extract object name from header if not provided
        if object_name is None:
            object_name = self.loader.get_header_value(header_a, 'OBJECT', 'UNKNOWN')
            if object_name == 'UNKNOWN':
                object_name = self.loader.get_header_value(header_b, 'OBJECT', 'UNKNOWN')
        
        return diff, (header_a, header_b), object_name, True, ""
    
    def process_all_pairs(self, object_names=None):
        """
        Process all A-B pairs in directory

        Parameters
        ----------
        object_names : list or None
            Object names to process. If None, process every A-B pair.
        
        Yields
        ------
        result : dict
            Dictionary with keys: 'object_name', 'data', 'headers', 
            'filepath_a', 'filepath_b', 'success', 'error'
        """
        pairs = self.find_ab_pairs()
        
        if not pairs:
            logger.warning("No A-B pairs found")
            return
        
        for filepath_a, filepath_b, object_name in pairs:
            if not self._object_is_selected(object_name, object_names):
                continue

            try:
                diff, headers, object_name_final, success, error_msg = \
                    self.process_pair(filepath_a, filepath_b, object_name)
                
                yield {
                    'object_name': object_name_final,
                    'data': diff,
                    'headers': headers,
                    'filepath_a': filepath_a,
                    'filepath_b': filepath_b,
                    'success': success,
                    'error': error_msg
                }
            except Exception as e:
                # Catch any unexpected errors and continue
                logger.error(f"Error processing pair {filepath_a.name} <-> {filepath_b.name}: {str(e)}")
                yield {
                    'object_name': object_name,
                    'data': None,
                    'headers': None,
                    'filepath_a': filepath_a,
                    'filepath_b': filepath_b,
                    'success': False,
                    'error': str(e)
                }

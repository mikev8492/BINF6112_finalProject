## Development Log:

Date 04/17/2026
- Updated `csv_output.py` module:
    - Added `CreateCSV` class to create the CSV output file for the user.
        - `create_csv_output` function writes the enzyme name, motif, cut site, observed count, and the start position locations as a row in the CSV file. 
- Refactored `main.py`:
    - Updated function calls in `main()` to create the CSV output file. 

Date: 04/09/2026
- Updated environment to include `matplotlib`
- Added `annotate_plasmid` function to  `output.py` module:
    - Uses results dictionary and plasmid sequence string to generate a circular map.
    - A circle is drawn with radius of 1.0, and the enzyme cut locations are converted from a bp to an angle on the circle, starting from the top:
        angle = 90° − (bp position / sequence length) × 360°
    - A tick mark is then drawn at the calculated angle. 
    - png graph is generated and saved to `results/plasmid_map.png`

Date: 04/01/2026
- Added `Motifs` class to `motif_locator.py` module:
    - `Motifs` class locates all instances of an enzyme motif in the plasmid sequence. 
    - Class functions added:
        - `array_set`: transforms the plasmid sequence to a numpy array.
        - `motif_search`: searches for all locations of an enzyme motif in the plasmid sequence
        - `get_motif_results`: Reports the results from `motif_search` in a dictionary containing enyzmes as keys and a list containing the corresponding motif sequence, cut site, observed motif counts, and the starting indeces of all observed motif counts as the values.
- Refactored `main.py`:
    - Updated function calls in `main()` to report the resulting dictionary from `Motifs.get_motif_results()`


Date: 03/25/2026
- Refactored `main.py`:
    - Create database function updated to conditional
    - User arguments updated to include: 
        - --interface: optional terminal interface
        - --enzymes: user list of enzymes to map with
        - --output: output filename
    - Added argument validation
    - Refactored main function calls to pass user arguments to input classes
- Refactored `input.py` module:
    - Removed error handling from `Sequence` class (handled with argparse now)
    - Replaced `Enzymes.load_RE()` function, database is read with csv.dictReader object at the filtering function `Enzymes.filter_enzymes()`instead 
    - Updated user interface to be optional
    - Updated default enzyme list to use constant defined in `main.py`


Date: 03/19/2026
- updated relative filepaths to run main.py from root project directory.
- Added `environment-alternative.yml` to resolve cross compatibility conflicts. 
- Refactored `main.py`: 
    1. Moved CLI header to `input.py` module to allow User interface as optional argument.
    2. Added `ArgParse` functionality to simplify CLI.
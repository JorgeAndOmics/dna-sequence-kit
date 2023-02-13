
# DNA Sequence Kit

This Python script contains various functions for analyzing DNA sequences. It requires the Biopython library to run. The functions in this script can be used for tasks such as calculating GC content, codon usage, performing multiple sequence alignment, removing introns, adding restriction enzyme sites, converting sequence formats, exporting sequences, batch processing sequences, and searching sequence databases.

## Dependencies

This script requires the following dependencies:

-   Biopython

## Functions

The following functions are provided:

-   `calculate_gc_content(sequence)`: Calculates the GC content of a DNA sequence.
-   `calculate_codon_usage(sequence)`: Calculates the codon usage of a DNA sequence.
-   `perform_restriction_analysis(sequence, enzyme_name)`: Performs a restriction enzyme analysis on a DNA sequence.
-   `perform_multiple_alignment(sequences, algorithm)`: Performs a multiple sequence alignment of DNA sequences.
-   `remove_introns(sequence)`: Removes introns from a DNA sequence.
-   `add_restriction_enzyme_sites(sequence, enzyme_name)`: Adds a specific restriction enzyme recognition site to a DNA sequence.
-   `convert_sequence_format(sequence, format)`: Converts a DNA sequence to a different file format.
-   `export_sequences(sequences, format)`: Exports DNA sequences in a format that can be easily imported into other tools.
-   `batch_process_sequences(sequences)`: Performs a batch process on multiple DNA sequences.
-   `search_sequence_databases(sequence, database_name)`: Performs a search against a DNA sequence database.

## `Use Example`
### `perform_restriction_analysis(sequence, enzyme_name)`

Performs a restriction enzyme analysis on a DNA sequence.

**Arguments:**

-   `sequence` (str): The DNA sequence to analyze.
-   `enzyme_name` (str): The name of the restriction enzyme to use.

**Returns:**

-   `fragments` (list): A list of the fragments produced by the restriction enzyme digestion of the sequence.

**Usage:**

    from DNA_sequence_analysis_script import perform_restriction_analysis
    
    sequence = "ATCGATCGATCGATCG"
    enzyme_name = "EcoRI"
    fragments = perform_restriction_analysis(sequence, enzyme_name)
    print(fragments)

## Contribution

If you would like to contribute to this script, please fork the repository, make your changes, and submit a pull request. Any contributions, bug fixes, and suggestions are welcome!

## License

This script is licensed under the MIT license.

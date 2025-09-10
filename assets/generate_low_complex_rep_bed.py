#!/usr/bin/env python3
"""
RepeatMasker Output Parser to low complex repetitive elements BED file

This script parses RepeatMasker files to extract genomic repetitive elements
and low complexity regions into a bed file for filtering.

Usage:
    python generate_low_complex_rep_bed.py INPUT_FILE OUTPUT_FILE [OPTIONS]

Examples:
    python generate_low_complex_rep_bed.py input.out.gz output.bed
    python generate_low_complex_rep_bed.py input.out.gz output.bed --valid-chromosomes chr1,chr2,chrX
    python generate_low_complex_rep_bed.py input.out.gz output.bed --skip-quality-checks
"""

# https://genome-blog.gi.ucsc.edu/blog/2016/12/12/the-ucsc-genome-browser-coordinate-counting-systems/

import click
import pandas as pd
import gzip
import logging
import sys
from pathlib import Path
from typing import Optional, TextIO, Union


# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)



def open_file(file_path: Union[str, Path]) -> TextIO:
    """
    Open a file, handling both regular and gzipped files.
    
    Args:
        file_path: Path to the file
        
    Returns:
        File handle
        
    Raises:
        FileNotFoundError: If file doesn't exist
        IOError: If file cannot be opened
    """
    file_path = Path(file_path)
    
    if not file_path.exists():
        raise FileNotFoundError(f"Input file not found: {file_path}")
    
    try:
        if file_path.suffix == '.gz':
            return gzip.open(file_path, 'rt')
        else:
            return open(file_path, 'r')
    except IOError as e:
        raise IOError(f"Cannot open file {file_path}: {e}")


def parse_repeatmasker_file(file_path: Union[str, Path], 
                          valid_chromosomes: Optional[list] = None) -> pd.DataFrame:
    """
    Parse RepeatMasker file.
    
    Args:
        file_path: Path to RepeatMasker file
        valid_chromosomes: List of valid chromosome names to filter
        
    Returns:
        DataFrame with parsed repetitive elements
        
    Raises:
        ValueError: If file format is invalid
        IOError: If file cannot be read
    """
    logger.info(f"Parsing RepeatMasker file: {file_path}")
    
    # Only read the 5 columns we need (0-indexed)
    # Col 4=chr, 5=start, 6=end, 9=rep_name, 10=rep_class
    usecols = [4, 5, 6, 9, 10]
    col_names = ['chr', 'start', 'end', 'repetition_name', 'repetition_class']
    
    # Memory-efficient data types
    dtype_dict = {
        'chr': 'category',           
        'start': 'int32',           
        'end': 'int32',
        'repetition_name': 'category',
        'repetition_class': 'category'
    }
    
    try:
        with open_file(file_path) as file:
            df = pd.read_csv(
                file,
                sep=r'\s+',              # Split on any whitespace
                skiprows=3,              # Skip header rows
                usecols=usecols,         # Only read columns we need -> less time
                names=col_names,
                engine='python',         # Required for regex separator
                dtype=dtype_dict,       
            )
            
        logger.info("Successfully parsed RepeatMasker file.")
                
    except Exception as e:
        raise ValueError(f"Failed to parse RepeatMasker file: {e}")
    
    # Remove empty rows (last row is often empty)
    df = df.dropna(subset=['chr', 'start', 'end'])
    
    # Filter valid chromosomes if specified
    if valid_chromosomes:
        initial_count = len(df)
        df = df[df['chr'].isin(valid_chromosomes)]
        filtered_count = initial_count - len(df)
        if filtered_count > 0:
            logger.info(f"Filtered out {filtered_count} rows with invalid chromosomes")
    
    # Substract 1 from start to convert to 0-based indexing
    df['start'] = df['start'] - 1

    logger.info(f"Parsed {len(df)} repetitive elements")
    return df


def validate_data_quality(df: pd.DataFrame) -> pd.DataFrame:
    """
    Perform data quality checks and cleaning.
    
    Args:
        df: DataFrame to validate
        
    Returns:
        Cleaned DataFrame
        
    Raises:
        ValueError: If critical data quality issues are found
    """
    logger.info("Performing data quality checks...")
    
    initial_count = len(df)
    
    # Check for missing critical values
    missing_coords = df[['chr', 'start', 'end']].isnull().any(axis=1)
    if missing_coords.any():
        missing_count = missing_coords.sum()
        logger.warning(f"Removing {missing_count} rows with missing coordinates")
        df = df[~missing_coords].copy()
    
    # Validate coordinates  
    invalid_coords = (df['start'] >= df['end']) | (df['start'] < 0) | (df['end'] < 0)
    if invalid_coords.any():
        invalid_count = invalid_coords.sum()
        logger.warning(f"Removing {invalid_count} rows with invalid coordinates")
        df = df[~invalid_coords].copy()
    
    # Check for duplicates
    duplicates = df.duplicated(subset=['chr', 'start', 'end', 'repetition_name'])
    if duplicates.any():
        duplicate_count = duplicates.sum()
        logger.warning(f"Removing {duplicate_count} duplicate entries")
        df = df[~duplicates].copy()
    
    # Validate chromosome names (convert to string for regex check)
    chr_str = df['chr'].astype(str)
    weird_chr = chr_str.str.contains(r'[^a-zA-Z0-9_]', na=False)
    if weird_chr.any():
        weird_count = weird_chr.sum()
        logger.warning(f"Found {weird_count} entries with unusual chromosome names")
    
    final_count = len(df)
    removed_count = initial_count - final_count
    
    if removed_count > 0:
        logger.info(f"Data quality checks removed {removed_count} rows "
                   f"({removed_count/initial_count*100:.1f}%)")
    
    if final_count == 0:
        raise ValueError("No valid data remaining after quality checks")
    
    return df


def save_results(df: pd.DataFrame, output_path: Union[str, Path]) -> None:
    """
    Save results to BED-like format file.
    
    Args:
        df: DataFrame to save
        output_path: Output file path
        
    Raises:
        IOError: If file cannot be written
    """
    output_path = Path(output_path)
    
    try:
        # Create output directory if it doesn't exist
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Save to BED-like format
        df[['chr', 'start', 'end']].to_csv(output_path, index=False, sep="\t", header=False)
        logger.info(f"Results saved to: {output_path}")
        
        # Log summary statistics
        logger.info(f"Output summary:")
        logger.info(f"  Total regions: {len(df):,}")
        logger.info(f"  Chromosomes: {df['chr'].nunique()}")
        
    except Exception as e:
        raise IOError(f"Failed to save results to {output_path}: {e}")


@click.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.argument('output_file', type=click.Path())
@click.option(
    '--valid-chromosomes',
    type=str,
    default='chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY',
    help="Comma-separated list of valid chromosome names to include. Default: chr1,chr2,...,chr22,chrX,chrY"
)
@click.option(
    '--skip-quality-checks',
    is_flag=True,
    help="Skip data quality validation checks"
)

def main(input_file, output_file, valid_chromosomes, skip_quality_checks):
    """Parse RepeatMasker .out files into structured format.
    
    INPUT_FILE: Input RepeatMasker .out file (can be gzipped)
    OUTPUT_FILE: Output BED-like format file path
    """
    logging.getLogger().setLevel(logging.DEBUG)
    
    try:
        # Parse comma-separated chromosome string
        if valid_chromosomes:
            valid_chr_list = [chrom.strip() for chrom in valid_chromosomes.split(',') if chrom.strip()]
        else:
            valid_chr_list = None
            
        logger.info(f"Using chromosomes: {valid_chr_list}")
        
        # Parse RepeatMasker file
        df = parse_repeatmasker_file(input_file, valid_chr_list)
        
        # Perform data quality checks unless skipped
        if not skip_quality_checks:
            df = validate_data_quality(df)
        
        # Save results
        save_results(df, output_file)
        
        logger.info("Processing completed successfully!")
        
    except Exception as e:
        logger.error(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()

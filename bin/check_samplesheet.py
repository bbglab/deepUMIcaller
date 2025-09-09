#!/usr/bin/env python


"""Provide a command line tool to validate and transform tabular samplesheets."""


import argparse
import csv
import logging
import sys
from collections import Counter
from pathlib import Path


logger = logging.getLogger()


requirementsDict = { "mapping": ["fastq_1" , "fastq_2", "read_structure"],
                    "groupbyumi": ["bam"],
                    "filterconsensus": ["bam"],
                    "calling": ["duplexbam", "csi"],
}

class RowChecker:
    """
    Define a service that can validate and transform each given row.

    Attributes:
        modified (list): A list of dicts, where each dict corresponds to a previously
            validated and transformed row. The order of rows is maintained.

    """

    VALID_FORMATS = (
        ".fq.gz",
        ".fastq.gz",
    )

    def __init__(
        self,
        sample_col="sample",
        init_cols = list(),
        step='mapping',
        **kwargs,
    ):
        """
        Initialize the row checker with the expected column names.

        Args:
            sample_col (str): The name of the column that contains the sample name
                (default "sample").
            first_col (str): The name of the column that contains the first (or only)
                FASTQ file path (default "fastq_1").
            second_col (str): The name of the column that contains the second (if any)
                FASTQ file path (default "fastq_2").
            read_structure (str): The name of the column that contains the read
                structure for the given sample.

        """
        super().__init__(**kwargs)
        self._sample_col = sample_col
        self._step = step       
        self._columns = init_cols
        self._read_structure_col = "read_structure"
        self._first_col = "fastq_1"
        self._second_col = "fastq_2"
        self._seen = set()
        self.modified = []

    def validate_and_transform(self, row):
        """
        Perform all validations on the given row and insert the read pairing status.

        Args:
            row (dict): A mapping from column headers (keys) to elements of that row
                (values).

        """
        self._validate_sample(row)
        if self._step == "mapping":
            self._validate_read_structure(row)
            self._validate_pair(row)

        self._seen.add(tuple(row[column] for column in self._columns))
        self.modified.append(row)

    def _validate_sample(self, row):
        """Assert that the sample name exists and convert spaces to underscores."""
        assert len(row[self._sample_col]) > 0, "Sample input is required."
        # Sanitize samples slightly.
        row[self._sample_col] = row[self._sample_col].replace(" ", "_")

    def _validate_pair(self, row):
        """Assert that read pairs have the same file extension. Report pair status."""
        assert (
            Path(row[self._first_col]).suffixes[-2:] == Path(row[self._second_col]).suffixes[-2:]
        ), "FASTQ pairs must have the same file extensions."

    def _validate_read_structure(self, row):
        """Assert that the read structure is a valid duplex read structure."""
        assert len(row[self._read_structure_col].split(' ')) == 2, (
            "Two read structures must be provided."
        )

    def validate_unique_samples(self):
        """
        Assert that the combination of sample name and FASTQ filename is unique.

        In addition to the validation, also rename the sample if more than one sample,
        FASTQ file combination exists.

        """
        assert len(self._seen) == len(self.modified), "The pair of sample name and FASTQ must be unique."
        counts = Counter(pair[0] for pair in self._seen)
        seen = Counter()
        for row in self.modified:
            sample = row[self._sample_col]
            seen[sample] += 1
            if counts[sample] > 1:
                row['id'] = f"{sample}_LPART{seen[sample]}"
            else:
                row['id'] = sample


def read_head(handle, num_lines=10):
    """Read the specified number of lines from the current position in the file."""
    lines = []
    for idx, line in enumerate(handle):
        if idx == num_lines:
            break
        lines.append(line)
    return "".join(lines)


def sniff_format(handle):
    """
    Detect the tabular format.

    Args:
        handle (text file): A handle to a `text file`_ object. The read position is
        expected to be at the beginning (index 0).

    Returns:
        csv.Dialect: The detected tabular format.

    .. _text file:
        https://docs.python.org/3/glossary.html#term-text-file

    """
    peek = read_head(handle)
    handle.seek(0)
    sniffer = csv.Sniffer()
    if not sniffer.has_header(peek):
        logger.critical(f"The given sample sheet does not appear to contain a header.")
        sys.exit(1)
    dialect = sniffer.sniff(peek)
    return dialect


def check_samplesheet(file_in, file_out, step = "mapping"):
    """
    Check that the tabular samplesheet has the structure expected by nf-core pipelines.

    Validate the general shape of the table, expected columns, and each row.

    Args:
        file_in (pathlib.Path): The given tabular samplesheet. The format can be either
            CSV, TSV, or any other format automatically recognized by ``csv.Sniffer``.
        file_out (pathlib.Path): Where the validated and transformed samplesheet should
            be created; always in CSV format.

    Example:
        This function checks that the samplesheet follows the following structure ::

            sample,fastq_1,fastq_2,read_structure
            sample,fastq_1,fastq_2,read_structure
            SAMPLE_DUPLEX_SEQ,SAMPLE_DUPLEX_SEQ.R1.fq.gz,SAMPLE_DUPLEX_SEQ.R2.fq.gz,10M1S+T 10M1S+T
            SAMPLE_SINGLE_UMI,SAMPLE_SINGLE_UMI.R1.fq.gz,SAMPLE_SINGLE_UMI.R2.fq.gz,12M+T +T

    """
    required_columns = {"sample"}
    required_columns.update(requirementsDict[step])

    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_in.open(newline="") as in_handle:
        reader = csv.DictReader(in_handle, dialect=sniff_format(in_handle))
        # Validate each row.
        checker = RowChecker(step = step, init_cols= required_columns)
        for i, row in enumerate(reader):
            try:
                checker.validate_and_transform(row)
            except AssertionError as error:
                logger.critical(f"{str(error)} On line {i + 2}.")
                sys.exit(1)
        checker.validate_unique_samples()
    # Remove 'original_sample' from the output header
    header = list(reader.fieldnames) + ['id']
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_out.open(mode="w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, header, delimiter=",")
        writer.writeheader()
        for row in checker.modified:
            writer.writerow(row)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv [--log-level INFO --step mapping]",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Tabular input samplesheet in CSV or TSV format.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Transformed output samplesheet in CSV format.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    parser.add_argument(
        "-s",
        "--step",
        help="The desired step (default WARNING).",
        choices=("mapping", "filterconsensus", "calling"),
        default="mapping",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    check_samplesheet(args.file_in, args.file_out, args.step)


if __name__ == "__main__":
    sys.exit(main())

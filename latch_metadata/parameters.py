
from dataclasses import dataclass
import typing
import typing_extensions

from flytekit.core.annotation import FlyteAnnotation

from latch.types.metadata import NextflowParameter
from latch.types.file import LatchFile
from latch.types.directory import LatchDir, LatchOutputDir

# Import these into your `__init__.py` file:
#
# from .parameters import generated_parameters

generated_parameters = {
    'input_genomes': NextflowParameter(
        type=LatchDir,
        default=None,
        section_title='Input/output options',
        description='Path to input directory of input genomes in FASTA format ending in .fa',
    ),
    'outdir': NextflowParameter(
        type=typing_extensions.Annotated[LatchDir, FlyteAnnotation({'output': True})],
        default=None,
        section_title=None,
        description='The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.',
    ),
    'genome_metadata': NextflowParameter(
        type=LatchFile,
        default=None,
        section_title=None,
        description='Path to genome metadata TSV file',
    ),
    'peptides_fasta': NextflowParameter(
        type=LatchFile,
        default=None,
        section_title='Databases',
        description='Path to FASTA file of peptide sequences database to compare input sequences to.',
    ),
    'models_dir': NextflowParameter(
        type=LatchDir,
        default=None,
        section_title=None,
        description='Path to directory where bioactivity machine learning models are.',
    ),
    'models_list': NextflowParameter(
        type=LatchFile,
        default=None,
        section_title=None,
        description='TXT file of list of models within the models_dir directory to make predictions for.',
    ),
}


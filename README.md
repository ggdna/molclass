# Designing disease diagnostic molecular classifiers

### Training diagnostic model

The classifier.py file contains functions that can train and validate models for gene expression based disease diagnostic tasks that could be feasibly implemented by molecular classifier frameworks.

### Sequence selection

The secstruc.py file contains wrappers for helpful NUPACK tools, taken from [Multistrand](https://github.com/DNA-and-Natural-Algorithms-Group/multistrand).

The sq.py file uses these functions to identify desirable sequences for probes to be used in molecular classifiers.

### Concentration determination for encoding weights in probe concentration

The conc_det.py file can be used to determine experimental procedures for an approach that encodes classifier weights in probe concentration rather than unique probe count.

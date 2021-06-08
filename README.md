# Genser - the GENome Size EstimatoR

This is a tool that will try to estimate a genome size based on
mapping reads to a draft assembly.  It collects coverage statistics
from a table or directly from a BAM file (if invoked with `-b` and
`samtools` is installed), and attempts to fit a series of negative
binomial distributions using expectation maximization.

The genome size estimate is split into diploid sequences (with
a normal coverage distribution) haploid sequences (having half the
expected coverage of diploid sequences), and coverages corresponding
to collapsed repeats with twice, thrice, or four times the diploid
coverage.  Use the `-p` option to display a plot of the coverage,
distributions, and residuals, and `-v` to show progress while the
distributions converge.

Haploid and diploid size are calculated by integrating the coverage
statistics and dividing by the expected coverages (to account for
collapsed repeats), while low coverage and zero coverage loci are
reported as is.

# FastCPA

FastCPA is a correlation power analysis attack different from the traditional CPA. While its accuracy is equivalent to the traditionnal CPA, it is orders of magnitude faster when the number of traces is large. This work has been published under the following reference:
Meunier, Quentin L. "FastCPA: Efficient correlation power analysis computation with a large number of traces." Proceedings of the Sixth Workshop on Cryptography and Security in Computing Systems. 2019.
[Link](https://dl.acm.org/doi/abs/10.1145/3304080.3304082)


## Usage

FastCPA is implemented in the fastcpa.c file, when the traditionnal CPA is provided in the cpa\_ref.c file. The implemented attacks target an AES. More specifically, the targeted internal value is the value of the state after the first SBox, although this could easily be adapted to any CPA. In order to compile the files, a Makefile is provided.

The characteristics of the trace files to use must be specified in a config.h file, with the following defines:
* START\_SAMPLE: index of the first sample to use for the analysis
* END\_SAMPLE: index of the first sample not to use
* FILE\_SAMPLES: number of samples per single trace in the trace file
* NB\_TRACES: number of traces in the trace file

For example, for a trace file with 2,000 samples per trace, setting START\_SAMPLE to 0 and END\_SAMPLE to 2,000 will make the attack consider all the samples in each trace.

An example trace file containing 10,000 traces of 1,600 samples is provided. The trace file is expected to contain all the samples values written consecutively in binary (trace after trace), either in float or double type (the typdef of float\_f must be defined accordingly in the common.h file).

The file key.raw must contain the 16-byte key in binay, while the file textin.raw must contain the 16-byte plaintexts for all traces written consecutively in binary.

After compiling the files, running the binary without argument will perform the attack.



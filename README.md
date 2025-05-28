# Genome Function Inference Pipeline Using KEGG KO Annotations

This pipeline analyzes KO annotations of draft microbial genomes (e.g. SAGs, MAGs) to infer their potential metabolic functions. 

## Installation
This pipeline requires Python 3.10 or higher
No additional Python packages are strictly required, as the pipeline mainly uses standard libraries and subprocess calls.

Clone the repository:

```bash
git clone https://github.com/kyuinHwang/ko2pathway.git
cd GenomeHabitatClassifier
```
## Usage

Run the entire pipeline:
``` bash
run_pipeline.py [-h] --ko_dir KO_DIR --module_file MODULE_FILE --keyenzyme_dir ./KeyEnzymes --output_dir OUTPUT_DIR
```

Or run each script individually for modular control:

## Scripts
1_track_steps.py: Checks KO presence for each step in KEGG modules.

2_evaluate_modules.py: Determines whether each module is likely functional, based on majority step coverage and key step presence.

3_check_reactions.py: Detects simple reactions based on presence of relevant KO(s).

## Input
KO_DIR : Directory including .ko files of queried genomes. .ko is a tab-delimited file listing each gene's locus tag and its corresponding KO (KEGG Orthology) assignment, one per line (e.g., EXCAFAC_00031 K00302). This file can be generated through BLAST searches as described in our publication (see Reference section below), or via the KEGG annotation web server:

MODULE_FILE : A file containing the definitions of KEGG modules. This can be downloaded from the KEGG FTP server (module/module.gz) by users with a valid KEGG subscription.


## Files
./KeyEnzymes/:
This directory contains files listing the KOs responsible for simple metabolic reactions. In contrast, the complex pathways of interest (based on KEGG modules) are defined within 1_track_steps.py, and the threshold conditions (e.g., minimum step coverage, key steps) are specified in 2_evaluate_modules.py.

> **Note on KEGG Version:**
> The selection of enzymes (KOs) and threshold were curated based on KEGG release from October 23, 2019. Since the module definition and participating KOs may have changed, users working with a different KEGG version should manually verify
>
> the list of target module IDs in 1_track_steps.py,
>
> threshold and key steps defined in 2_evalulate_modules.py, and 
>
> the KO lists in files under ./KeyEnzymes

./examples/:
Contains KO profiles of three public genomes belonging to the genus g__SURF-13.

./extended_examples/:
This directory provides example input files used in our full research pipeline, which integrates additional sources of evidence (e.g., FeGenie results, BLAST hits) that are not handled by the core scripts in this repository. These files were used in the published study but are not part of the reproducible module.Notably, the presence of form I coxL was manually determined based on a phylogenetic tree and is not  part of either the core scripts or the extended example files.


## Output
By default (run_pipeline.py), the results are summarized into two tab-delimited tables:

./OUTPUT_DIR/Pathway.txt: Summarizes the presence (Y) or absence (N) of complex pathways (based on KEGG modules) for each genome.

./OUTPUT_DIR/Reactions.txt: Summarizes the presence (Y) or absence (N) of simple reactions for each genome.
Each row corresponds to a genome, and each column to a specific pathway or reaction.

## License

This repository is released under the MIT License.

## Contact
If you encounter a problem or have a question, please open an issue on this repository:
👉 [Submit an issue](https://github.com/kyuinHwang/ko2pathway/issues)

For direct inquiries, you may contact the maintainer at: rbdls77@gmail.com

## Reference

This pipeline was developed as part of a research project on single-cell genomics analysis of Mercer Subglacial Lake in Antarctica.

The preprint is available on Research Square:
https://www.researchsquare.com/article/rs-4392950/v1

(This link will be updated upon journal publication.)

The key enzymes and thresholds are available in the Supplementary Tables.

Parts of this README were written or revised with the help of AI to enhance clarity and precision.

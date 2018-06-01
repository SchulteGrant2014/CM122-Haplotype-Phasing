# CM124 - Haplotype Phasing
### *This project is for UCLA CS CM124 (Computational Genetics).*

This is a bioinformatics project for phasing the correct haplotypes from genotype information derived from DNA sequencing. Specifically, we strive to accurately generate the pair of chromosomes that forms an individual's genotype.

The main challenge of this project is that the number of possible chromosome pairs that make up an individual's genotype grows exponentially with each heterozygous position (value 1) of that genotype. This means that the number of possible haplotype phases doubles with each heterozygous position. And because we are dealing with 50 individuals, the number of possible genotypes becomes much too large for any computer to handle.

As a result, with nearly 40,000 SNPs and 50 individuals, it is impossible to phase all genotypes at once with current technology. Therefore, a modified method of phasing is required, as opposed to more general techniques. In this project, we attempt to implement a "windowed" approach, in which we phase sections of the genotype individually and then stitch them together to form the final haplotypes.

Collaborators:
* Andy Hsu
* Ames Ma
* Grant Schulte


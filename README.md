# Ahemp
Acropora hemprichii genome structural and functional annotation
- The quality of the assembled genome was assessed using BlobToolKit (BTK), full handout for how to run it can be found here (https://github.com/blobtoolkit/tutorials/tree/main/futurelearn).Also, a course on BTK is avaliable as well https://www.futurelearn.com/courses/eukaryotic-genome-assembly-how-to-use-blobtoolkit-for-quality-assessment

- After first BTK check, ~90 contigs were removed because of low coverage or no-related taxa.
- BTK re-runed with the filtered assembly and another 528 contigs were removed for low coverage or no-related taxa contigs.
## identifying and maksing Repeats
- Repeats in Ahemp and avliable Acropora's coral genomes using RepeatModeler

````bash
singularity pull dfam-tetools-latest.sif docker://dfam/tetools:latest
singularity run dfam-tetools-latest.sif
singularity run ../../dfam-tetools-latest.sif BuildDatabase -name Ahemp_genome ../Ahemp_genome/Ahemp.gapclosed.fasta
singularity run ../../dfam-tetools-latest.sif RepeatModeler -database Ahemp_genome -LTRStruct -threads 40
````

# Ahemp
_Acropora hemprichii_ genome structural and functional annotation
- The quality of the assembled genome was assessed using BlobToolKit (BTK), full handout for how to run it can be found here (https://github.com/blobtoolkit/tutorials/tree/main/futurelearn). Also, a course on BTK is available as well https://www.futurelearn.com/courses/eukaryotic-genome-assembly-how-to-use-blobtoolkit-for-quality-assessment
- For the second round, i used the new nextflow pipeline for BTK (https://pipelines.tol.sanger.ac.uk/blobtoolkit/0.2.0).
     - it requires:
          - download BUSCO database.
          - converts sequencing data to .CRAM
          ````bash
          java -jar build/libs/picard.jar FastqToSam F1=/home/fiesingera/proj/Ahemp_genome/Ahemp_raw_data/Ahem_1.fastq.gz F2=/home/fiesingera/proj/Ahemp_genome/Ahemp_raw_data/Ahem_2.fastq.gz  O=Ahemp_unaligned.bam SM=Ahemp
          samtools view -@ 30 -T Ahemp_final.fasta -C -o Ahemp.cram Ahemp_unaligned.bam
          ````
          - prepare a full metadata.yaml file, [example](https://github.com/blobtoolkit/blobtoolkit/tree/main/src/blobtoolkit-pipeline/src#configuration)
          - prepare .csv input sheet
          ````
          sample,datatype,datafile
          Ahemp,ont,Ahmep.cram
          ````
          - run the nextflow command
          ````bash
          nextflow run sanger-tol/blobtoolkit --input Ahemp.csv --fasta Ahemp_final.fasta –-accession Ahemp --taxon 213635 --yaml Ahemp.yaml --taxdump /share/databases/taxdump --blastp /share/databases/taxdump --blastn /share/databases/uniprot --blastx /share/databases/uniprot
          ````
- After the first BTK check, ~90 contigs were removed because of low coverage or non-related taxa (Ahemp.gapclosed_f1.fasta).
- BTK was re-runed with the filtered assembly and another 528 contigs were removed for low coverage or non-related taxa contigs (Ahemp.gapclosed_f2.fasta).
- After submitting the genome to NCBI, you have been asked to remove several contigs, which we did using "funannotate clean" and manually. we end up with 33983 contigs/scaffolds
# clean assembly 
````bash
funannotate clean -i Ahemp.gapclosed_f2.fasta -m 200 -o Ahemp_clean.fasta
````
# Identify rRNA using [barrnap](https://github.com/tseemann/barrnap)

````bash
./barrnap -q -k euk Ahemp_final.fasta --threads 50 --outseq Ahemp_rrna.fasta > Ahemp_rrna.gff 
````

## Identifying and masking Repeats
- We Identify repeats in Ahemp using RepeatModeler and EDTA, to confirm the repeats percentages.
### [EDTA](https://github.com/oushujun/EDTA) 

````bash
docker pull oushujun/edta:2.0.0
docker run -v $PWD:/in -w /in oushujun/edta:2.0.0 EDTA.pl --genome Ahemp_final.fasta --threads 30
````


### RepeatModeler



````bash
docker run -v $PWD:/in -w /in dfam/tetools:latest BuildDatabase -name Ahemp_genome Ahemp_final.fasta
docker run -v $PWD:/in -w /in dfam/tetools:latest RepeatModeler -database Ahemp_genome -LTRStruct -threads 30
````

- Then merged the indetified repeats from RepeatModeler and EDTA

- available Repeats in  Acropora's coral genomes
  
````bash

for i in `ls *.fna|sed 's/_genomic.fna//g`;
do
    docker run -v $PWD:/in -w /in dfam/tetools:latest BuildDatabase -name $i ${i}_genomic.fna
    docker run -v $PWD:/in -w /in dfam/tetools:latest RepeatModeler -database $i -LTRStruct -threads 40;
done
````


- Repeats database assembly

````bash
cat *-families.fa > Acropora_RE_DB.fsa
unsearch -fastx_uniques Acropora_RE_DB.fsa -fastaout Acropora_RE_DB.faa
````


- Repeats masking in Ahemp using Acropora repeats DB

````bash
docker run -v $PWD:/in -w /in dfam/tetools:latest RepeatMasker Ahemp.gapclosed_f2.fasta -lib Acropora_RE_DB.faa -pa 8 -norna -xsmall 
````

- How is the distribution of repeats by types

````bash
grep '>' Ahemp_genome_f2-families.fa | sed -r 's/.+#//' | sed -r 's/\s+.+//' | sort | uniq -c
````
- the full report of the repeats percentage will outputted after masking in ```Ahemp.gapclosed_f2.fasta.tbl``` file
## Preparing RNASeq evidence 

- Indexing the assembled genome
````bash

STAR --runThreadN 50 --runMode genomeGenerate --genomeDir Ahemp_index --genomeFastaFiles Ahemp.gapclosed_f2.fasta --genomeSAindexNbases 10
````

- Mapping the RNA-Seq reads to the assembled genome
  
````bash
STAR --runThreadN 30 --genomeDir Ahemp_index --readFilesIn M_19_2595_HE1-33-T1_D703-D505_L008_R1_001.fastq.gz M_19_2595_HE1-33-T1_D703-D505_L008_R2_001.fastq.gz --readFilesCommand "gunzip -c" --outSAMtype  BAM SortedByCoordinate --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix M_19_2595_HE1-33-T1_D703-D505_ --limitBAMsortRAM 10000000000

STAR --runThreadN 30 --genomeDir Ahemp_index --readFilesIn M_19_2596_HE1-36-T1_D703-D506_L008_R1_001.fastq.gz M_19_2596_HE1-36-T1_D703-D506_L008_R2_001.fastq.gz --readFilesCommand "gunzip -c" --outSAMtype  BAM SortedByCoordinate --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix M_19_2596_HE1-36-T1_D703-D506_ --limitBAMsortRAM 10000000000

STAR --runThreadN 30 --genomeDir Ahemp_index --readFilesIn M_19_2597_Ahem_D704-D505_L008_R1_001.fastq.gz M_19_2597_Ahem_D704-D505_L008_R2_001.fastq.gz --readFilesCommand "gunzip -c" --outSAMtype  BAM SortedByCoordinate --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix M_19_2597_Ahem_D704-D505_ --limitBAMsortRAM 10000000000
````

- Merge all mapping files "bam" into one

````bash
samtools merge Ahemp_RNASeqAll.STAR.bam M_19_2595_HE1-33-T1_D703-D505_sortedByCoord.out.bam M_19_2596_HE1-36-T1_D703-D506_sortedByCoord.out.bam M_19_2597_Ahem_D704-D505_sortedByCoord.out.bam
````
## Preparing expression evidence based on StringTie  

````bash
stringtie -p 30 -o Ahemp_RNASeqAll.Stringtie.gtf Ahemp_RNASeqAll.STAR.bam

#to get the gtf file details "optional" 
grep -v "#" Ahemp_RNASeqAll.Stringtie.gtf  | cut -f3 | sort | uniq -c
````

## Preparing transcripts evidence

````bash
gtf_genome_to_cdna_fasta.pl Ahemp_RNASeqAll.Stringtie.gtf Ahemp.gapclosed_f2.fasta > Ahemp_RNASeqAll.transcripts.fasta
````

## Preparing protein evidence

 - Search for curated proteins sequences for "Acropora" on UniProt database and save them as "uniprot_Acropora.faa"


## funannotate

- download/pull the image from docker hub

````bash
docker pull nextgenusfs/funannotate
````

- download bash wrapper script (optional)

````bash
wget -O funannotate-docker https://raw.githubusercontent.com/nextgenusfs/funannotate/master/funannotate-docker
````

- you might need to make this executable on your system
````bash
chmod +x /path/to/funannotate-docker
````
### Install GeneMark
- Download GeneMark-ES/ET/EP+ ver 4.71_lic and key from http://topaz.gatech.edu/GeneMark/license_download.cgi  
- extract folders and copy the key as "cp -r gm_key/ ~/.gm_key/"
- locate the perl "which perl", change all script perl tag from gmes folder "perl change_path_in_perl_scripts.pl 'perl PATH'"
- set the path for GeneMark using "export GENEMARK_PATH=~/gmes/"
- making sure that "gmes_petap.pl" in the PATH using "export PATH=$PATH:~/gmes/"
- Finally, check that everything is working using "funannotate check"

### SNAP 
-It's recomended to install SNAP locally not through conda using:
git clone https://github.com/KorfLab/SNAP.git && cd SNAP/ && make &&  cp forge /opt/conda/envs/myenv/bin/

- Use "funannotate test -t predict", to ensure that all predictors are working.
   
## Gene prediction
  
````bash

funannotate predict -i Ahemp.f2.fasta.masked -s "Acropora hemprichii" -o funannotate_predict_Abdo --name Ahemp --rna_bam Ahemp_RNASeqAll.STAR.bam --stringtie Ahemp_RNASeqAll.Stringtie.gtf --protein_evidence uniprot_Acropora.faa.fasta --transcript_evidence Ahemp_RNASeqAll.transcripts.fasta  --organism other --busco_db metazoa --min_protlen 100 --cpus 50
````
- prediction can be updated by adding UTR, for this all RNASeq data were merged and running the following command:
````bash
funannotate update -i funannotate_predict_Abdo/ --species "Acropora hemprichii" -l Ahemp_RNASeqAll_1.fastq.gz -r Ahemp_RNASeqAll_2.fastq.gz --cpus 50

````
-We found 87 problematic gene models after the update that we dropped.

````bash
funannotate fix -i funannotate_final/update_results/Acropora_hemprichii.gbk -t funannotate_final/update_results/Acropora_hemprichii.tbl -d funannotate_final/drop_model.ls
````

### QC of the prediction


- get the prediction details from gff3 file

````bash
grep -v "#" funannotate_predict/predict_results/Acropora_hemprichii.gff3  | cut -f3 | sort | uniq -c
````
- BUSCO scores

````bash  
#against eukaryota_odb10
busco -i funannotate_predict/predict_results/Acropora_hemprichii.proteins.fa -m proteins -l eukaryota_odb10 -c 30 -o Ahemp_busco_eukaryota
#against metazoa_odb10
busco -i funannotate_predict/predict_results/Acropora_hemprichii.proteins.fa -m proteins -l metazoa_odb10 -c 30 -o Ahemp_busco_metazoa
````

## Functional annotation

-  transmembrane topology and signal peptide predictor
- install phobius
````bash
#download the program from:
https://phobius.sbc.su.se/data.html
#decompress
cat phobius101_linux.tgz| tar xz
#Edit phobius.pl L25 to my $DECODEANHMM = "$PHOBIUS_DIR/decodeanhmm.64bit"
#run the analysis
phobius/phobius.pl -short Acropora_hemprichii.proteins.fa > phobius.results.txt

````

- IterProScan and eggnog-mapper analyses will be computed separately
- IterProScan
- i used the funnannotate command
  
````bash
funannotate iprscan -i funannotate_final -m docker -c 30
````
- but also can be installed locally using:
  
````bash

mkdir interproscan

cd interproscan

wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.63-95.0/interproscan-5.63-95.0-64-bit.tar.gz

wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.63-95.0/interproscan-5.63-95.0-64-bit.tar.gz.md5

md5sum -c interproscan-5.63-95.0-64-bit.tar.gz.md5

tar -pxvzf  interproscan-5.63-95.0-64-bit.tar.gz
 
python3 setup.py -f interproscan.properties
````

- eggnog-mapper
````bash

mamba install -c bioconda -c conda-forge eggnog-mapper
#download_eggnog_data.py --data_dir /share/databases/
````

- interproscan and eggnog-mappe searches
  
````bash
mkdir Ahemp_funano_iprosc
/share/databases/interproscan/interproscan-5.63-95.0/interproscan.sh -t p --cpu 30 -goterms -pa -i funannotate_predict/predict_results/Acropora_hemprichii.proteins.fa -d Ahemp_funano_iprosc

emapper.py --cpu 30 -m diamond --data_dir /share/databases/eggnog/ -i funannotate_predict/predict_results/Acropora_hemprichii.proteins.fa -o Ahemp_eggnog
````

- Implement annotation using funannotate
````bash
funannotate-docker annotate -i funannotate_predict/ -s "Acropora hemprichii" -o funannotate_anno --busco_db  metazoa --eggnog  Ahemp_eggnog.emapper.annotations --iprscan Ahemp_funano_iprosc.xml --phobius phobius.results.txt  --cpus 40
````

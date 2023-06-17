# Ahemp
Acropora hemprichii genome structural and functional annotation
- The quality of the assembled genome was assessed using BlobToolKit (BTK), full handout for how to run it can be found here (https://github.com/blobtoolkit/tutorials/tree/main/futurelearn).Also, a course on BTK is avaliable as well https://www.futurelearn.com/courses/eukaryotic-genome-assembly-how-to-use-blobtoolkit-for-quality-assessment

- After first BTK check, ~90 contigs were removed because of low coverage or no-related taxa (Ahemp.gapclosed_f1.fasta).
- BTK re-runed with the filtered assembly and another 528 contigs were removed for low coverage or no-related taxa contigs (Ahemp.gapclosed_f2.fasta).

# Identifiy rRNA 

````bash
./barrnap -q -k euk Ahemp.gapclosed_f2.fasta --threads 50 --outseq Ahemp_rrna.fasta > Ahemp_rrna..gff 
````

## Identifying and maksing Repeats

- Repeats in Ahemp and avliable Acropora's coral genomes using RepeatModeler

````bash

singularity pull dfam-tetools-latest.sif docker://dfam/tetools:latest
singularity run dfam-tetools-latest.sif BuildDatabase -name Ahemp_genome Ahemp.gapclosed_f2.fasta
singularity run dfam-tetools-latest.sif RepeatModeler -database Ahemp_genome -LTRStruct -threads 40
````

- Repeats in avliable Acropora's coral genomes
````bash

for i in `ls *.fna|sed 's/_genomic.fna//g`;
do
    singularity run ../../dfam-tetools-latest.sif BuildDatabase -name $i ${i}_genomic.fna
    singularity run ../../dfam-tetools-latest.sif RepeatModeler -database $i -LTRStruct -threads 40;
done
````


- Repeats database assembly
cat *-families.fa > Acropora_RE_DB.fsa
unsearch ...

- Repeats masking in Ahemp using Acropora repeats DB
 
singularity run dfam-tetools-latest.sif RepeatMasker Ahemp.gapclosed_f2.fasta -lib Acropora_RE_DB.fsa -pa 8 -norna -nolow -xsmall

## How is the distribution of repeats by types

grep '>' Dehan101_genome-families.fa | sed -r 's/.+#//' | sed -r 's/\s+.+//' | sort | uniq -c


#RNA

STAR --runThreadN 50 --runMode genomeGenerate --genomeDir Ahemp_index --genomeFastaFiles
Ahemp.gapclosed_f2.fasta --genomeSAindexNbases 10

for i in `ls *.gz|sed 's/.fastq.gz//g'`;
do
    STAR --runThreadN 30 --genomeDir Ahemp_index --readFilesIn $i.fastq.gz --readFilesCommand "gunzip -c" --outSAMtype  BAM SortedByCoordinate --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix $i --limitBAMsortRAM 10000000000;
done

samtools merge Ahemp_RNASeqAll.STAR.bam KB_S01.Aligned.sortedByCoord.out.bam KB_S02.Aligned.sortedByCoord.out.bam KB_S03.Aligned.sortedByCoord.out.bam
samtools view -c
stringtie -p 30 -o Ahemp_RNASeqAll.Stringtie.gtf Ahemp_RNASeqAll.STAR.bam
grep -v "#" Ahemp_RNASeqAll.Stringtie.gtf  | cut -f3 | sort | uniq -c

gtf_genome_to_cdna_fasta.pl Ahemp_RNASeqAll.Stringtie.gtf Ahemp.gapclosed_f2.fasta > Ahemp_RNASeqAll.transcripts.fasta
gtf_to_alignment_gff3.pl Ahemp_RNASeqAll.Stringtie.gtf > Ahemp_RNASeqAll.Stringtie.gff3
TransDecoder.LongOrfs -t Ahemp_RNASeqAll.transcripts.fasta
diamond blastp -d /DBS/uniref90  -q Ahemp_RNASeqAll_longest_orfs.pep --max-target-seqs 1 --outfmt 6 --evalue 1e-5 --threads 30 > Ahemp_RNASeqAll_longest_orfs.out

TransDecoder.Predict -t Ahemp_RNASeqAll.transcripts.fasta --retain_blastp_hits Ahemp_RNASeqAll_longest_orfs.out
cdna_alignment_orf_to_genome_orf.pl Ahemp_RNASeqAll.transcripts.fasta.transdecoder.gff3 Ahemp_RNASeqAll.Stringtie.gff3 Ahemp_RNASeqAll.transcripts.fasta > Ahemp.StringtieTransdecoder.gff3


TSEBRA/bin/rename_gtf.py --gtf Dqua_braker/augustus.hints_utr.gtf --prefix Dqua --translation_tab translation.tab --out tsebra_result_renamed.gtf
gffread Dqua_renamed.gtf -g Dqua_filtered.fna -w Dqua1.braker.mRNA.fasta
gffread Dqua_renamed.gtf -g Dqua_filtered.fna -x Dqua1.braker.CDS.fasta
gffread Dqua_renamed.gtf -g Dqua_filtered.fna -y Dqua1.braker.protein.fasta

diamond blastp -q Dqua1.braker.protein.fasta -d /DBS/uniref90 --outfmt 6 --sensitive --max-target-seqs 1 --evalue 1e-25 --out Dqua1.braker.protein.dmd.out --threads


#funannotate

~/funannotate-docker predict -i Dqua_filtered.fna -s "Desmodesmus quadricauda" -o funannotate_predict --name Dqua --rna_bam DquaRNASeqAll.Dqua1.sort.bam --stringtie DquaRNASeqAll.Dqua1.gtf --protein_evidence uniprot_Scenedesmaceae.fasta --transcript_evidence DquaRNASeqAll.transcripts.fasta  --cpus 30

#others 
funannotate predict -i GCA_001630525.fna.masked -s "Coelastrella sp M60" --augustus_species Scenedesmaceae -o funannotate_GCA_001630525  --name CspM60 --protein_evidence uniprot_Scenedesmaceae.fasta --header_length 50 --cpus 30


funannotate remote -m phobius -e abdoallah.sharaf@umbr.cas.cz -i funannotate_predict/ -o funannotate_phobius

iterproscan and emapper separatly

~/funannotate-docker annotate -i funannotate_predict/ -s "Desmodesmus quadricauda" -o funannotate_anno   --eggnog  Dqua_eggnog.emapper.annotations --iprscan Dqua_funano_iprosc.xml --phobius funannotate_predict/annotate_misc/phobius.results.txt  --cpus 40


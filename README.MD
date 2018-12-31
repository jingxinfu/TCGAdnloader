## TCGA 

> The most important changes between HG19 and HG38 are the immune genes

### Hg19

#### RNA-Seq 

There are several kinds of RNA-Seq level 3 result available on [firebrowse](http://firebrowse.org/#). They are:

- .junction_quantification.txt
- .rsem.genes.results
- .rsem.isoforms.results
- .rsem.genes.normalized_results
- .rsem.isoforms.normalized_results
- .bt.exon_quantification.txt

Among those results, I only focus on gene level results. It's the `.rsem.genes.normalized_results` and `.rsem.genes.results`.

The transcript expression level are quantified by two different method in the `.rsem.genes.results`. For each samples, we can see two columns herer, which are **raw_count** and **scaled_estimate**. **Normlized count** result is stored in the `.rsem.genes.results` table.

> NOTE: At gene level, expression is simply the sum of the expression of all possible isoforms. [well explained](https://gitlab.com/zyxue/understanding-firebrowse-data-format/blob/master/confirm-relationship-between-gene-level-and-isoform-level-scaled-estimates.ipynb)

- **raw_count**: The estimated number of fragments derived from a given gene. (This value is not an integer because RSEM only reports a guess of how many reads belong to a transcript/gene)

- **normalized count**: The raw count divided by the 75th percentile of the column (after removing 0) and then multiply by 1000. 

  - > Most TCGA papers use the log-transformed "normalized_count" in their analyses.

- **scale_estimate**: The fraction of transcripts made up by a given gene.

  - > **Li B, et al., *Bioinformatics*. 2009**
    >
    > 1. For isofrom $i$, we will denote the fraction of transcripts as $\tau_i$, and  $\tau_i = \frac{v_i}{l_i}(\sum_j \frac{v_j}{l_j})^{-1}$where $l_i$ is the length, which is transcripts per million **(TPM)** by multiplying by $10^6$ . 
    > 2. The $\frac{c_i}{N}$ approaches $v_i$ as $N \to \infin$, where $c_i$ is the number of reads from isoform $i$ and $N$ is the total number of sequence reads.
    > 3. Reads per killobase per million mapped reads (**RPKM**; Mortazavi et al., 2008): $10^9\times \frac{v_i}{l_i}$, which is an unnormalized value of $\tau_i$

$TPM=\frac{count_i}{N\cdot l_i}(\sum_j \frac{count_j}{N\cdot  l_j})^{-1} \cdot 10^6$

#### Copy Number Variation

##### [GISTIC (Mermel CH, et al., *Genome Biol*. 2011)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3218867/) 

*On the notion that regions containing driver events should be altered at higher frequencies than regions containing only passengers, this tool identifies likely **driver SCNAs (Somatic copy-number alterations)** by evaluating the frequency and amplitude of observed events.*

> **Challenges:**
>
> 1. The abundance of SCNAs in humam cancers varies accoridng to their size. Chromosome-arm length SCNAs occurrs much more fequently than SCNAs of smaller size. ([Beroukhim R, et al.m *Nature*. 2010.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2826709/))
> 2. Above observation indicates that analysis methods need to model complex cancer genomes with SCNA types occurring at distinct background rates.







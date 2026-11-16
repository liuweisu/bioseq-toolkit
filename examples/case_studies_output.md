# bioseq-toolkit case studies

## 1. Bacterial genome assembly from short reads

18000 reads of 100 bp, 0.5% substitution errors, 30x coverage
k-mer spectrum: 348184 k-mers seen once (errors), peak near 17x (genomic)
k=21: 4 contigs, N50 17,424 bp, NG50 17,424 bp, total 60,080 bp (1.2s)
k=31: 1 contigs, N50 60,030 bp, NG50 60,030 bp, total 60,030 bp (1.2s)
k=51: 8 contigs, N50 11,400 bp, NG50 11,400 bp, total 60,332 bp (1.1s)

## 2. Plant genome: LTR retrotransposons and DNA transposons

| finder | precision | recall | exact boundaries |
| --- | --- | --- | --- |
| LTR | 1.00 | 0.93 | 73/80 |
| TIR | 0.74 | 0.61 | 35/80 |

LTR insertion ages (seed 7), molecular clock r = 1.3e-8 subst/site/year:
    39877-44277   RLC_Tekay   LTR/Copia  LTR identity 0.993  ~0.26 Myr
    52926-57326   RLC_Tekay   LTR/Copia  LTR identity 0.996  ~0.17 Myr
    64031-69931   RLG_Retand  LTR/Gypsy  LTR identity 0.991  ~0.33 Myr
    76080-81980   RLG_Retand  LTR/Gypsy  LTR identity 0.999  ~0.05 Myr
    85185-89585   RLC_Tekay   LTR/Copia  LTR identity 0.996  ~0.17 Myr
    99318-103718  RLC_Tekay   LTR/Copia  LTR identity 0.974  ~1.04 Myr
   115637-121537  RLG_Retand  LTR/Gypsy  LTR identity 0.984  ~0.61 Myr

## 3. Human DNA: a new Alu-like insertion in a patient's gene

Insertion of 318 bp at reference position 2887 (truth 2888), TSD ATTTTTGTTTTTCT (14 bp)
Location: exon 3; frameshift: False
Best library match: AluLike (SINE/Alu), identity 0.96, poly(A) tail: True
Interpretation: an exonic SINE insertion with a TSD and poly(A) tail is the signature of L1-mediated retrotransposition; in a real patient this would be confirmed by PCR across both junctions.

## 4. Bacterial DNA: origin of replication and genes

Cumulative GC-skew minimum at 29,996 (true oriC 30,000); maximum at 91,244 (true ter 90,000)
DnaA boxes (<=1 mismatch) within 500 bp of the predicted origin: 5 vs 0.30 expected (17x enrichment)
ORFs >= 100 aa: 126; planted genes recovered: 59/60 (precision 0.47: random sequence also contains long ORFs, which is why real gene finders model codon usage and ribosome binding sites)

## 5. Proteomics: peptide mass fingerprinting

24 MALDI peaks searched against 200 proteins, 20 ppm
  PROT0113: 9 peptides, coverage 42%, score 8.958
  PROT0189: 2 peptides, coverage 13%, score 1.969
  PROT0115: 2 peptides, coverage 5%, score 1.951
Top hit is the true protein: True. MW 34.0 kDa, pI 9.92, GRAVY -0.227

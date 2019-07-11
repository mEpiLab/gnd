# Introduction

Current culture-based methods for measuring *Escherichia coli* (*E. coli*) diversity are slow and laborious and often fail to identify *E. coli* types present at low abundance in communities dominated by one or two other distinct clones. Routinely the 16S rRNA gene is targeted using amplicon sequencing methods to provide information on the overall
bacterial community associated with a sample, but this method fails to provide any intra-species resolution. Thus, 16S rRNA community profiling can provide information on the relative abundance of *E. coli* but provides no information on the diversity of *E. coli* present in the sample, nor whether it is *E. coli* or Shigella spp. present as neither can be separated using the 16S rRNA gene.

## *gnd* has a hypervariable nucleotide sequence

The *E. coli* genome is a highly impacted by horizontal gene transfer leading to a mosaic structure of prophage insertions, acquisition and loss of genetic traits, and genome rearrangement through recombination events. A region of significant recombination in the *E. coli* genome is the O-antigen biosynthesis gene cluster (O-AGC) that encodes the
O-antigen component of lipopolysaccharide, an antigenic extracellular component anchored within the outer membrane of the *E. coli* cell wall. The immunogenic nature of the O-antigen causes purifying selection of various *E. coli* in the gut environment resulting in extensive antigenic variability of the O-antigen due to recombination. At least 184 antigenically distinct *E. coli* O-groups have been identified, many of which are also genetically distinct through recombination events where O-group gene cluster exchange has taken place. *gnd* encodes 6-phophogluconate dehydrogenase, an enzyme that catalyses the third stage of the pentose phosphate pathway and is found in all *E. coli* and many Enterobacteriaceae. Although playing no role in O-antigen biosynthesis, *gnd* is often found close by to the O-AGC and has a hypervariable nucleotide sequence brought about by recombination events impacting the O-antigen. As a result, *gnd* has been described as a 'hitch hiker' of recombination.

## *E. coli* community analysis using *gnd* metabarcoding

We have targeted the *gnd* gene as a suitable *E. coli* barcode gene that provides enhanced resolution to identify *E. coli* community diversity using culture-independent methods (Cookson et al. 2017; <https://www.nature.com/articles/s41598-017-00890-6>), akin to targeting 16S rRNA for bacterial community analysis . Like the 16S rRNA gene, the
complete *gnd* gene at 1407bp is unable to be used as a target for routine amplicon sequencing methods using the Illumina platform. Therefore, by using whole genome sequence data to compare *gnd* sequences from hundreds of *E. coli* isolates, we have been able to identify unique PCR primer sites that permit the amplification of a
284bp region of *gnd* allowing the differentiation of many *E. coli*. Furthermore, the use of WGS data available on public repositories and described in studies to evaluate the genomic diversity of *E. coli*, has provided valuable resources from which we have been able to assess the hypervariability of the *gnd* gene and develop a database of unique
*gnd* partial amplicons. Previous studies employing sequencing of the full O-AGC region have targeted the *wzx*/*wzy* or *wzm*/*wzt* gene
pairs located within the O-AGC to provide targets for O-group-specific PCR amplification and the identification of novel O-groups. However as
additional WGS data from *E. coli* becomes available, subsequent
analysis of O-AGC regions has provided evidence of genetically distinct
O-AGC, *wzx*/*wzy* or *wzm*/*wzt* gene pairs and unique *gnd* sequence
variants. For example, to extend our preliminary analysis of *gnd*
allelic variation obtained from WGS data submitted to GenBank, IMG and
PATRIC databases, we also examined genomes from 2244 *E. coli* isolates
consisting of a species-wide collection assembled for the purposes of
examining *E. coli* genomic diversity. This collection included multiple
animal and environmental isolates and isolates from six different
phylogroups (A, B1, B2, D, E, F). To date, 534 unique 284bp partial
*gnd* sequences identified from *E. coli* have been included in the
*gnd* database together with corresponding *E. coli* O-group/serotype
data and a representative genome or nucleotide accession number. O and
H-groups were identified from WGS *in silico* using SerotypeFinder
(<https://cge.cbs.dtu.dk/services/SerotypeFinder/>). OXY (*wzx*/*wzy*)
and OMT (*wzm*/*wzt*) are provisional O-group designations given to the
O-AGC region of untyped (&lt;95% similar at the nucleotide level) *E.
coli* not in the serotyping panel.

Differentiation of *gnd* SNP variants

Unlike 16S rRNA amplicon sequencing where operational taxonomic units
are arbitrarily clustered at the 97% homology level to represent
distinct ‘species’, and to negate the effect of sequencing errors
associated with the Illumina platform, 284bp partial *gnd* sequences
representative of *E. coli* isolates sometimes differ by only one
nucleotide/single nucleotide polymorphism. Thus, the resolution required
for *E. coli* community analysis to distinguish genuine *gnd* sequence
types that differ by only one nucleotide from those generated from
sequencing error has necessitated the development of an ‘error
correction’ model where a probability threshold can be established that
provides a mechanism where the relative abundance of two *gnd* sequence
types that differ by a single nucleotide, can be compared to determine
whether they are likely to be genuine or not. For example, two *gnd*
sequence types that are equally abundant in a single sample and that
differ by only a single nucleotide are more likely to be genuine than
two *gnd* sequence types that differ by only a single nucleotide where
one is represented by a high number of reads, and the other is
represented by only a small number of reads. Application of the error
correction method in parallel with the inclusion of ‘mock’ *E. coli*
communities with distinct *gnd* sequence types, allows many of the less
abundant *gnd* sequences types to be clustered with a ‘parent’ type.

Enhanced resolution of *E. coli* populations using *gnd* amplicon
sequencing methods

This *E. coli* community analysis technique is widely applicable to
understand the impact of various interventions on *E. coli* populations
associated with animal or environmental samples. DNA sample extractions
or sample enrichments provide a source of genetic material from which
the 284bp partial *gnd* sequence region from all *E. coli* can be
amplified with indexed PCR primers to generate a unique amplicon. By
using unique combinations of indexed forward and reverse primers with
Illumina tags, unique *gnd* amplicons from separate samples can then be
pooled and sequenced using the Illumina MiSeq platform. Subsequent
demultiplexing, pairing of forward and reverse sequence reads and
comparison with the *gnd* database provides information on the overall
community diversity and read numbers associated with each *gnd* sequence
type. The development of such a metabarcoding technique has provided a
much greater resolution of *E. coli* community diversity to be
established from complex matrices such as faecal material where
frequently, an animal may be colonised by &gt;10 distinct *E. coli*,
where one or two *gnd* sequence types are present in a sample together
with a larger number of sub-dominant *gnd* sequence types.

Target primer sites

Our observations suggest that the hypervariability of the *gnd* locus
and lack of conserved regions precludes the opportunity to develop PCR
primers that do not contain degeneracies. Current primers 2*gnd*F
5′-TCYATYATGCCWGGYGGVCAGAAAGAAG (*gnd* coordinates 415 to 442) and
2*gnd*R 5′-CATCAACCARGTAKTTACCSTCTTCATC (*gnd* coordinates 754 to 726)
generate a 284bp amplicon, excluding primer sequences. Allelic
variability outside of this 340bp total amplicon results in some
redundancy where *E. coli* with distinct O-AGC cannot be distinguished.
However, other *E. coli* O-groups (O37, O111, O113) or even serotypes
(O157:H7), in the current iteration of the database, possess unique
*gnd* sequence types.

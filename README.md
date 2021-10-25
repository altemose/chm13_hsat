#chm13_hsat\n
Assembly_HSat2and3_v2.pl\n
Nicolas Altemose 2021\n
  purpose: to annotate HSat2 and HSat3 arrays in a set of input DNA sequences\n
  inputs: a fasta file containing DNA sequences, and two provided text files containing HSat2 and HSat3 specific 24-mers\n
  output: a BED9 file listing all contiguous regions in the reference likely to be HSat2 or HSat3, along with their strand orientation\n
  note: merges all adjacent regions within 5 kb (same strand, same type)\n
  HSat2_kmers.txt and HSat3_kmers.txt must be in the same directory in which this script is executed\n
    the kmers in these input files were defined using HSat2/3 HuRef reads identified in Altemose et al. PLoS Comp Bio 2014\n
  usage: perl Assembly_HSat2and3_v2.pl /path/to/reference.fasta\n
  runtime and memory usage: annotates the full chm13 assembly in about 2 minutes on a macbook pro, with a memory footprint around ~1 GB\n

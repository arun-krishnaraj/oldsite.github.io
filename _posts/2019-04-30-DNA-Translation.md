---
title: "Computational Biology: Translation"
date: 2019-04-30
tags: [biology, python, applied data science]
header:
  image: "/images/"
excerpt: "Computational Biology, DNA Translation and Processing"
---


DNA is a self replicating material that is present in almost all living cells, and serves as a set of hereditary instructions. It's used to direct synthesis of proteins that eventually contribute to organism-level body structure or function. In this post we'll be exploring how DNA is used to make these proteins and then come up with some code to replicate the process.

DNA generally exists in a double-helix, with parallel strands running antiparallel to each other. Each DNA (deoxyribonucleic acid) strand is comprised of ribose sugars, nucleic acids, and phosphate groups. Each strand has a 5' and 3' end, referring to the numbered carbon in the ribose sugar.

Our main interest is the nucleic acids, since their order determines the information the DNA sequence holds. There are 4 nucleic acids in DNA, Adenine, Thymine, Cytosine, and Guanine. With 4 possible nucleic acids in each position along a strand, the number of possible sequences is 4^n, where n is the number of bases in the sequence. Each base pairs with another base on the opposing strand to form the ladder-rung structure. Bases pair under the following rules: A-T, and C-G (also A-U when talking about transcription)

Compare this to binary, which can generate 2^n possible sequences with base 0 and 1, where n is the number of bits. We can see that a byte (8 bits) has 256 configurations, while an 8 base DNA sequence has 256^2 configurations. From this we can appreciate how efficient evolution was in selecting a data storage system.

The first step in extracting this information is Transcription. DNA is generally unable to leave the nucleus, or too big to be efficiently used to manufacture proteins directly. So an intermediate mRNA is created by reading the sequence from 3'-5' along one strand and synthesizing mRNA from the 5' to 3'.

mRNA is a messenger form of RNA that has a few key differences from DNA.
1. It is single stranded
2. It replaces the base Thymine with Uracil
3. It's often not the entire length of the sequence.

I'll touch more on the third point later on, but for now we have all the information we need to code a DNA transcription mechanism.

```python
blank = str()
x: str
for x in string:
  if x == "A" or x == "a":
      blank = blank + "U"
  elif x == "T" or x == "t":
      blank = blank + "A"
  elif x == "C" or x == "c":
      blank = blank + "G"
  elif x == "G" or x == "g":
      blank = blank + "C"
  else:
      print('Error, invalid base pair')
```

This code captures the basic process of DNA transcription. By inputting a string of bases in sequence, we can generate the resulting mRNA sequence. Notice that we haven't yet specified the directionality of this sequence, so our code is assuming its runs 3' to 5' from left to right; we will add this flexibility later on.

Now that we have an mRNA sequence, we can move outside the nucleus and begin talking about translation. Translation occurs in the ribosomes, where an mRNA strand is used to specify the addition of amino acids in a sequence, called a polypeptide chain. Amino acids are the fundamental building blocks of proteins, and bond in a linear manner in the polypeptide chain. There are typically 20 amino acids in eukaryotic organisms.

Each of these amino acids is coded for by one or more 3 base sequences on the mRNA, called codons. For example the mRNA sequence UUU would direct the ribosomal machinery to add a Phenylalanine to the current polypeptide chain.

Using this knowledge, we can code this mechanism like so:

```python
for x in range(0, len(blank)+1, 3):
    codon = (seq1[x-3:x])
    if codon == "UUU" or codon == "UUC":
        amino_seq = amino_seq + "Phe"
    elif codon == "UUA" or codon == "UUG" or codon == "CUU" or codon == "CUC" or codon == "CUA" or codon == "CUG":
        amino_seq = amino_seq + "Leu"
    elif codon == "AUG":
        amino_seq = amino_seq + "Met"
```

This code iterates through the mRNA sequence we generated previously by groups of 3 characters. These 3 characters can then be compared to codons for each amino acid (three of which I have listed) and used to append the polypeptide chain amino_seq. After repeating this process for all 20 amino acids and all 64 codons we should be able to reliably generate amino acid sequences from DNA base sequences right? Well not exactly.

Our model so far is ignoring a key feature of the 3 base codon system: reading frames. A reading frame refers to which index we begin this iteration through. For example, in the string "AAAUUUGGG", with codon size 4 we can generate the following codon sequences:
1. By indexing at position 0, "AAA"+"UUU"+"GGG"
2. By indexing at position 1, "AAU"+"UUG"+*"GG"*
3. By indexing at position 2, "AUU"+"UGG"+*"G"*
4. By indexing at position 3, "UUU"+"GGG"

There's three takeaways from this feature that we need to update our code to reflect.
* With codon size n, the sequence we get from indexing at position k+N is equivalent to the sequence by indexing at k, where k is a nonnegative integer, and N is some positive multiple of n, but excludes the first N/n codons
* When indexing some sequence of size j, n-1 possible reading frames will include partial or remainder sequences that cannot be used
* Most importantly, the 3 possible codon sequences we can generate are generally not the same, and can produce vastly different polypeptide chains.

Let's look at what polypeptide chains our mRNA sequence "AAAUUUGGG" could produce:
    1. "Lysine" + "Phenylalanine" + "Glycine"
    2. "Asparagine" + "Leucine"
    3. "Isoleucine" + "Tryptophan"
    4. "Phenylalanine" + "Glycine"

This small mRNA sequence is capable of producing 3 completely different polypeptide chains, and this variation grows rapidly as we start using larger sequence data.

We can account for this in our mechanism by creating 3 uniquely indexed sequences within our for loop and iterating through each one to generate 3 separate polypeptide chains.

```python
seq1 = blank
seq2 = blank[1:]
seq3 = blank[2:]
for x in range(0, len(blank)+1, 3):
    codon1 = (seq1[x-3:x])
      .....
    codon2 = (seq2[x - 3:x])
      .....
    codon3 = (seq3[x - 3:x])
      .....
```

Where the code for each codon index is just the above translation mechanism with a separate amino acid sequence for each index.

Now that we have a basic mechanism for going from DNA base sequences to polypeptide chains, we can go back and refine our model to account for some additional complexity.

As mentioned above, the entire DNA sequence is unfeasible to transcribe into a single mRNA. Instead, DNA sequences are comprised of two different regions:
1. Exons: regions that will be included in the final mRNA
2. Introns: regions that will be removed from the mRNA before it moves to the ribosome

For humans, about 50% of our the genome is comprised of intronic regions, meaning that they are usually not included in any mRNA and subsequently not used to manufacture any proteins. For our model to efficiently process realistic sequence data, we want to exclude these intronic regions from the sequence prior to running the translation mechanism.

**Exon-Intron Recognition**

Intron positions can be determined experimentally, by using known DNA sequences and observing the resulting mRNA or polypeptide sequence. These intron positions and sequences can then be aggregated in a database for future use.

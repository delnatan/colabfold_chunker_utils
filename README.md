# colabfold_chunker_utils
various scripts to break up a long polypeptide sequence into overlapping sequence. For use with colabfold and ChimeraX for stitching the segments back together.

# About this repo
Inspired by this ChimeraX recipe (https://rbvi.github.io/chimerax-recipes/big_alphafold/bigalpha.html), I wanted to work out a way to "stitch" together large proteins (that can't fit in one GPU when using Colabfold/Alphafold). The first one is a python program "chunker.py" that takes in protein sequences in FASTA format and breaks them up into overlapping chunks. The result is another fasta file where the protein has been divided into multiple overlapping segments. The default segment length is 1400 amino acids with 200 amino-acid overlaps. 1400-a.a. is small enough to fit in a 3080 Ti GPU so you can run it using a local install of colabfold (https://github.com/YoshitakaMo/localcolabfold). Each segment's name has a '_seg\d_' in it where the '\d' is the fragment number ("my-protein_seg1.pdb", "my-protein_seg2.pdb", ...). I modified Tom's `bigalpha` command so that each segment is actually connected. This way, when you straighten out some loops, the connection between segments are preserved.

1) Run `colabfold_batch` with the "chunked" FASTA file to get your protein segments.
2) Open "bigalpha_v2.py" with ChimeraX to register the command "bigalpha"
3) Open "straighten.py" with ChimeraX to register the command "straighten"
4) run `bigalpha my-protein directory /path/to/segmented_pdbs`
5) run `straighten` on disordered residues. For example, you can do `straighten #1:21-30` to straighten out the residues 21-30 in model #1.



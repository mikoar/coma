# COMA

COMA (**C**ross-correlation **O**ptical **M**ap **A**lignment) is a software for aligning optical mapping sequences in
CMAP format. Optical mapping is a technique for imaging DNA molecules along which specific labels are captured. These
labels form distinct patterns along DNA molecules based on their nucleotide sequences. Optical mapping has been widely
applied in assisted-scaffolding in sequence assemblies, microbial strain typing and detection of structural variations.

## Features

- Supports input files in CMAP format from OpGen and BioNano platforms
- Uses cross-correlation to find approximate alignment of sequences
- Allows for flexible alignment parameters such as gap penalty and mismatch penalty
- Outputs alignment results in XMAP format

## Installation

To install the software, you need to have Python 3.9 or higher installed on your system. You can download the software
using the following commands:

```
git clone https://github.com/mikoar/coma
cd ./coma/
python3 -m venv env
source env/bin/activate
pip install -e .
coma --help
```

## Usage

To use the software, you need to provide two input files: a reference file and a query file. The reference file contains
the optical mapping sequences of the reference genome, and the query file contains the optical mapping sequences of the
sample genome. The files should be in CMAP format, with each line containing the sequence ID, label position, label
channel and label intensity. For example:

```
# CMAP File Version:	0.1
# Label Channels:	1
# Nickase Recognition Site 1:	GCTCTTC
# Number of Consensus Maps:	2
#h CMapId	ContigLength	NumSites	SiteID	LabelChannel	Position	StdDev	Coverage	Occurrence
#f int	float	int	int	int	float	float	float	float
1	1000000.0	100	1	1	10000.0	0.0	1.0	1.0
1	1000000.0	100	2	1	20000.0	0.0	1.0	1.0
...
2	800000.0	80	1	1	12000.0	0.0	1.0	1.0
2	800000.0	80	2	1	24000.0	0.0	1.0	1.0
...
```

The software also accepts optional arguments to customize the alignment process. You can use `-h` or `--help` to see the
full list of arguments and their descriptions.

To run the software with default parameters and [sample data](https://bionano.com/public-datasets/) , you can use the
following command:

`
coma -r ./data/NA12878_BSPQI/alignmolvref_contig24_r.cmap -q ./data/NA12878_BSPQI/alignmolvref_contig24_q.cmap -o ./alignment24.xmap
`

## Indels detection

All scripts used during indels detection are stored in the "sv" folder. There are two possible paths used during this
process. The first one analyzes places where the molecules were joined in the "joined" output version files.

The `molecule_indels.py` wrapper present in this repository includes alignment of molecules using the comma "all" option and then analyzes its output files looking for indels:

`
python molecule_indels.py -r reference.cmap -q query.cmap -j aligned_all.xmap -o joined_indels.txt
`

The second possible path focuses on the conflicts observed when more than one segment is created for a single molecule
and they are joined to form one alignment. To identify indels in those places use created `segments_indels.py` wrapper
which includes aligning molecules with COMA:

`
python segment_indels.py -r reference.cmap -q query.cmap -a aligned_output.xmap -s found_segments_output.csv -o segment_indels.txt
`

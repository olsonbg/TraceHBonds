TraceHBonds
=====

A program for finding strings of hydrogen bonded atoms in a trajectory file
generated from the _Discover_ and _LAMMPS_ molecular dynamics program.

[![Build Status](https://travis-ci.com/olsonbg/TraceHBonds.svg?branch=master)](https://travis-ci.com/olsonbg/TraceHBonds)

![Hydrogen bond strings, colored by chain length](/images/HydrogenBondStringsB.png?raw=true "Hydrogen Bond strings")

# Contents

  * [Installation](#installation)
  * [Usage](#usage)
  * [Input](#input)
      * [Discover](#discover)
      * [LAMMPS](#lammps)
  * [Output](#output)
      * [Hydrogen bond strings](#sizehist)
          * [Generating Images](#generating-images)
      * [Neighbor distance in chains](#neighborhist)
      * [Hydrogen bond lengths](#lengths)
      * [Hydrogen bond angles](#angles)
      * [Hydrogen bond list](#list)
      * [Hydrogen bond lifetime correlations](#lifetime)

# Installation

To compile the program

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=release ..
make
src/TraceHBonds
```

To generate the documentation as html, which will be available at
`docs/html/index.html`, use

```bash
make docs
```

To cross-compile, use your toolchain cmake file as:

```bash
cmake -DCMAKE_BUILD_TYPE=release -DCMAKE_TOOLCHAIN_FILE=<Toolchain cmake file> ..
```

# Usage

A typical command line would look like this when using CAR/MDF files:

```bash
TraceHBonds --input molecule.arc -p HBonds -s .dat -H h1o -A o2h -r 2.5 -a 90.0 --verbose --all
```

or this when using _LAMMPS_ trajectory files:

```bash
TraceHBonds --input molecule.data --trajectory molecule.lammpstrj --molecule molecule.dat -p HBonds -s .dat -H h1o -A o2h -r 2.5 -a 90.0 --verbose --all
```

Below is a table of all options available from the command line. The long
form options are preceded with `--`, and short form a single `-`. For
options with both a long and short form, either one may be used on the
command line.



|Long form                                                   | Short form |Option Type    | Required? | Description |
|:-----------------------------------------------------------|:----------:|:--------------|:---------:|:------------|
|input                          <a name="input-t"></a>       |   i        | string        | yes       | The archive file generated from _Discover_ (without .arc), or the _LAMMPS_ data file with the extension. |
|[trajectory](#trajectory)      <a name="trajectory-t"></a>  |   t        | string        | no        | The trajectory file generated from _LAMMPS_, required when load _LAMMPS_ data. |
|[molecule](#molecule)          <a name="molecule-t"></a>    |   m        | string        | no        | The molecule file used for defining molecules in LMMPS data file. |
|outprefix                                                   |   p        | string        | yes       | All output will have this string as a prefix to the filenames. For example, to save data as `HBonds1.dat`, use `-p HBonds` as the prefix|
|outsuffix                                                   |   s        | string        | yes       | All output will have this string as a suffix to the filenames. For example, to save data as 'HBonds1.dat', use `-s .dat` as the suffix|
|rcutoff                                                     |   r        | real number   | yes       | Set the cutoff length, in angstroms, for the determination of a hydrogen bond (e.g. `-r 2.5`). |
|anglecutoff                                                 |   a        | real number   | yes       | Set the cutoff angle, in degrees, for the determination of a hydrogen bond (e.g. `-a 90.0`).|
|hydrogen                       <a name="hydrogen-t"></a>    |   H        | string        | yes       | Set the force field of donor hydrogens for hydrogen bonding (e.g. `-H h1o`). More than one force field may be used by specifying this option multiple times.  **NOTE** the short option is a capital 'H.'|
|acceptor                       <a name="acceptor-t"></a>    |   A        | string        | yes       | Set the force field of acceptor atoms for hydrogen bonding. More than one force field may be used by specifying this option multiple times (e.g. `-A o2h -A o1=`). **NOTE** the short option is a capital 'A.'|
|bins                                                        |   b        | integer       | no        | Minimum number of bins to show in histograms (e.g. `-b 20`).|
|[povray](#povray)              <a name="povray-t"></a>      |            |               | no        | Output in povray format, relevant for [--sizehist](#sizehist-t) only.|
|[json](#blender)               <a name="json-t"></a>        |            |               | no        | Output in json format, relevant for [--sizehist](#sizehist-t) only. Useful for processing with python, and blender scripts in the blender/ directory|
|jsonall                                                     |            |               | no        | Saves the chemical structure for all frames in json format. Useful for processing with python, and blender scripts in the blender/ directory.|
|incell                                                      |            |               | no        | Apply PBC to all hydrogen bond chains (each chain will start inside the PBC cell). Relevant for [--sizehist](#sizehist-t) only.|
|verbose                                                     |            |               | no        | Show verbose messages while running. |
|brief                                                       |            |               | no        | Show brief messages while running. |
|[lifetime](#lifetime)          <a name="lifetime-t"></a>    |            |               | no        | Calculate hydrogen bond lifetime correlations. |
|[lengths](#lengths)            <a name="lengths-t"></a>     |            |               | no        | Save length of all hydrogen bonds.|
|[angles](#angles)              <a name="angles-t"></a>      |            |               | no        | Save angle of all hydrogen bonds. |
|[list](#list)                  <a name="list-t"></a>        |            |               | no        | Save list of all hydrogen bonds. |
|[sizehist](#sizehist)          <a name="sizehist-t"></a>    |            |               | no        | Save hydrogen bond strings and histograms. |
|[neighborhist](#neighborhist)  <a name="neighborhist-t"></a>|            |               | no        | Save neighbor length lists. |
|all                                                         |            |               | no        | Do all calculations and save all data. |
|help                                                        |   h        |               | no        | This help screen |

# Input

The standard _CAR/MDF_ and _LAMMPS_ trajectory files are recognized. If
compiled with support for _LZMA_, _GZIP_, and/or _BZIP2_, compressed files
may be read and written.

## Discover

The standard _CAR/MDF_ files are read.

## LAMMPS

### Data file (--input)

The [--input](#input-t) option, when loading a _LAMMPS_ file, defines the name
of the data file. The masses section of the data file must have a comment,
it will be used as the forefield type used for atom selection with the
[--hydrogen](#hydrogen-t) and [--acceptor](#acceptor-t) flags.

And example masses section of a _LAMMPS_ data file:

~~~~~~~~~~~~~
Masses

   1  15.999400 # o
   2  12.011150 # c2
   3  12.011150 # c
   4   1.007970 # h
   5  12.011150 # cPrime
   6  12.011150 # c3
   7  15.999400 # oPrime
   8  15.999400 # oh
   9   1.007970 # ho
~~~~~~~~~~~~~

### <a name="trajectory"></a>Trajectory (--trajectory)

The [--trajectory](#trajectory-t) option defines the name of a _LAMMPS_
trajectory file. The trajectory file must be saved by a dump custom or
custom/gz command with atom attributes starting with 'id mol x y z', any
other atom attribute may follow.

### <a name="molecule"></a>Molecule (--molecule)

The [--molecule](#molecule-t) option is used when loading _LAMMPS_ trajectory
and data files, and loads a file that defines which atoms belong to which
molecule. While the _LAMMPS_ data file has a column, `Molecule-Id` for
defining molecules, the `msi2lmp` script uses residue IDs there. This file
is required when using _LAMMPS_ files, whether `msi2lmp` is used or not.

The example file below defines 4 molecules and their associated atoms. The
atoms must be in sequence for each molecule.

~~~~~~~~~~~~~
# Molecule FirstAtom LastAtom
     1      1    248
     2    249    496
     3    497    744
     4    745    992
~~~~~~~~~~~~~

# Output

Depending upon the command line options, many data files may be generated.
The description of all calculations are listed below

  - [Hydrogen bond strings (--sizehist)](#sizehist)
    - [List of Files Created](#sizehist-files)
  - [Neighbor distance in chains (--neighborhist)](#neighborhist)
    - [List of Files Created](#neighborhist-files)
  - [Hydrogen bond lengths (--length)](#lengths)
    - [List of Files Created](#lengths-files)
  - [Hydrogen bond angles (--angles)](#angles)
    - [List of Files Created](#angles-files)
  - [Hydrogen bond list (--list)](#list)
    - [List of Files Created](#list-files)
  - [Hydrogen bond lifetime correlations (--lifetime)](#lifetime)
    - [List of Files Created](#lifetime-files)

## <a name="sizehist"></a>Hydrogen bond strings (--sizehist)

The [--sizehist](#sizehist-t) option calculates hydrogen bonds, traces the
hydrogen bonds into connected strings, and tabulates the sizes.

### <a name="sizehist-files"></a>Files Created

 - \<prefix\>#\<suffix\>

Where _#_ indicates the frame of the trajectory. For a trajectory containing
100 frames, 100 files will be generated.

### Description of files

The files generated from the [--sizehist](#sizehist-t) option consist of two
parts: the [individual chains](#individual-chains) and their atoms, and
[histograms of chain sizes](#chain-histograms).

#### Individual Chains

As an example, a few lines from a data file follows, showing 3 chains.  The
output will look a little different if either the [--povray](#povray-t), or
[--json](#json-t) options are used.

~~~~~~~~~~~~~
# Current Element : 589
# Atoms in Chain : 3
# Molecules : 2
# Unique forcefields : 3
# Times chain switched between Molecules (switching) : 1
# Periodic boundary conditions applied.
  21.6276   27.0261   44.9056 [O]  Molecule  O87  o2h
  21.3585   27.9405   44.8265 [H]  Molecule  H592  h1o
  20.1321   29.5776   43.8281 [O]  Molecule2  O1493  o1=
# Chain end-to-end distance: 3.147660


# Current Element : 590
# Atoms in Chain : 9
# Molecules : 2
# Unique forcefields : 3
# Times chain switched between Molecules (switching) : 1
# Periodic boundary conditions applied.
  24.6835   23.4992   41.1162 [O]  Molecule  O95  o2h
  24.2312   22.9546   41.6786 [H]  Molecule  H601  h1o
  23.0622   23.9573   43.4825 [O]  Molecule  O88  o2h
  22.4585   23.9183   44.2600 [H]  Molecule  H593  h1o
  21.5617   24.1262   45.5633 [O]  Molecule9  O8207  o2h
  21.6362   24.7840   46.3132 [H]  Molecule9  H8711  h1o
  22.6499   25.4540   47.5703 [O]  Molecule9  O8208  o2h
  23.5335   25.7281   47.6493 [H]  Molecule9  H8712  h1o
  25.0893   24.1132   46.6706 [O]  Molecule9  O8209  o1=
# Chain end-to-end distance: 5.602967


# Current Element : 591
# Atoms in Chain : 5
# Molecules : 2
# Unique forcefields : 2
# Times chain switched between Molecules (switching) : 1
# Periodic boundary conditions applied.
  28.8885   22.9984   40.5006 [O]  Molecule  O96  o2h
  29.6141   22.3548   40.3725 [H]  Molecule  H602  h1o
  30.7792   21.1779   38.8075 [O]  Molecule6  O5419  o2h
  31.6526   21.1038   39.2359 [H]  Molecule6  H5925  h1o
  33.3777   21.2051   38.4486 [O]  Molecule6  O5420  o2h
# Chain end-to-end distance: 5.251614
~~~~~~~~~~~~~

#### Chain Histograms

This histogram shows how many chains of a specific length there are in this
frame. It also shows how many chains form closed loops (begin and end at the
same hydrogen bond). This particular data had no closed loops.

```
# Atoms/HBonds |Count| (For all Chains, including Closed Loops)
#    3 /   1   |  175|**************************************************************
#    5 /   2   |   72|**************************
#    7 /   3   |   19|*******
#    9 /   4   |   16|******
#   11 /   5   |    7|**
#   13 /   6   |    1|
#   15 /   7   |    1|
#   17 /   8   |    1|
#
# Atoms/HBonds |Count| (For Closed Loops)
#
```

The next series of histograms show how many times a hydrogen bond chain
switches to another molecule, for each chain length. Switching of 0 means
the chain was on a single molecule. Switching of 1 means it started on a
molecule and ended on another. The switching number does not indicate how
many molecules a chain is composed of, since it may switch back and forth
between two molecules multiple times.

```
# Switching |Count| (For Chain length of 3)
#    0      |  107|**************************************************************
#    1      |   68|***************************************
#
# Switching |Count| (For Chain length of 5)
#    0      |   18|************************
#    1      |   46|**************************************************************
#    2      |    8|***********
#
# Switching |Count| (For Chain length of 7)
#    0      |    8|**************************************************************
#    1      |    6|***********************************************
#    2      |    5|***************************************
```

The following histograms shows how many molecules a chain is composed of,
for each chain length.

```
# Molecules |Count| (For Chain length of 3)
#    1      |  107|**************************************************************
#    2      |   68|***************************************
#
# Molecules |Count| (For Chain length of 5)
#    1      |   18|**********************
#    2      |   50|**************************************************************
#    3      |    4|*****
#
# Molecules |Count| (For Chain length of 7)
#    1      |    8|**************************************************************
#    2      |    8|**************************************************************
#    3      |    3|***********************
```

### Generating Images

#### <a name="blender"></a>Using Blender

The program [Blender](http://www.blender.org) can be used to convert the
output of [--sizehist](#sizehist-t) [--json](#json-t) to an image or even a
movie of the image rotating. Here is an example image using the
[hbchain.py](blender/hbchain.py) script in the [blender](blender/)
directory (I need to add a legend):

![Figure 1: Hydrogen bond strings, colored by chain length. Rendered with Blender](/images/HydrogenBondStringsB.png?raw=true "Hydrogen Bond strings (Bender)")

The image was created by Blender using

```bash
blender -P blender/hbchain.py -- -f output/HBonds1.json
```

where `output/HBonds1.json` is the output of `TraceHBonds` with the appropriate
options.

For help with the Blender script, try

```bash
blender -P hbchain.py -- --help
```

#### <a name="povray"></a>Using POV-Ray

The program [POV-Ray](http://www.povray.org) can be used to convert the
output of [--sizehist](#sizehist-t) [--povray](#povray-t) to an image or
even a movie of the image rotating. Here is an example image using the
[prettybox.pov](povray/prettybox.pov) script in the [povray](povray/)
directory:

![Figure 2: Hydrogen bond strings, colored by chain length. Rendered with POV-Ray.](/images/HydrogenBondStrings.png?raw=true "Hydrogen Bond strings (POV-Ray")


The legend text was added with [GIMP](http://www.gimp.org).


## <a name="neighborhist"></a>Neighbor Distance in chains (--neighborhist)

The [--neighborhist](#neighborhist-t) option calculates the distance between
non-hydrogen atoms in the hydrogen bond chains. For a chain consisting of
only oxygen and hydrogen atoms, this would calculate the neighbor distances
between oxygen atoms in the chain.

### <a name="neighborhist-files"></a>Files Created

  - \<prefix\>-NN-AllFrames\<suffix\>
  - \<prefix\>-NN-Combined\<suffix\>
  - \<prefix\>-NN-only\<suffix\>

All files contain tab delimited text.

### Description of files

The data start with statistics of each individual frame in a trajectory
[NN-AllFrames](#nn-allframes), then all frames together
[NN-Combined](#nn-combined), followed by combining all chain lengths
together for a list of only neighbor distances [NN-only](#nn-only).

#### NN-AllFrames

The table in this file consists of a Count, Average, and StdDev column for
each frame of the trajectory, so the table can have many columns if a
trajectory with many frames is used.  It's not uncommon to read 100, 1000,
or more frames. Sample output is shown below, in table format for easy
viewing.  **Note**: Only the first 2 off 1000 frames of a large trajectory
file are shown, and the first 13 of 27 atoms in chain for this particular
sample trajectory file.

n.n. = Nearest neighbor, and f is the frame number.

|Atoms in chain | Nth n.n. |Count(f=0) |Average(f=0)  |StdDev(f=0)   |Count(f=1) |Average(f=1)  |StdDev(f=1)   |
|:--------------|:--------:|:----:|--------:|--------:|:----:|--------:|--------:|
|3              |    1     |  175 |2.76054  |0.193364 |181   |2.74879  |0.189819 |
|               |          |      |         |         |      |         |         |
|5              |    1     |  144 |2.75077  |0.197673 |154   |2.76424  |0.204004 |
|               |    2     |  72  |4.75806  |0.635903 |77    |4.64895  |0.732237 |
|               |          |      |         |         |      |         |         |
|7              |    1     |  57  |2.7841   |0.201227 |60    |2.75164  |0.22737  |
|               |    2     |  38  |5.01798  |0.581104 |40    |4.78955  |0.619429 |
|               |    3     |  19  |6.88486  |0.948441 |20    |6.28423  |1.4411   |
|               |          |      |         |         |      |         |         |
|9              |    1     |  64  |2.7419   |0.195794 |52    |2.75299  |0.217351 |
|               |    2     |  48  |4.70959  |0.743358 |39    |4.64655  |0.670621 |
|               |    3     |  32  |6.49354  |1.32184  |26    |6.21771  |1.27282  |
|               |    4     |  16  |8.07276  |1.88673  |13    |7.66016  |2.06277  |
|               |          |      |         |         |      |         |         |
|11             |    1     |  35  |2.69268  |0.130133 |25    |2.69679  |0.159238 |
|               |    2     |  28  |4.6182   |0.587221 |20    |4.7146   |0.696855 |
|               |    3     |  21  |6.43492  |0.672441 |15    |6.60273  |1.07617  |
|               |    4     |  14  |8.2164   |1.09446  |10    |8.3314   |1.70064  |
|               |    5     |  7   |9.65632  |1.52155  |5     |9.95754  |2.07065  |
|               |          |      |         |         |      |         |         |
|13             |    1     |  6   |2.76977  |0.150354 |12    |2.78778  |0.272727 |
|               |    2     |  5   |4.82386  |0.386727 |10    |4.65379  |0.546978 |
|               |    3     |  4   |6.49485  |0.499925 |8     |6.08429  |0.413477 |
|               |    4     |  3   |8.00401  |0.466091 |6     |7.38176  |1.13275  |
|               |    5     |  2   |9.38291  |0        |4     |8.29956  |2.52054  |
|               |    6     |  1   |9.95807  |0        |2     |9.34289  |0        |

#### NN-Combined

The table in this file combines data from all the frames into single
Count, Average, and StdDev columns. Sample data is shown below, in table
format for easy viewing.  **Note**: Only the first 13 of 27 atoms in a chain
are shown for this particular sample trajectory file.

|Atoms in chain | Nth n.n. | Count|  Average|StdDev  |
|:--------------|:--------:|:----:|--------:|-------:|
|3              |    1     |22655 |2.86888  |0.159418|
|               |          |      |         |        |
|5              |    1     |22230 |2.86208  |0.159733|
|               |    2     |11115 |4.68886  |0.70739 |
|               |          |      |         |        |
|7              |    1     |19182 |2.85872  |0.159745|
|               |    2     |12788 |4.67906  |0.703396|
|               |    3     |6394  |6.07321  |1.31757 |
|               |          |      |         |        |
|9              |    1     |15888 |2.85741  |0.158921|
|               |    2     |11916 |4.68648  |0.702248|
|               |    3     |7944  |6.08516  |1.26822 |
|               |    4     |3972  |7.19271  |1.89908 |
|               |          |      |         |        |
|11             |    1     |12900 |2.85904  |0.159156|
|               |    2     |10320 |4.68125  |0.705725|
|               |    3     |7740  |6.08836  |1.27751 |
|               |    4     |5160  |7.21112  |1.88283 |
|               |    5     |2580  |8.16827  |2.45298 |
|               |          |      |         |        |
|13             |    1     |10452 |2.85698  |0.15854 |
|               |    2     |8710  |4.67941  |0.690471|
|               |    3     |6968  |6.06397  |1.27943 |
|               |    4     |5226  |7.13685  |1.88829 |
|               |    5     |3484  |8.06878  |2.4206  |
|               |    6     |1742  |8.9083   |2.88292 |

#### NN-only

This file combines all frames as in [NN-Combined](#nn-combined), and also
all the 'Atoms in chain' column for a complete nearest neighbor table.
Sample output is shown below, in table format for easy viewing. Only first 6
of 13 shown for this particular sample trajectory file.

|Nth n.n.   | Count | Average |StdDev  |
|:----------|:-----:|--------:|-------:|
|1          |145342 |2.85969  |0.159263|
|2          |92595  |4.68331  |0.700193|
|3          |62503  |6.08867  |1.27689 |
|4          |43526  |7.22511  |1.85633 |
|5          |30943  |8.2102   |2.36345 |
|6          |22332  |9.11149  |2.79441 |

## <a name="lengths"></a>Hydrogen bond lengths (--lengths)

The [--lengths](#lengths-t) option calculates the length of all hydrogen
bonds (hydrogen-acceptor distance), in every frame.

### <a name="lengths-files"></a>File Created

  - \<prefix\>-lengths\<suffix\>

### Description of file

Single column of data listing the hydrogen bond lengths in angstroms.

## <a name="angles"></a>Hydrogen bond angles (--angles)

The [--angles](#angles-t) option calculates the angle of all hydrogen bonds,
in every frame.

### <a name="angles-files"></a>File Created

  - \<prefix\>-angles\<suffix\>

### Description of file

Single column of data listing the hydrogen bond angles in degrees.

## <a name="list"></a>Hydrogen bond list (--list)

The [--list](#list-t) option calculates the list of all hydrogen bonds,
in every frame.

### <a name="list-files"></a>File Created

  - \<prefix\>-list\<suffix\>

### Description of file

Seven (7) column data, in tab delimited format. The columns are:

 1. Frame number, starting with 0.
 1. Molecule of atom connected to donor hydrogen.
 1. Name of atom connected to donor hydrogen.
 1. Molecule of donor hydrogen atom.
 1. Name of donor hydrogen atom.
 1. Molecule of acceptor atom.
 1. Name of acceptor atom.

## <a name="lifetime"></a> Hydrogen bond lifetime correlations (--lifetime)

The [--lifetime](#lifetime-t) option calculates the continuous and
intermittent hydrogen bond lifetime autocorrelation.

### <a name="lifetime-files"></a>File Created

  - \<prefix\>-lifetimes\<suffix\>

### Description of file

Three column data, in tab delimited format. The columns are:

 1. Frame number, starting with 0.
 1. Continuous hydrogen bond lifetime correlation
 1. Intermittent hydrogen bond lifetime correlation


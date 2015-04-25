# TraceHBonds

###USAGE
        TraceHBonds -i <arc file> -p <prefix> -s <suffix> -r <distance
        cutoff> -a <angle cutoff> -H <hydrogen forcefield> -A <acceptor
        forcefield> [-b <number>] [--povray]

###OPTIONS
        --input <arc file>
        -i <arc file>
                <Arc file> is the archive file generated from Discover.

        --outprefix <prefix>
        -p <prefix>
                <Prefix> of the output filename. The text before the
                integer in the filename. For a filename of 'HBonds1.dat'
                the <prefix> would be HBonds'

        --outsuffix <suffix>
        -s <suffix>
                <Suffix> of the output filename. The text after the integer
                in the filename. For a filename of 'HBonds1.dat' the
                <suffix> would be '.dat'

        --rcutoff <Rc>
        -r <Rc>
                Set the cutoff length to <Rc> for the determination of a
                hydrogen bond.

        --anglecutoff <Ac>
        -a <Ac>
                Set the cutoff angle to <Ac> for the determination of a
                hydrogen bond.

        --hydrogen <force field>
        -H <force field>
                Set the <force field> of donor hydrogens for hydrogen
                bonding (e.g. -H h1o). More than one <force field> may be
                used by specifying additional -H <force field> parameters.
                NOTE: the short option is a capital 'H.'

        --acceptor <force field>
        -A <force field>
                Set the <force field> of acceptor atoms for hydrogen
                bonding. More than one <force field> may be used by
                specifying additional -A <force field> parameters (e.g. -A
                o2h -A o1=). NOTE: the short option is a capital 'A.'

        --bins <number>
        -b <number>
                Minimum <number> of bins to show in histograms.

        --povray
                Output in povray format.

        --verbose
                Show some extra information while running.

        -h
                This help screen



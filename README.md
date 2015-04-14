# TraceHBonds

###USAGE 
        TraceHBonds -p <prefix> -s <suffix> -f <number> -l <number> -P
                    <prefix> -S <suffix> [-b <number>] [--povray]

###OPTIONS
        --inprefix <prefix>
        -p <prefix>
                <Prefix> of the input filename. The text before the integer
                in the filename. For a filename of 'HBonds1.dat' the
                <prefix> would be 'HBonds'

        --insuffix <suffix>
        -s <suffix>
                <Suffix> of the input filename. The text after the integer
                in the filename. For a filename of 'HBonds1.dat' the
                <suffix> would be '.dat'

        --outprefix <prefix>
        -P <prefix>
                <Prefix> of the output filename. The text before the
                integer in the filename. For a filename of 'HBonds1.txt'
                the <prefix> would be HBonds'

        --outsuffix suffix
        -S suffix
                <Suffix> of the output filename. The text after the integer
                in the filename. For a filename of 'HBonds1.txt' the
                <suffix> would be '.dat'

        --first <number>
        -f <number>
                The first <number> to start the processing at. <number> is
                the integer in the filename. For filenames of
                'HBonds1.dat,' 'HBonds2.dat,' ..., 'HBonds1000.dat,' the
                first <nmber> would be 1.

        --last <number>
        -l <number>
                The last <number> to start the processing at. <number> is
                the integer in the filename. For a filenames of
                'HBonds1.dat,' 'HBonds2.dat,' ..., 'HBodns1000.dat,' the
                last <number> would be 1000.

        --bins <number>
        -b <number>
                Minimum <number> of bins to show in histograms.

        --povray
                Output in povray format.

        --verbose
                Show some extra information while running.

        -h
                This help screen



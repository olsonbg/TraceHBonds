/**
 * \file
 * \author Brian G. Olson
 * \date   30 April 2015
 * \brief  Make histograms of chain length, molecules, and correlations
 *
 **/
#ifndef _Histograms_h
#define _Histograms_h

#include "ListOfHBonds.h"
#include "queue.h"
#include "WorkerThreads.h"
#include "cpu.h"
// Macros

/** Vector of unsigned integers */
typedef std::vector<unsigned int> vui;
/** Vector of vector of unsigned integers */
typedef std::vector< vui > vvui;
/** Vector of doubles */
typedef std::vector< double > vd;
/** Vector of vector of doubles */
typedef std::vector< vd > vvd;


/**
 * metafunction to allocate a vector
 *
 * Allocates vector with \p nelem elements of \p val. If \p v is already \p
 * nelem or larger, the original vector is returned. If \p v is smaller than \p
 * nelem, initialize the needed number of elements to \p val. Values already in
 * \p v are never touched.
 *
 * \param[in,out] v      vector
 * \param[in]     val    elements to initialize \p v with
 * \param[in]     nelem  number of elements to initialize \p v
 *
 * \tparam        T      Type of vector.
 *
 * \return \c TRUE on success
 * \return \c FALSE on failure.
 */
template<class T> bool alloc_vector( std::vector<T> *v,
                                     T val,
                                     unsigned int nelem);

/**
 * metafunction to allocate a vector of vectors
 *
 * Allocates vector with \p melem vectors of \p nelem elements of \p val. For
 * each melem, alloc_vector(std::vector<T>*, T, unsigned int) is called with \p
 * nelem.
 *
 * Values already in \p v are never touched.
 *
 * @param[in,out] v      vector of vectors
 * @param[in]     val    element to initialize v with
 * @param[in]     nelem  number of elements to initialize vector
 * @param[in]     melem  number of elements to initialize vector of vector
 *
 * \tparam        T      Type of vector
 *
 * \return TRUE on success
 * \return FALSE on failure
 */
template<class T> bool alloc_vector(std::vector< std::vector<T> > *v,
                                    T val,
                                    unsigned int nelem,
                                    unsigned int melem);

/**
 * Structure to hold data in histogram form
 */
struct Histograms_s
{
	unsigned int TrjIdx; /**< Trajectory index (Frame number) */

	/**
	 * \name Chain lengths
	 */
	/**@{*/
	vui ChainLength;     /**< Chain Lengths */
	vui ClosedLoop;      /**< Chain Lengths of only Closed Loops */

	unsigned int MaxChainLength; /**< Maximum Chain Length */
	unsigned int MaxClosedLoop;  /**< Maximum Chain length of closed loops */
	/**@}*/

	/**
	 * \name Molecules in Chains
	 */
	/**@{*/
	vvui SwitchesInChain;    /**< Molecule Switches for each Chain Length */
	vvui MoleculesInChain;   /**< Molecules in Chain, for each Chain Length */
	vui MaxSwitchesInChain;  /**< Maximum molecule switches */
	vui MaxMoleculesInChain; /**< Maximum number of molecules in chain **/
	/**@}*/

	/**
	 * \name Nearest Neighbors
	 */
	/**@{*/
	/**
	 * NearestNeighbors (NN):
	 *   - 1st element, NN[0]:vector of nearest neighbor distances.
	 *   - 2nd element, NN[1]:vector of next nearest neighbor distances.
	 *   - 3rd element, NN[2]:vector of nextnext nearest neighbor distances.
	 *   - ...
	 *   - nth element, NN[n]:vector of nth nearest neighbor distances.
	 */
	std::vector< std::vector<double> >NearestNeighbors;
	/**@}*/
};

/**
 * Go through the vector of hydrogen bond strings and:
 *  - Bin Chain Lengths                             (1D)
 *  - Bin Chain Lengths of only Closed Loops        (1D)
 *  - Bin Molecule Switches for each Chain Length   (2D)
 *  - Bin Molecules in Chain, for each Chain Length (2D)
 *
 * \param[in] HBStrings   Hydrogen bond strings
 * \param[in] TrjIdx      Trajectory index to operate on (Frame number)
 *
 * \return Histograms as #Histograms_s struct.
 */
struct Histograms_s
makeHistograms( std::vector<ListOfHBonds *>* HBStrings,
                unsigned int TrjIdx);

/**
 * Generate neighbor distances of non hydrogen atoms.
 *
 * Build up the Histograms_s::NearestNeighbors vector with neighbor distances
 * between all non-hydrogen atoms in hydrogen bond strings (\p HBStrings), for
 * the frame corresponding to Histograms_s::TrjIdx. HBStrings of all lengths
 * are considered together.
 *
 * \param[in,out] Histograms  Histograms
 * \param[in]     HBStrings   List of hydrogren bond strings
 * \param[in]     Cell        Dimension of periodic cell
 */
void getNeighbors( struct Histograms_s *Histograms,
                   std::vector<ListOfHBonds *> *HBStrings,
                   struct PBC *Cell);
/**
 * Print out histograms in various formats.
 *
 * The default output format is plain text, however, the output of hydrogen
 * bond strings may be in povray format if \p POVRAY is TRUE. If \p NumBins is
 * non-zero, then each histogram will print a maximum of \p NumBins bins,
 * otherwise each histogram will output all their bins. With \p NumBins zero,
 * each histogram will be of different length (bins), therefore to aid in
 * subsequent data processing, all histogram can be sent to \p NumBins bins.
 * This will lead to many bins containing data of zero, but all histograms will
 * be of the same number of bins. The ouput consists of
 *
 *   - Hybrogen bond strings with coordinate and type of each constituent atom.
 *   - Histogram of string lengths
 *   - Histogram of closed loop lengths
 *   - Histogram of number of molecule switches
 *     - One Histogram for each chain length
 *   - Histogram of Molecules in chains
 *     - One histogram for each chain length
 *
 * \todo Check the order of output.
 *
 * \param[in] out         Stream to send results to.
 * \param[in] HBStrings   Hydrogen bond strings
 * \param[in] Histogram   Histograms
 * \param[in] CC          String to use for comments in output.
 * \param[in] NumBins     Number of bins to output.
 * \param[in] Cell        Dimension of periodic cell
 * \param[in] TrjIdx      Trajectory index to operate on (Frame number)
 * \param[in] POVRAY      Output in POVRAY format
 *
 */
void
prntHistograms( std::ostream *out,
                std::vector<ListOfHBonds *> *HBStrings,
                struct Histograms_s *Histogram,
                std::string CC, unsigned int NumBins,
                struct PBC *Cell, unsigned int TrjIdx,
                bool POVRAY);
/**
 * Add to the counts in element \p bin of vector \p h, making sure enough space
 * is allocated for \p h by calling
 * alloc_vector(std::vector<T>*, T, unsigned int).  Also record the maximum bin
 * as \p max noting that it may have already been assigned a value from a
 * previous call.
 *
 * \param[in,out] h    Vector for histogram
 * \param[in,out] max  Maximum bins in historam
 * \param[in]     bin  Bin of histogram to add a single count to.
 */
bool Bin(vui *h, unsigned int *max, unsigned int bin);

/**
 * Add to the counts in element \p h[\p hb][\p c], making sure enough space is
 * allocated for \p h by calling
 * alloc_vector(std::vector< std::vector<T> > *, T, unsigned int,unsigned int).
 * and for \p hmax by calling
 * alloc_vector(std::vector< std::vector<T> > *, T, unsigned int,unsigned int).
 * Also, record the maximum \p c for each \p hb in \p hmax, noting that it may
 * have already been assigned a value from a previous call.
 *
 * \param[in,out] h     Vector of vectors for histogram
 * \param[in,out] hmax  Maximum \p c for each \p hb bin
 * \param[in]     hb    Bin corresponding to hydrogen bond
 * \param[in]      c    Bin corresponding to chain length
 *
 */
bool Bin(vvui *h, vui *hmax, unsigned int hb, unsigned int c);

/** \todo Correlations(), CorrelationsThread(), and CorrelationsTableThread()
 * should be moves to their own file.
 */

/**
 * Calculate autocorrelation.
 *
 * Puts jobs for both CorrelationsTableThread() and CorrelationsThread() in the Queue, and
 * combines results from all jobs.
 *
 * \param[in] out   Stream to send results to.
 * \param[in] v     Boolean indicating which hydrogen bonds are formed, and
 *                  in which frame
 */
void Correlations( std::ostream *out,
                   std::vector< std::vector<bool> > *v );

/**
 * Combine autocorrelation results for all hydrogen bonds.
 *
 * For both continuous and intermittent hydrogen bond autocorrelations, sum
 * results for all hydrogen bonds in a frame (time slice), then divide by
 * number total number.
 *
 * \param[out] C            Average continuous hydrogen bond autocorrelation
 * \param[out] I            Average intermittent hydrogen bond autocorrelation
 * \param[in]  continuous   Continuous autocorrelation of each hydrogen bond,
 *                          from CorrelationsTableThread()
 * \param[in]  intermittent Intemittent autocorrelation of each hydrogen bond,
 *                          from CorrelationsTableThread()
 * \param[in]  NumThreads   Number of jobs this calculation has been split
 *                          into, determined by NumberOfCPUs().
 * \param[in]  ThreadID     Job ID of this call 0 <=\p ThreadID <\p NumThreads.
 *
 */
void CorrelationsThread(vd *C, vd *I,
                        vvui *continuous, vvui *intermittent,
                        unsigned int NumThreads, unsigned int ThreadID );

/**
 * Calculate autocorrelation of hydrogen bonds.
 *
 * Calculate both the continuous and intermittent hydrogen bond
 * autocorrelations for each hydrogen bond in the system. Uses a sliding
 * window, therefore \p fcutoff should be half the total number of frames.
 *
 * The average over all hydrogen bonds in the system can be obtained by
 * subsequently calling CorrelationsThread().
 *
 * \param[in]  v            Boolean indicating which hydrogen bonds are formed,
 *                          and in which frame
 * \param[out] continuous   Continuous hydrogen bond autocorrelation
 * \param[out] intermittent Intermittent hydrogen bond autocorrelation
 * \param[in]  numHBs       Number of hydrogen bonds.
 * \param[in]  fcutoff      Time cutoff (in number of frame units).
 * \param[in]  NumThreads   Number of jobs this calculation has been split
 *                          into, determined by NumberOfCPUs().
 * \param[in]  ThreadID     Job ID of this call 0 <=\p ThreadID <\p NumThreads.
 *
 */
void CorrelationsTableThread( std::vector< std::vector<bool> > *v,
                              vvui *continuous, vvui *intermittent,
                              unsigned int numHBs,
                              unsigned int fcutoff,
                              unsigned int NumThreads,
                              unsigned int ThreadID);

#endif // _Histograms_h

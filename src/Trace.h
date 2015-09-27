/**
 * \file
 * \brief Trace connections between hydrogen bond atoms
 */
#ifndef _Trace_h
#define _Trace_h

extern bool THB_VERBOSE;

/**
 *
 * Find all chains of hydrogen bonds that are chemically linked together
 * (including hydrogen bonds)
 *
 * Called as a thread, put in the Queue by Thread()
 *
 * \todo This function always returns \c TRUE, therefore change to \c void
 *
 * \param[in,out] HBStrings   Hydrogen bonds which connect to each other
 * \param[in]     HBit        Iterator for all hydrogen bonds in a frame
 *
 * \return \c TRUE
 */
bool TraceThread( std::vector<ListOfHBonds *> *HBStrings,
                  struct HydrogenBondIterator_s *HBit);

/**
 * Puts NumberOfCPUs() TraceThread() jobs in Queue.
 *
 * \param[in,out] HBStrings    Hydrogen bonds which connect to each other
 * \param[in]     TrjIdx_iter  Iterator for all hydrogen bonds in a frame
 *
 * \return \c TRUE
 */
void Trace( std::vector<ListOfHBonds *> *HBStrings,
            std::vector<struct HydrogenBondIterator_s> *TrjIdx_iter );
#endif // _Trace_h

#ifndef _tracehbonds_lzma_h
#define _tracehbonds_lzma_h
#include <boost/iostreams/char_traits.hpp>
#include <boost/iostreams/concepts.hpp>
#include <boost/iostreams/operations.hpp>
#include <boost/iostreams/read.hpp>
#include <boost/asio/basic_streambuf.hpp>
#include <iostream>
#include <lzma.h>

class lzma_input_filter : public boost::iostreams::input_filter {
	private:
		lzma_stream strm;
		uint8_t inbuf[BUFSIZ];
		uint8_t outbuf[BUFSIZ];
		lzma_ret retcode;
		size_t out_size;
		size_t out_curr;
		bool lzma_done;


		template<typename Source> bool getraw(Source &src)
		{
			if (strm.avail_in == 0 && !lzma_done)
			{
				int c;
				strm.next_in = inbuf;

				int i = 0;
				while( i < sizeof(inbuf) )
				{
					if ( (c = boost::iostreams::get(src)) == EOF ||
					     c == boost::iostreams::WOULD_BLOCK)
					{
						lzma_done = true;
						break;
					}

					inbuf[i] = static_cast<uint8_t>(std::string::traits_type::to_char_type(c));
					++i;
				}
				// for(i=0; (i< sizeof(inbuf)) && 
				//     ( ((c=boost::iostreams::get(src)) != EOF) || ( c == boost::iostreams::WOULD_BLOCK )); ++i)
				//     inbuf[i] = static_cast<uint8_t>(std::string::traits_type::to_char_type(c));

				strm.avail_in = i;

				return true;
			}
			return false;
		}
		template<typename Source> bool decompress(Source &src)
		{
			lzma_action action = LZMA_RUN;

			getraw(src);

			if ( lzma_done ) { action = LZMA_FINISH; }

			// std::cout << "strm.avail_in : " << strm.avail_in  <<  inbuf[0] << "\n";
			// std::cout << "strm.avail_out: " << strm.avail_out << outbuf[0] << "\n";

			lzma_ret retcode = lzma_code(&strm, action);

			// std::cout << "strm.avail_in : " << strm.avail_in  <<  inbuf[0] << "\n";
			// std::cout << "strm.avail_out: " << strm.avail_out << outbuf[0] << "\n";

			if ( strm.avail_in == 0 && strm.avail_out == sizeof(outbuf) )
				return false;

			if (strm.avail_out < sizeof(outbuf) || retcode == LZMA_STREAM_END) {
			// if ( strm.avail_out < sizeof(outbuf) ) {
				out_curr = 0;
				out_size = sizeof(outbuf) - strm.avail_out;

				// FILE *outfile = fopen("asd.out","w");
				// fwrite(outbuf, 1, out_size, outfile);
				// fclose(outfile);

				// std::string aa((char *)outbuf, out_size);
				// std::cout << ">" << aa << "]\n";
				// exit(1);

				return true;
			}

			std::cout << "done? " << lzma_done;
			if (retcode == LZMA_OK) {
			std::cout << "Why am I here?\n"; }
			//     size_t write_size = sizeof(outbuf) - strm->avail_out;

			//     if (fwrite(outbuf, 1, write_size, outfile)
			//             != write_size) {
			//         fprintf(stderr, "Write error: %s\n",
			//                 strerror(errno));
			//         return false;
			//     }

			//     strm->next_out = outbuf;
			//     strm->avail_out = sizeof(outbuf);

			if (retcode != LZMA_OK) {
				// Once everything has been decoded successfully, the
				// return value of lzma_code() will be LZMA_STREAM_END.
				//
				// It is important to check for LZMA_STREAM_END. Do not
				// assume that getting ret != LZMA_OK would mean that
				// everything has gone well or that when you aren't
				// getting more output it must have successfully
				// decoded everything.
				if (retcode == LZMA_STREAM_END)
				{
					return true;
				}
				// It's not LZMA_OK nor LZMA_STREAM_END,
				// so it must be an error code. See lzma/base.h
				// (src/liblzma/api/lzma/base.h in the source package
				// or e.g. /usr/include/lzma/base.h depending on the
				// install prefix) for the list and documentation of
				// possible values. Many values listen in lzma_ret
				// enumeration aren't possible in this example, but
				// can be made possible by enabling memory usage limit
				// or adding flags to the decoder initialization.
				const char *msg;
				switch (retcode) {
					case LZMA_MEM_ERROR:
						msg = "Memory allocation failed";
						break;

					case LZMA_FORMAT_ERROR:
						// .xz magic bytes weren't found.
						msg = "The input is not in the .xz format";
						break;

					case LZMA_OPTIONS_ERROR:
						// For example, the headers specify a filter
						// that isn't supported by this liblzma
						// version (or it hasn't been enabled when
						// building liblzma, but no-one sane does
						// that unless building liblzma for an
						// embedded system). Upgrading to a newer
						// liblzma might help.
						//
						// Note that it is unlikely that the file has
						// accidentally became corrupt if you get this
						// error. The integrity of the .xz headers is
						// always verified with a CRC32, so
						// unintentionally corrupt files can be
						// distinguished from unsupported files.
						msg = "Unsupported compression options";
						break;

					case LZMA_DATA_ERROR:
						msg = "Compressed file is corrupt";
						break;

					case LZMA_BUF_ERROR:
						// Typically this error means that a valid
						// file has got truncated, but it might also
						// be a damaged part in the file that makes
						// the decoder think the file is truncated.
						// If you prefer, you can use the same error
						// message for this as for LZMA_DATA_ERROR.
						msg = "Compressed file is truncated or "
							"otherwise corrupt";
						break;

					default:
						// This is most likely LZMA_PROG_ERROR.
						msg = "Unknown error, possibly a bug";
						break;
				}

				fprintf(stderr, "Decoder error: "
				        "%s (error code %u)\n",
				        msg, retcode);
				return false;
			}

			return true;
		}

	public:
		explicit lzma_input_filter(lzma_stream ini=LZMA_STREAM_INIT)
			: strm(ini), out_size(0), out_curr(0), lzma_done(false) {
			// strm = LZMA_STREAM_INIT;


			lzma_ret retcode = lzma_stream_decoder( &strm, UINT64_MAX, LZMA_CONCATENATED);

			if ( retcode == LZMA_OK)
			{
				strm.next_in = NULL;
				strm.avail_in = 0;
				strm.next_out = outbuf;
				strm.avail_out = sizeof(outbuf);
				return;
			}

			// Something went wrong. The possible errors are documented in
			// lzma/container.h (src/liblzma/api/lzma/container.h in the source
			// package or e.g. /usr/include/lzma/container.h depending on the
			// install prefix).
			//
			// Note that LZMA_MEMLIMIT_ERROR is never possible here. If you
			// specify a very tiny limit, the error will be delayed until
			// the first headers have been parsed by a call to lzma_code().
			const char *msg;
			switch (retcode) {
				case LZMA_MEM_ERROR:
					msg = "Memory allocation failed";
					break;

				case LZMA_OPTIONS_ERROR:
					msg = "Unsupported decompressor flags";
					break;

				default:
					// This is most likely LZMA_PROG_ERROR indicating a bug in
					// this program or in liblzma. It is inconvenient to have a
					// separate error message for errors that should be impossible
					// to occur, but knowing the error code is important for
					// debugging. That's why it is good to print the error code
					// at least when there is no good error message to show.
					msg = "Unknown error, possibly a bug";
					break;
			}

			fprintf(stderr, "Error initializing the decoder: %s (error code %u)\n",
			        msg, retcode);

		}

		template<typename Source> int get(Source& src)
		{
			if ( out_curr == out_size )
			{
				strm.next_out = outbuf;
				strm.avail_out = sizeof(outbuf);

				if ( !decompress(src) )
					return EOF;

				out_curr = 0;
			}

			return outbuf[out_curr++];
		}


		template<typename Source> void close(Source &) {
			lzma_end(&strm); }

		// template<typename Source>
		// static uint8_t read_uint8(Source& src, int error)
		// {
		//     int c;
		//     if ((c = boost::iostreams::get(src)) == EOF || c == WOULD_BLOCK)
		//         throw gzip_error(error);
		//     return static_cast<uint8_t>(traits_type::to_char_type(c));
		// }
};
#endif // _tracehbonds_lzma_h


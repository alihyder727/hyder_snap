// C/C++ headers
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <algorithm>
#include <cstring>
#include <string>
#include <vector>

// Athena++ headers
#include "../defs.hpp"
#include "../inversion/mcmc.hpp"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

void usage() {
  printf("Usage:  (1) fismerge -i inputfile1 inputfile2 ... -o outputfile\n");
  printf("        (2) fismerge -i *.fits -o outputfile\n");
  printf("\n");
  printf("Merge fits files generated from MCMC inversion\n");
  printf("\n");
  printf("Examples:\n");
  printf("\n");
  printf("fitsmerge -i mwr.out1.00001.fits mwr.out1.00002.fits -o mwr-0101a.fits\n");
  printf("fitsmerge -i mwr.out1.*.fits -o mwr-0101a.fits\n");
  printf("\n");
};

int main(int argc, char **argv) {
  int c, num = 0;
  char outfile[256];
  bool input = false;
  bool output = false;
  int status = 0;

#ifdef MPI_PARALLEL
  if (MPI_SUCCESS != MPI_Init(&argc, &argv)) exit(1);
#endif

  std::vector<std::string> infiles;

  while ((c = getopt (argc, argv, "i:o:")) != -1) {
    switch (c) {
      case 'i':
        infiles.push_back(optarg);
        input = true;
        break;
      case 'o':
        strcpy(outfile, optarg);
        output = true;
        break;
      case 'h':
        usage();
        exit(0);
      default:
        usage();
        exit(1);
    }
  }

  if (input && output) {
    for (int i = optind; i < argc; ++i)
      infiles.push_back(argv[i]);

    // sort files
    std::sort(infiles.begin(), infiles.end());

    mcmc_opts opts;
    mcmc_recs recs, new_recs;

#ifdef MPI_PARALLEL
    opts.mpi_comm = MPI_COMM_WORLD;
#endif

    // read first, determine dimension
    mcmc_load_fits(infiles[0].c_str(), &opts, &recs);
    mcmc_alloc(&new_recs, infiles.size()*recs.nstep, recs.nwalker, recs.ndim, recs.nvalue);
    mcmc_append_recs(&new_recs, &recs);

    // open input files, alloc = false
    for (int i = 1; i < infiles.size(); ++i) {
      mcmc_load_fits(infiles[i].c_str(), &opts, &recs, false);
      mcmc_append_recs(&new_recs, &recs);
    }

    // reset accepted states
    new_recs.accept = recs.accept;
    new_recs.reset = 0;

    // create output file, include_last = true
    mcmc_save_fits(outfile, &opts, &new_recs, true);

    mcmc_free(&recs);
    mcmc_free(&new_recs);
  } else {
    printf("ERROR: Input/Output files missing\n");
    printf("\n");
    usage();
    exit(1);
  }

#ifdef MPI_PARALLEL
  MPI_Finalize();
#endif
}

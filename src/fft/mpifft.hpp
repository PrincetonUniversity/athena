#ifdef MPI_PARALLEL
#if defined(FFTW_MPI)
#include "fftw3-mpi.h"

typedef fftw_complex AthenaFFTComplex;
typedef ptrdiff_t AthenaFFTInt;
typedef struct AthenaFFTPlan{
  fftw_plan plan;
  enum AthenaFFTDirection dir;
  int dim;
} AthenaFFTPlan;

#elif defined(PLIMPTON)
#include "mpi.h"
#include "fft_3d.h"
#include "fft_2d.h"
typedef fftw_complex AthenaFFTComplex;
typedef int AthenaFFTInt;
typedef struct AthenaFFTPlan{
  struct fft_plan_3d *plan3d;
  struct fft_plan_2d *plan2d;
  fftw_plan plan;
  int dir;
  int dim;
} AthenaFFTPlan;

#elif defined(ACCFFT)
#include "accfft.h"
typedef fftw_complex AthenaFFTComplex;
typedef int AthenaFFTInt;
typedef struct AthenaFFTPlan{
  accfft_plan *plan3d;
  fftw_plan plan;
  int dir;
  int dim;
} AthenaFFTPlan;

#endif
#else // MPI_PARALLEL
#include "fftw3.h"

typedef fftw_complex AthenaFFTComplex;
typedef int AthenaFFTInt;
typedef struct AthenaFFTPlan{
  fftw_plan plan;
  enum AthenaFFTDirection dir;
  int dim;
} AthenaFFTPlan;

#endif // MPI_PARALLEL

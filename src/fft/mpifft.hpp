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
#include "plimpton/fft_3d.h"
#include "plimpton/fft_2d.h"
typedef fftw_complex AthenaFFTComplex;
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
typedef struct AthenaFFTPlan{
  accfft_plan *plan3d;
  fftw_plan plan;
  int dir;
  int dim;
} AthenaFFTPlan;

#elif defined(FFTW)
#include "fftw3.h"

typedef fftw_complex AthenaFFTComplex;
typedef int AthenaFFTInt;
typedef struct AthenaFFTPlan{
  fftw_plan plan;
  enum AthenaFFTDirection dir;
  int dim;
} AthenaFFTPlan;
#else //no FFT
typedef Real AthenaFFTComplex[2];
typedef struct AthenaFFTPlan{
  void *plan;
  enum AthenaFFTDirection dir;
  int dim;
} AthenaFFTPlan;
#endif

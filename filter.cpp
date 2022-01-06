#include "macros.h"
#include "initialize.h"
#include <stdlib.h>


// constants for filtering (centered uniform)
// f_i(filtered) = c * f_i + b / 2 * (f_{ i - 1 }+f_{ i + 1 })
//               + a/2*(f_{i-2}+f_{i+2}) 
// One can derive from Taylor expansion
// a = (c - 1) / 3, b = -4(c - 1) / 3

constexpr real C_FILTER_2P_0 = 5./8.;
constexpr real C_FILTER_2P_1 = 1./4.;
constexpr real C_FILTER_2P_2 = -1./16.;

// constants for filtering (Lele 1992 C.2.11)
constexpr real C_FILTER_BOUND0_0 = 15. / 16.;
constexpr real C_FILTER_BOUND0_1 = 4. / 16.;
constexpr real C_FILTER_BOUND0_2 = -6. / 16.;
constexpr real C_FILTER_BOUND0_3 = 4. / 16.;
constexpr real C_FILTER_BOUND0_4 = -1. / 16.;

constexpr real C_FILTER_BOUND1_0 = 1. / 16.;
constexpr real C_FILTER_BOUND1_1 = 3. / 4.;
constexpr real C_FILTER_BOUND1_2 = 6. / 16.;
constexpr real C_FILTER_BOUND1_3 = -4. / 16.;
constexpr real C_FILTER_BOUND1_4 = 1. / 16.;


real *filter1d(real *arr, int nx)
{
	real* arr_f = (real *)malloc(sizeof(real) * nx);

	// arr_f[0] = arr[0];
	arr_f[0] = C_FILTER_BOUND0_0 * arr[0] +
		C_FILTER_BOUND0_1 * arr[1] +
		C_FILTER_BOUND0_2 * arr[2] +
		C_FILTER_BOUND0_3 * arr[3] +
		C_FILTER_BOUND0_4 * arr[4];

	// arr_f[1] = arr[1];
	arr_f[1] = C_FILTER_BOUND1_0 * arr[0] +
		C_FILTER_BOUND1_1 * arr[1] +
		C_FILTER_BOUND1_2 * arr[2] +
		C_FILTER_BOUND1_3 * arr[3] +
		C_FILTER_BOUND1_4 * arr[4];

	for (int i = 2; i < nx-2; i++)
	{
		arr_f[i] = C_FILTER_2P_0 * arr[i] +
			C_FILTER_2P_1 * (arr[i - 1] + arr[i + 1]) +
			C_FILTER_2P_2 * (arr[i - 2] + arr[i + 2]);
	}

	// arr_f[nx - 2] = arr[nx - 2];
	arr_f[nx - 2] = C_FILTER_BOUND1_0 * arr[nx - 1] +
		C_FILTER_BOUND1_1 * arr[nx - 2] +
		C_FILTER_BOUND1_2 * arr[nx - 3] +
		C_FILTER_BOUND1_3 * arr[nx - 4] +
		C_FILTER_BOUND1_4 * arr[nx - 5];

	// arr_f[nx - 1] = arr[nx - 1];
	arr_f[nx - 1] = C_FILTER_BOUND0_0 * arr[nx - 1] +
		C_FILTER_BOUND0_1 * arr[nx - 2] +
		C_FILTER_BOUND0_2 * arr[nx - 3] +
		C_FILTER_BOUND0_3 * arr[nx - 4] +
		C_FILTER_BOUND0_4 * arr[nx - 5];

	return arr_f;
}


void filter(real *uu, int nx, int nvar)
{
	real *arr_1d = (real *)malloc(sizeof(real) * nx);
	real *arr_1d_filtered;

	for (int iv = 0; iv < nvar; iv++)
	{
		for (int i = 0; i < nx; i++)
		{
			arr_1d[i] = uu[IDX(i, iv, nvar)];
		}

		arr_1d_filtered = filter1d(arr_1d, nx);


		for (int i = 0; i < nx; i++)
		{
			uu[IDX(i, iv, nvar)] = arr_1d_filtered[i];
		}

		free(arr_1d_filtered);
	}

	free(arr_1d);
}
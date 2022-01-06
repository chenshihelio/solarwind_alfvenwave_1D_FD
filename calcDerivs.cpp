#include <stdlib.h>
#include <string.h>
#include "initialize.h"
#include "macros.h"
#include "support.h"
#include <iostream>

using namespace std;

real* derivX(real *uu, real *xgrid, real *Br, int nx, int nvar)
{
	real *du;
	du = (real *)malloc(sizeof(real) * nx * nvar);
	memset(du, 0, sizeof(real) * nx * nvar);

	real *utmp;
	utmp = (real *)malloc(sizeof(real) * nvar);
	memset(utmp, 0, sizeof(real) * nvar);

	real *waveSpeeds;

	for (int i = 0; i < nx; i++)
	{
		
		
		if (i == 0)// left most point -- right scheme
		{
			for (int iv = 0; iv < nvar; iv++)
			{
				for (int j = 0; j < NPOINTFD; j++)
				{
					du[IDX(i, iv, nvar)] = du[IDX(i, iv, nvar)] + 
						FD_Coeff_Right_1st[i*NPOINTFD + j] * uu[IDX(i + j, iv, nvar)];
				}
				/*du[IDX(i, iv, nvar)] = (uu[IDX(i + 1, iv, nvar)] - uu[IDX(i, iv, nvar)])
					/ (xgrid[i + 1] - xgrid[i]);*/
			}
		}
		else if (i == 1)// 2nd point -- centerR scheme (points i-1,i,i+1,i+2)
		{
			for (int iv = 0; iv < nvar; iv++)
			{
				for (int j = 0; j < NPOINTFD; j++)
				{
					du[IDX(i, iv, nvar)] = du[IDX(i, iv, nvar)] +
						FD_Coeff_CenterR_1st[i*NPOINTFD + j] * uu[IDX(i - 1 + j, iv, nvar)];
				}
			}
		}
		else if (i == 2) // 3rd point -- centerL scheme : but centerR scheme is also okay
		{
			for (int iv = 0; iv < nvar; iv++)
			{
				for (int j = 0; j < NPOINTFD; j++)
				{
					du[IDX(i, iv, nvar)] = du[IDX(i, iv, nvar)] + 
						FD_Coeff_CenterL_1st[i*NPOINTFD + j] * uu[IDX(i - 2 + j,iv,nvar)];

					/*du[IDX(i, iv, nvar)] = du[IDX(i, iv, nvar)] +
						FD_Coeff_CenterR_1st[i*NPOINTFD + j] * uu[IDX(i - 1 + j, iv, nvar)];*/
				}
			}
		}
		else if (i == nx - 1) // right boundary point -- Left scheme
		{
			for (int iv = 0; iv < nvar; iv++)
			{
				for (int j = 0; j < NPOINTFD; j++)
				{
					du[IDX(i, iv, nvar)] = du[IDX(i, iv, nvar)] +
						FD_Coeff_Left_1st[i*NPOINTFD + j] * uu[IDX(i - j, iv, nvar)];
				}
			}
		}
		else // inner points -- need the wave speed to determine schemes
		{
			// copy fields at one point
			for (int iv = 0; iv < nvar; iv++)
			{
				utmp[iv] = uu[IDX(i, iv, nvar)];
			}

			// calculate wave speeds
			waveSpeeds = calcEigens(utmp, Br[i]);

			if (waveSpeeds[0] > 0) // all wave speeds are to the right -- use left scheme
			{
				for (int iv = 0; iv < nvar; iv++)
				{
					for (int j = 0; j < NPOINTFD; j++)
					{
						du[IDX(i, iv, nvar)] = du[IDX(i, iv, nvar)] +
							FD_Coeff_Left_1st[i*NPOINTFD + j] * uu[IDX(i - j, iv, nvar)];
					}
				}
			}
			else if (waveSpeeds[NEigen - 1] < 0) // all wave speeds are to the left -- use right scheme
			{
				if (i == nx - 3) // cannot use right scheme (out of bound), use centerR
				{
					for (int iv = 0; iv < nvar; iv++)
					{
						for (int j = 0; j < NPOINTFD; j++)
						{
							du[IDX(i, iv, nvar)] = du[IDX(i, iv, nvar)] +
								FD_Coeff_CenterR_1st[i*NPOINTFD + j] * uu[IDX(i - 1 + j, iv, nvar)];
						}
					}
				}
				else if (i == nx - 2) // use centerL
				{
					for (int iv = 0; iv < nvar; iv++)
					{
						for (int j = 0; j < NPOINTFD; j++)
						{
							du[IDX(i, iv, nvar)] = du[IDX(i, iv, nvar)] +
								FD_Coeff_CenterL_1st[i*NPOINTFD + j] * uu[IDX(i - 2 + j, iv, nvar)];
						}
					}
				}
				else  // i<=nx-4
				{
					for (int iv = 0; iv < nvar; iv++)
					{
						for (int j = 0; j < NPOINTFD; j++)
						{
							du[IDX(i, iv, nvar)] = du[IDX(i, iv, nvar)] +
								FD_Coeff_Right_1st[i*NPOINTFD + j] * uu[IDX(i + j, iv, nvar)];
						}
					}
				}
			}
			else // waves go to both directions -- use centerL scheme, but centerR is also okay
			{
				for (int iv = 0; iv < nvar; iv++)
				{
					for (int j = 0; j < NPOINTFD; j++)
					{
						du[IDX(i, iv, nvar)] = du[IDX(i, iv, nvar)] +
							FD_Coeff_CenterL_1st[i*NPOINTFD + j] * uu[IDX(i - 2 + j, iv, nvar)];
					}
				}
			}
			free(waveSpeeds);
		}
	}

	free(utmp);

	return du;
}


real* derivXX(real *uu, real *xgrid, real *Br, int nx, int nvar)
{
	real *du;
	du = (real *)malloc(sizeof(real) * nx * nvar);
	memset(du, 0, sizeof(real) * nx * nvar);

	real *utmp;
	utmp = (real *)malloc(sizeof(real) * nvar);
	memset(utmp, 0, sizeof(real) * nvar);

	real *waveSpeeds;

	for (int i = 0; i < nx; i++)
	{
		if (i == 0)// left most point -- right scheme
		{
			for (int iv = 0; iv < nvar; iv++)
			{
				for (int j = 0; j < NPOINTFD; j++)
				{
					du[IDX(i, iv, nvar)] = du[IDX(i, iv, nvar)] +
						FD_Coeff_Right_2nd[i*NPOINTFD + j] * uu[IDX(i + j, iv, nvar)];
				}
			}
		}
		else if (i == 1)// 2nd point -- centerR scheme (points i-1,i,i+1,i+2)
		{
			for (int iv = 0; iv < nvar; iv++)
			{
				for (int j = 0; j < NPOINTFD; j++)
				{
					du[IDX(i, iv, nvar)] = du[IDX(i, iv, nvar)] +
						FD_Coeff_CenterR_2nd[i*NPOINTFD + j] * uu[IDX(i - 1 + j, iv, nvar)];
				}
			}
		}
		else if (i == 2) // 3rd point -- centerL scheme : but centerR scheme is also okay
		{
			for (int iv = 0; iv < nvar; iv++)
			{
				for (int j = 0; j < NPOINTFD; j++)
				{
					du[IDX(i, iv, nvar)] = du[IDX(i, iv, nvar)] +
						FD_Coeff_CenterL_2nd[i*NPOINTFD + j] * uu[IDX(i - 2 + j, iv, nvar)];
				}
			}
		}
		else if (i == nx - 1) // right boundary point -- Left scheme
		{
			for (int iv = 0; iv < nvar; iv++)
			{
				for (int j = 0; j < NPOINTFD; j++)
				{
					du[IDX(i, iv, nvar)] = du[IDX(i, iv, nvar)] +
						FD_Coeff_Left_2nd[i*NPOINTFD + j] * uu[IDX(i - j, iv, nvar)];
				}
			}
		}
		else // inner points -- need the wave speed to determine schemes
		{
			// copy fields at one point
			for (int iv = 0; iv < nvar; iv++)
			{
				utmp[iv] = uu[IDX(i, iv, nvar)];
			}

			// calculate wave speeds
			waveSpeeds = calcEigens(utmp, Br[i]);

			if (waveSpeeds[0] > 0) // all wave speeds are to the right -- use left scheme
			{
				for (int iv = 0; iv < nvar; iv++)
				{
					for (int j = 0; j < NPOINTFD; j++)
					{
						du[IDX(i, iv, nvar)] = du[IDX(i, iv, nvar)] +
							FD_Coeff_Left_2nd[i*NPOINTFD + j] * uu[IDX(i - j, iv, nvar)];
					}
				}
			}
			else if (waveSpeeds[NEigen - 1] < 0) // all wave speeds are to the left -- use right scheme
			{
				if (i == nx - 3) // cannot use right scheme (out of bound), use centerR
				{
					for (int iv = 0; iv < nvar; iv++)
					{
						for (int j = 0; j < NPOINTFD; j++)
						{
							du[IDX(i, iv, nvar)] = du[IDX(i, iv, nvar)] +
								FD_Coeff_CenterR_2nd[i*NPOINTFD + j] * uu[IDX(i - 1 + j, iv, nvar)];
						}
					}
				}
				else if (i == nx - 2) // use centerL
				{
					for (int iv = 0; iv < nvar; iv++)
					{
						for (int j = 0; j < NPOINTFD; j++)
						{
							du[IDX(i, iv, nvar)] = du[IDX(i, iv, nvar)] +
								FD_Coeff_CenterL_2nd[i*NPOINTFD + j] * uu[IDX(i - 2 + j, iv, nvar)];
						}
					}
				}
				else  // i<=nx-4
				{
					for (int iv = 0; iv < nvar; iv++)
					{
						for (int j = 0; j < NPOINTFD; j++)
						{
							du[IDX(i, iv, nvar)] = du[IDX(i, iv, nvar)] +
								FD_Coeff_Right_2nd[i*NPOINTFD + j] * uu[IDX(i + j, iv, nvar)];
						}
					}
				}
			}
			else // waves go to both directions -- use centerL scheme, but centerR is also okay
			{
				for (int iv = 0; iv < nvar; iv++)
				{
					for (int j = 0; j < NPOINTFD; j++)
					{
						du[IDX(i, iv, nvar)] = du[IDX(i, iv, nvar)] +
							FD_Coeff_CenterL_2nd[i*NPOINTFD + j] * uu[IDX(i - 2 + j, iv, nvar)];
					}
				}
			}
			free(waveSpeeds);
		}

	}

	free(utmp);

	return du;
}
#ifndef DJANGO_FDM_1D_H_
#define DJANGO_FDM_1D_H_

#include "FDM.h"

#include "acquisition.h"
#include "grid_1D_float.h"
#include "data.h"
#include "snapshot.h"
#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class FDM_1D: public FDM

{

public:

	FDM_1D(void) ;

	// initialize
	virtual Rtn_code initialize(Grid* pGrid, Myfloat fmax) ;

	// finalize
	virtual Rtn_code finalize(void) ;

	// allocate grid for frequency wavefield extraction
	Rtn_code allocate_grid_for_freq_wavefield(Grid**) ;

	// deallocate grid for frequency wavefield extraction
	Rtn_code deallocate_grid_for_freq_wavefield(Grid*) ;

protected:

	// nb of grid points in medium along axis z
	Myint nz ;

	// grid point spacing along axis z
	Myfloat dz ;

	// inverse of grid spacing
	Myfloat inv_dz ;

	// inverse of grid spacing square
	Myfloat inv_dz2 ;

	// type of surrounding layer
	Boundary_type boundary_type_zBeg, boundary_type_zEnd ;

	// nb of grid points in surrounding layer along axis z
	Myint nlayer_zBeg, nlayer_zEnd ;

	// usefull grid indices
	Myint izBeg, izBeg1, izBeg2, izEnd2, izEnd1, izEnd ;

	// **** acquisition ****
	// number of grid points for the source excitation
	Myint npoint_src ;

	// position of source in the grid along axis z (grid point index)
	Myint* iz_src ;

	// weigth of source in the grid along axis z (grid point index)
	Myfloat* weight_src ;

	// position of receivers in the grid along axis z (grid point index)
	Myint* iz_rec ;

	// position of pixel in the grid along axis z (grid point index)
	Myint* iz_pixel ;

	// compute energy
	void compute_energy(Myfloat*) ;

	// locate src and rec in the grid
	virtual Rtn_code locate_src_and_rec_in_grid(Acquisition*) ;

	// locate snapshot pixel in the grid
	virtual Rtn_code locate_pixel_in_grid(Snapshot*) ;

	// write time solution at receivers
	Rtn_code write_trace(Variable*, Myint) ;

	// free src and rec in the grid
	Rtn_code free_position_arrays(void) ;

	// write time solution snapshot
	Rtn_code write_time_sol_snapshot(Myfloat*, const char*, Myint) ;

	// write freq solution at receivers
	Rtn_code write_freq_sol_rec(Mycomplex**) ;

	// store freq solution at receivers
	Rtn_code store_freq_sol_rec(Mycomplex**, Data*) ;

	// write freq solution at all grid points
	Rtn_code write_freq_sol_grid(Mycomplex**) ;

	// write time solution at all grid points
	Rtn_code write_time_sol_grid(Myfloat*, ofstream*) ;

	// extract frequency solution from the time solution of the FDM grid
	Rtn_code update_freq_sol_grid(Myfloat*, Mycomplex**, Myint) ;

	// source excitation
	Rtn_code source_excitation(Myfloat*, Myint, Wavefield_type, Data*) ;
	Rtn_code source_excitation(Myfloat*, Myint, Wavefield_type, Data*, Myfloat) ;

	// check time step
	Myfloat compute_optimal_time_step(Grid_1D_float *vp) ;

} ;

} // namespace django

#endif

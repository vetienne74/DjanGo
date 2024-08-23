#ifndef DJANGO_FDM_2D_H_
#define DJANGO_FDM_2D_H_

#include "FDM_1D.h"

#include "acquisition.h"
#include "data.h"
#include "grid.h"
#include "grid_2D_float.h"
#include "snapshot.h"
#include "type_def.h"

namespace django {

//------------------------------------------------------------------------------------

class FDM_2D: public FDM_1D

{

public:

	FDM_2D(void) ;

	// initialize
	virtual Rtn_code initialize(Grid* pGrid, Myfloat fmax) ;

	// finalize
	virtual Rtn_code finalize(void) ;

	// allocate grid for frequency wavefield extraction
	Rtn_code allocate_grid_for_freq_wavefield(Grid**) ;

	// deallocate grid for frequency wavefield extraction
	Rtn_code deallocate_grid_for_freq_wavefield(Grid*) ;

	// Write snapshot
	virtual Rtn_code write_snapshot(Variable*, Myint) ;

protected:

	// nb of grid points along axis x
	Myint nx ;

	// grid point spacing along axis x
	Myfloat dx ;

	// inverse of grid spacing
	Myfloat inv_dx ;

	// inverse of grid spacing square
	Myfloat inv_dx2 ;

	// nb of grid points in surrounding layer along axis x
	Myint nlayer_xBeg, nlayer_xEnd ;

	// usefull grid indices
	Myint ixBeg, ixBeg1, ixBeg2, ixEnd2, ixEnd1, ixEnd ;

	// **** acquisition ****
	// position of source in the grid along axis x
	Myint* ix_src ;

	// position of receivers in the grid along axis x
	Myint* ix_rec ;

	// position of pixel in the grid along axis x
	Myint* ix_pixel ;

	// compute energy
	void compute_energy(Myfloat**) ;

	// compute error for eigen mode and write on disk
	Rtn_code compute_eigen_error(Myfloat**, Myint, Acquisition*) ;

	// locate src and rec in the grid
	virtual Rtn_code locate_src_and_rec_in_grid(Acquisition*) ;

	// locate snapshot pixel in the grid
	virtual Rtn_code locate_pixel_in_grid(Snapshot*) ;

	// free src and rec in the grid
	Rtn_code free_position_arrays(void) ;

	// write time solution at receivers
	Rtn_code write_trace(Variable*, Myint) ;

	// write time solution snapshot
	Rtn_code write_time_sol_snapshot(Myfloat**, const char*, Myint) ;

	// write freq solution at receivers
	Rtn_code write_freq_sol_rec(Data*) ;

	// write freq solution at all grid points
	Rtn_code write_freq_sol_grid(Mycomplex***) ;

	// write time solution at all grid points
	Rtn_code write_time_sol_grid(Myfloat**, ofstream*) ;

	// extract frequency solution at all grid points
	Rtn_code update_freq_sol_grid(Myfloat**, Mycomplex***, Myint) ;

	// extract frequency solution at receivers
	Rtn_code update_freq_sol_rec(Myfloat**, Data*, Myint) ;

	// source excitation
	Rtn_code source_excitation(Myfloat**, Myint, Wavefield_type, Data*) ;
	Rtn_code source_excitation(Myfloat**, Myint, Wavefield_type, Data*, Myfloat) ;

	// check time step
	Myfloat compute_optimal_time_step(Grid_2D_float *vp) ;

} ;

} // namespace django

#endif

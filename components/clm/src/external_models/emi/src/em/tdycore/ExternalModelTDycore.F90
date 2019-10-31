module ExternalModelTDycore

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscdm.h"
#include "finclude/tdycore.h"

  use abortutils                   , only : endrun
  use shr_kind_mod                 , only : r8 => shr_kind_r8
  use shr_log_mod                  , only : errMsg => shr_log_errMsg
  use EMI_DataMod                  , only : emi_data_list, emi_data
  use clm_varctl                   , only : iulog
  use ExternalModelBaseType        , only : em_base_type
  use decompMod                    , only : bounds_type
  use ExternalModelConstants
  use EMI_Atm2LndType_Constants
  use EMI_CanopyStateType_Constants
  use EMI_ColumnType_Constants
  use EMI_EnergyFluxType_Constants
  use EMI_Filter_Constants
  use EMI_Landunit_Constants
  use EMI_SoilHydrologyType_Constants
  use EMI_SoilStateType_Constants
  use EMI_TemperatureType_Constants
  use EMI_WaterFluxType_Constants
  use EMI_WaterStateType_Constants
  use elm_tdycore_interface_data

  !use petscvec
  use petscdm
  use petscts
  use Mapping_module

#ifdef USE_TDYCORE
  use tdycore
#endif

  implicit none

  type, public, extends(em_base_type) :: em_tdycore_type
     ! ----------------------------------------------------------------------
     ! Indicies required tdyduring the initialization
     ! ----------------------------------------------------------------------
     integer :: index_l2e_init_col_active
     integer :: index_l2e_init_col_type
     integer :: index_l2e_init_col_landunit_index
     integer :: index_l2e_init_col_gridcell_index
     integer :: index_l2e_init_col_zi
     integer :: index_l2e_init_col_dz
     integer :: index_l2e_init_col_z
     integer :: index_l2e_init_col_area

     integer :: index_l2e_init_state_wtd
     integer :: index_l2e_init_state_soilp

     integer :: index_l2e_init_h2osoi_liq
     integer :: index_l2e_init_h2osoi_ice

     integer :: index_e2l_init_state_h2osoi_liq
     integer :: index_e2l_init_state_h2osoi_ice
     integer :: index_e2l_init_state_h2osoi_vol
     integer :: index_e2l_init_state_wtd
     integer :: index_e2l_init_parameter_watsatc
     integer :: index_e2l_init_parameter_hksatc
     integer :: index_e2l_init_parameter_bswc
     integer :: index_e2l_init_parameter_sucsatc

     integer :: index_e2l_init_flux_mflx_snowlyr_col
     integer :: index_l2e_init_flux_mflx_snowlyr_col

     integer :: index_l2e_init_landunit_type
     integer :: index_l2e_init_landunit_lakepoint
     integer :: index_l2e_init_landunit_urbanpoint

     integer :: index_l2e_init_parameter_watsatc
     integer :: index_l2e_init_parameter_hksatc
     integer :: index_l2e_init_parameter_bswc
     integer :: index_l2e_init_parameter_sucsatc
     integer :: index_l2e_init_parameter_effporosityc

     ! ----------------------------------------------------------------------
     ! Indicies required during timestepping
     ! ----------------------------------------------------------------------

     integer :: index_l2e_state_tsoil
     integer :: index_l2e_state_h2osoi_liq
     integer :: index_l2e_state_h2osoi_ice

     integer :: index_e2l_state_h2osoi_liq
     integer :: index_e2l_state_h2osoi_ice
     integer :: index_e2l_state_wtd
     integer :: index_e2l_state_soilp

     integer :: index_l2e_flux_infil
     integer :: index_l2e_flux_et
     integer :: index_l2e_flux_dew
     integer :: index_l2e_flux_snow_sub
     integer :: index_l2e_flux_snowlyr
     integer :: index_l2e_flux_drainage

     integer :: index_e2l_flux_qrecharge
     integer :: index_e2l_flux_drain_perched
     integer :: index_e2l_flux_drain
     integer :: index_e2l_flux_qrgwl
     integer :: index_e2l_flux_rsub_sat

     integer :: index_l2e_filter_hydrologyc
     integer :: index_l2e_filter_num_hydrologyc

     integer :: index_l2e_column_zi
     integer :: index_l2e_col_active
     integer :: index_l2e_col_gridcell
     integer :: index_l2e_col_dz

#ifdef USE_TDYCORE
     TDy, pointer :: tdy
     DM :: dm
     TS :: ts
#endif
    type(mapping_type), pointer :: map_elm_sub_to_tdy_sub
    type(mapping_type), pointer :: map_tdy_sub_to_elm_sub

   contains

     procedure, public :: Populate_L2E_Init_List  => EM_TDycore_Populate_L2E_Init_List
     procedure, public :: Populate_E2L_Init_List  => EM_TDycore_Populate_E2L_Init_List
     procedure, public :: Populate_L2E_List       => EM_TDycore_Populate_L2E_List
     procedure, public :: Populate_E2L_List       => EM_TDycore_Populate_E2L_List
     procedure, public :: PreInit                 => EM_TDycore_PreInit
     procedure, public :: Init                    => EM_TDycore_Init
     procedure, public :: Solve                   => EM_TDycore_Solve

  end type em_tdycore_type

contains

  !------------------------------------------------------------------------
  subroutine EM_TDycore_Populate_L2E_Init_List(this, l2e_init_list)
    !
    implicit none
    !
    class(em_tdycore_type) :: this
    class(emi_data_list), intent(inout) :: l2e_init_list

    !
    ! !LOCAL VARIABLES:
    integer              , pointer       :: em_stages(:)
    integer                              :: number_em_stages
    integer                              :: id
    integer                              :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_INITIALIZATION_STAGE

    id                                         = L2E_STATE_WTD
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_state_wtd              = index

    id                                         = L2E_STATE_VSFM_PROGNOSTIC_SOILP
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_state_soilp            = index

    id                                         = L2E_FLUX_RESTART_SNOW_LYR_DISAPPERANCE_MASS_FLUX
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_flux_mflx_snowlyr_col  = index

    id                                         = L2E_COLUMN_ACTIVE
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_active             = index

    id                                         = L2E_COLUMN_TYPE
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_type               = index

    id                                         = L2E_COLUMN_LANDUNIT_INDEX
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_landunit_index     = index

    id                                         = L2E_COLUMN_GRIDCELL_INDEX
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_gridcell_index     = index

    id                                         = L2E_COLUMN_ZI
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_zi                 = index

    id                                         = L2E_COLUMN_DZ
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_dz                 = index

    id                                         = L2E_COLUMN_Z
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_z                  = index

    id                                         = L2E_COLUMN_AREA
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_col_area               = index

    id                                         = L2E_LANDUNIT_TYPE
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_landunit_type          = index

    id                                         = L2E_LANDUNIT_LAKEPOINT
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_landunit_lakepoint     = index

    id                                         = L2E_LANDUNIT_URBANPOINT
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_landunit_urbanpoint    = index

    id                                         = L2E_PARAMETER_WATSATC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_watsatc      = index

    id                                         = L2E_PARAMETER_HKSATC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_hksatc       = index

    id                                         = L2E_PARAMETER_BSWC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_bswc         = index

    id                                         = L2E_PARAMETER_SUCSATC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_sucsatc      = index

    id                                         = L2E_PARAMETER_EFFPOROSITYC
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_parameter_effporosityc = index

    id                                        = L2E_STATE_H2OSOI_LIQ_NLEVGRND
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_h2osoi_liq            = index

    id                                        = L2E_STATE_H2OSOI_ICE_NLEVGRND
    call l2e_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_init_h2osoi_ice            = index

    deallocate(em_stages)

  end subroutine EM_TDycore_Populate_L2E_Init_List

  !------------------------------------------------------------------------
  subroutine EM_TDycore_Populate_E2L_Init_List(this, e2l_init_list)
    !
    implicit none
    !
    class(em_tdycore_type) :: this
    class(emi_data_list), intent(inout) :: e2l_init_list
    !
    ! !LOCAL VARIABLES:
    integer              , pointer       :: em_stages(:)
    integer                              :: number_em_stages
    integer                              :: id
    integer                              :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_INITIALIZATION_STAGE

    id                                        = E2L_STATE_H2OSOI_LIQ
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_h2osoi_liq      = index

    id                                        = E2L_STATE_H2OSOI_ICE
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_h2osoi_ice      = index

    id                                        = E2L_STATE_H2OSOI_VOL_NLEVGRND
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_h2osoi_vol      = index

    id                                        = E2L_STATE_WTD
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_state_wtd             = index

    id                                        = E2L_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_flux_mflx_snowlyr_col = index

    id                                         = E2L_PARAMETER_WATSATC
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_parameter_watsatc      = index

    id                                         = E2L_PARAMETER_HKSATC
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_parameter_hksatc       = index

    id                                         = E2L_PARAMETER_BSWC
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_parameter_bswc         = index

    id                                         = E2L_PARAMETER_SUCSATC
    call e2l_init_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_init_parameter_sucsatc      = index

    deallocate(em_stages)

  end subroutine EM_TDycore_Populate_E2L_Init_List

  !------------------------------------------------------------------------
  subroutine EM_TDycore_Populate_L2E_List(this, l2e_list)
    !
    implicit none
    !
    class(em_tdycore_type) :: this
    class(emi_data_list), intent(inout) :: l2e_list
    !
    !
    ! !LOCAL VARIABLES:
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_TDYCORE_SOIL_HYDRO_STAGE

    id                                   = L2E_STATE_TSOIL_NLEVGRND
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_tsoil           = index

    id                                   = L2E_STATE_H2OSOI_LIQ_NLEVGRND
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osoi_liq      = index

    id                                   = L2E_STATE_H2OSOI_ICE_NLEVGRND
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osoi_ice      = index

    id                                   = L2E_FLUX_INFIL_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_infil            = index

    id                                   = L2E_FLUX_VERTICAL_ET_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_et               = index

    id                                   = L2E_FLUX_DEW_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_dew              = index

    id                                   = L2E_FLUX_SNOW_SUBLIMATION_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_snow_sub         = index

    id                                   = L2E_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_snowlyr          = index

    id                                   = L2E_FLUX_DRAINAGE_MASS_FLUX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_drainage         = index

    id                                   = L2E_FILTER_HYDROLOGYC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_filter_hydrologyc     = index

    id                                   = L2E_FILTER_NUM_HYDROLOGYC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_filter_num_hydrologyc = index

    id                                   = L2E_COLUMN_ZI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_column_zi             = index

    id                                   = L2E_COLUMN_ACTIVE
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_col_active            = index

    id                                   = L2E_COLUMN_GRIDCELL_INDEX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_col_gridcell          = index

    id                                    = L2E_COLUMN_DZ
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_col_dz                 = index

    deallocate(em_stages)


  end subroutine EM_TDycore_Populate_L2E_List

  !------------------------------------------------------------------------
  subroutine EM_TDycore_Populate_E2L_List(this, e2l_list)
    !
    implicit none
    !
    class(em_tdycore_type) :: this
    class(emi_data_list), intent(inout) :: e2l_list
    !
    ! !LOCAL VARIABLES:
    integer              , pointer       :: em_stages(:)
    integer                              :: number_em_stages
    integer                              :: id
    integer                              :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))
    em_stages(1) = EM_TDYCORE_SOIL_HYDRO_STAGE

    id                              = E2L_STATE_H2OSOI_LIQ
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_h2osoi_liq = index

    id                              = E2L_STATE_H2OSOI_ICE
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_h2osoi_ice = index

    id                              = E2L_STATE_WTD
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_wtd        = index

    id                              = E2L_STATE_VSFM_PROGNOSTIC_SOILP
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_state_soilp      = index

    id                              = E2L_FLUX_AQUIFER_RECHARGE
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_qrecharge   = index

    id                              = E2L_FLUX_DRAIN_PERCHED
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_drain_perched   = index

    id                              = E2L_FLUX_DRAIN
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_drain       = index

    id                              = E2L_FLUX_QRGWL
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_qrgwl   = index

    id                              = E2L_FLUX_RSUB_SAT
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_rsub_sat   = index

    deallocate(em_stages)


  end subroutine EM_TDycore_Populate_E2L_List

  !------------------------------------------------------------------------
  subroutine EM_TDycore_PreInit(this,bounds)
    !
    use spmdMod     , only : mpicom
    use clm_varpar  , only : nlevgrnd
    !
    implicit none
    !
    class(em_tdycore_type) :: this
    character(len=4096)    :: mesh_filename
    type(bounds_type), intent(in)  :: bounds
    PetscErrorCode         :: ierr
    DM                     :: dmDist
    PetscBool              :: flg
    PetscInt               :: dim, faces(3)
    PetscReal              :: lower(3), upper(3)

    allocate(this%tdy)

    call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
         "-mesh_filename",mesh_filename,flg,ierr);CHKERRQ(ierr)

    if (.not.flg) then
       !call endrun(msg='ERROR EM_TDycore_PreInit: '//&
       !   'Need to specify a mesh file for TDycore via runtime option -tdy_mesh_filename')
       faces(1) = bounds%endg-bounds%begg+1;
       faces(2) = 1;
       faces(3) = nlevgrnd;
       lower(:) = 0.d0;
       upper(:) = 1.d0;
       dim = 3

       call DMPlexCreateBoxMesh(PETSC_COMM_WORLD, dim, PETSC_FALSE, faces, lower, upper, &
            PETSC_NULL_INTEGER, PETSC_TRUE, this%dm, ierr);
       CHKERRA(ierr);
    else
       call DMPlexCreateFromFile(PETSC_COMM_WORLD, mesh_filename, PETSC_TRUE, this%dm, ierr); CHKERRQ(ierr);
    endif
    
    call DMSetFromOptions(this%dm, ierr); CHKERRA(ierr);
    call DMPlexDistribute(this%dm, 1, PETSC_NULL_SF, dmDist, ierr);CHKERRA(ierr);
    if (dmDist /= PETSC_NULL_DM) then
       call DMDestroy(this%dm, ierr);CHKERRA(ierr);
       this%dm = dmDist;
    end if
    call DMView(this%dm,PETSC_VIEWER_STDOUT_WORLD,ierr); CHKERRA(ierr)

    call TDyCreate(this%dm, this%tdy, ierr);
    CHKERRA(ierr);

  end subroutine EM_TDycore_PreInit

  !------------------------------------------------------------------------
  subroutine EM_TDycore_Init(this, l2e_init_list, e2l_init_list, iam, bounds_clump)
    !
    implicit none
    !
    class(em_tdycore_type)               :: this
    class(emi_data_list) , intent(in)    :: l2e_init_list
    class(emi_data_list) , intent(inout) :: e2l_init_list
    integer              , intent(in)    :: iam
    type(bounds_type)    , intent(in)    :: bounds_clump
    !
    ! LOCAL VARAIBLES:
    integer :: clm_npts
    integer :: clm_surf_npts
    PetscInt g, j, count, ngridcells
    Vec :: U
    PetscReal, pointer :: p_loc(:)
    PetscReal, pointer :: hksat_x(:), hksat_y(:), hksat_z(:), watsat(:), thetares(:)
    PetscErrorCode :: ierr

    ! Create ELM-TDycore interface data
    call CreateELMTDycoreInterfaceData(this, bounds_clump)
    
    ! ==========================================================================
    ! Create maps
    nullify(this%map_elm_sub_to_tdy_sub)
    nullify(this%map_tdy_sub_to_elm_sub)

    this%map_elm_sub_to_tdy_sub => MappingCreate()
    this%map_tdy_sub_to_elm_sub => MappingCreate()

    call InitializeMap(this, this%map_elm_sub_to_tdy_sub)
    call InitializeMap(this, this%map_tdy_sub_to_elm_sub)
    ! ==========================================================================

    ! ==========================================================================
    ! Set material properties
    call InitializationData_ELM2TDycore(this)
    ! ==========================================================================

    ! ==========================================================================
    ! Finish setup of TDycore
    call TDySetDiscretizationMethod(this%tdy,MPFA_O,ierr);
    CHKERRA(ierr);

    call TDySetFromOptions(this%tdy,ierr);
    CHKERRA(ierr);
    ! ==========================================================================

    ! ==========================================================================
    ! Set initial condition
    call DMCreateGlobalVector(this%dm,U,ierr);
    CHKERRA(ierr);

    call VecSet(U,91325.d0,ierr);
    CHKERRA(ierr);

    call VecGetArrayF90(U,p_loc,ierr);
    CHKERRA(ierr);

    call TDyUpdateState(this%tdy,p_loc,ierr);
    CHKERRA(ierr);

    call VecRestoreArrayF90(U,p_loc,ierr);
    CHKERRA(ierr);
    ! ==========================================================================

    ! ==========================================================================
    ! Set material properties
    call InitializationData_TDycore2ELM(this, l2e_init_list, e2l_init_list, bounds_clump)
    ! ==========================================================================

    ! Save the data need by ELM
    !call extract_data_for_elm(this, l2e_init_list, e2l_init_list, bounds_clump)

  end subroutine EM_TDycore_Init

  !-----------------------------------------------------------------------

  subroutine InitializationData_TDycore2ELM(this,l2e_init_list, e2l_init_list, bounds_clump)

    implicit none

    class(em_tdycore_type)  , intent(inout) :: this
    class(emi_data_list) , intent(in)    :: l2e_init_list
    class(emi_data_list) , intent(inout) :: e2l_init_list
    type(bounds_type)    , intent(in)    :: bounds_clump

    call PackInitializationData_TDycore2ELM(this)
    call MapInitializationData_TDycore2ELM(this)
    call UnpackInitializationData_TDycore2ELM(this,l2e_init_list, e2l_init_list, bounds_clump)

  end subroutine InitializationData_TDycore2ELM

  !-----------------------------------------------------------------------
  subroutine PackInitializationData_TDycore2ELM(this)

    implicit none
    class(em_tdycore_type)  , intent(inout) :: this

    PetscInt :: nvalues, ngridcells, g
    PetscReal, pointer :: liquid_mass(:), liquid_sat(:), matprop_m(:), matprop_alpha(:), porosity(:), blockPerm(:)
    PetscReal, pointer :: mass(:), sat(:), sucsat(:), bsw(:), watsat(:), hksat(:)
    PetscReal :: den, vis, grav
    PetscErrorCode :: ierr
    character(len= 128)   :: subname = 'PackInitializationData_TDycore2ELM'

    call VecGetSize(elm_tdycore_idata%hksat_x_tdycore_svec,ngridcells,ierr);
    CHKERRA(ierr);

    allocate(liquid_sat(ngridcells))
    allocate(liquid_mass(ngridcells))
    allocate(matprop_m(ngridcells))
    allocate(matprop_alpha(ngridcells))
    allocate(porosity(ngridcells))
    allocate(blockPerm(ngridcells*9))

    call TDyGetBlockPermeabilityValuesLocal(this%tdy,nvalues,blockPerm,ierr)
    CHKERRA(ierr);

    call TDyGetSaturationValuesLocal(this%tdy,nvalues,liquid_sat,ierr)
    CHKERRA(ierr);

    call TDyGetLiquidMassValuesLocal(this%tdy,nvalues,liquid_mass,ierr)
    CHKERRA(ierr);

    call TDyGetMaterialPropertyMValuesLocal(this%tdy,nvalues,matprop_m,ierr)
    CHKERRA(ierr);

    call TDyGetMaterialPropertyAlphaValuesLocal(this%tdy,nvalues,matprop_alpha,ierr)
    CHKERRA(ierr);

    call TDyGetPorosityValuesLocal(this%tdy,nvalues,porosity,ierr)
    CHKERRA(ierr);

    if (ngridcells /= nvalues) then
      write(iulog,*)'Number of local grid cells in TDycore does not match the grid cells in the map'
      call endrun( trim(subname)//' ERROR: Land Unit type not supported' )
    endif

    call VecGetArrayF90(elm_tdycore_idata%mass_tdycore_mvec, mass, ierr)
    CHKERRA(ierr);
    call VecGetArrayF90(elm_tdycore_idata%sat_tdycore_mvec, sat, ierr)
    CHKERRA(ierr);
    call VecGetArrayF90(elm_tdycore_idata%sucsat_tdycore_mvec, sucsat, ierr)
    CHKERRA(ierr);
    call VecGetArrayF90(elm_tdycore_idata%bsw_tdycore_mvec, bsw, ierr)
    CHKERRA(ierr);
    call VecGetArrayF90(elm_tdycore_idata%watsat_tdycore_mvec, watsat, ierr)
    CHKERRA(ierr);
    call VecGetArrayF90(elm_tdycore_idata%hksat_x_tdycore_mvec, hksat, ierr)
    CHKERRA(ierr);

    den = 998.2d0       ! [kg/m^3]  @ 20 degC
    vis = 0.001002d0    ! [N s/m^2] @ 20 degC
    grav = 9.81d0       ! [m/s^2]

    do g = 1, ngridcells
      mass(g) = liquid_mass(g)
      sat(g) = liquid_sat(g)
      bsw(g) = 1.0/matprop_m(g)
      sucsat(g) = 1.0/(matprop_alpha(g) * grav)
      watsat(g) = porosity(g)
      hksat(g) = blockPerm((g-1)*9+1)/vis*(den*grav)*1000.d0
    enddo

    call VecRestoreArrayF90(elm_tdycore_idata%mass_tdycore_mvec, mass, ierr)
    CHKERRA(ierr);
    call VecRestoreArrayF90(elm_tdycore_idata%sat_tdycore_mvec, sat, ierr)
    CHKERRA(ierr);
    call VecRestoreArrayF90(elm_tdycore_idata%sucsat_tdycore_mvec, sucsat, ierr)
    CHKERRA(ierr);
    call VecRestoreArrayF90(elm_tdycore_idata%bsw_tdycore_mvec, bsw, ierr)
    CHKERRA(ierr);
    call VecRestoreArrayF90(elm_tdycore_idata%watsat_tdycore_mvec, watsat, ierr)
    CHKERRA(ierr);
    call VecRestoreArrayF90(elm_tdycore_idata%hksat_x_tdycore_mvec, hksat, ierr)
    CHKERRA(ierr);

    deallocate(liquid_sat)
    deallocate(liquid_mass)
    deallocate(matprop_m)
    deallocate(matprop_alpha)
    deallocate(porosity)

  end subroutine PackInitializationData_TDycore2ELM

  !-----------------------------------------------------------------------
  subroutine MapInitializationData_TDycore2ELM(this)

    implicit none
    class(em_tdycore_type)  , intent(inout) :: this

    call MappingSourceToDestination(this%map_tdy_sub_to_elm_sub, &
        elm_tdycore_idata%mass_tdycore_mvec, &
        elm_tdycore_idata%mass_elm_svec)

    call MappingSourceToDestination(this%map_tdy_sub_to_elm_sub, &
        elm_tdycore_idata%sat_tdycore_mvec, &
        elm_tdycore_idata%sat_elm_svec)

    call MappingSourceToDestination(this%map_tdy_sub_to_elm_sub, &
        elm_tdycore_idata%sucsat_tdycore_mvec, &
        elm_tdycore_idata%sucsat_elm_svec)

    call MappingSourceToDestination(this%map_tdy_sub_to_elm_sub, &
        elm_tdycore_idata%bsw_tdycore_mvec, &
        elm_tdycore_idata%bsw_elm_svec)

    call MappingSourceToDestination(this%map_tdy_sub_to_elm_sub, &
        elm_tdycore_idata%watsat_tdycore_mvec, &
        elm_tdycore_idata%watsat_elm_svec)

    call MappingSourceToDestination(this%map_tdy_sub_to_elm_sub, &
        elm_tdycore_idata%hksat_x_tdycore_mvec, &
        elm_tdycore_idata%hksat_x_elm_svec)

  end subroutine MapInitializationData_TDycore2ELM

  !-----------------------------------------------------------------------
  subroutine UnpackInitializationData_TDycore2ELM(this,l2e_init_list, e2l_init_list, bounds_clump)

    use landunit_varcon           , only : istsoil
    use clm_varcon                , only : denice, denh2o
    use clm_varpar                , only : nlevgrnd

    implicit none
    class(em_tdycore_type)  , intent(inout) :: this
    class(emi_data_list) , intent(in)    :: l2e_init_list
    class(emi_data_list) , intent(inout) :: e2l_init_list
    type(bounds_type)    , intent(in)    :: bounds_clump

    ! !LOCAL VARIABLES:
    integer               :: c,j,g,l,pf_j  ! do loop indices
    integer               :: nlevmapped
    integer               :: gcount
    integer               :: bounds_proc_begc, bounds_proc_endc

    integer     , pointer :: col_active(:)
    integer     , pointer :: col_landunit(:)
    integer     , pointer :: col_gridcell(:)
    integer     , pointer :: lun_type(:)

    real(r8)    , pointer :: l2e_h2osoi_ice(:,:)

    real(r8)    , pointer :: e2l_h2osoi_liq(:,:)
    real(r8)    , pointer :: e2l_h2osoi_ice(:,:)
    real(r8)    , pointer :: e2l_h2osoi_vol(:,:)
    real(r8)    , pointer :: e2l_zwt(:)
    real(r8)    , pointer :: e2l_mflx_snowlyr_col(:)
    real(r8)    , pointer :: e2l_watsatc(:,:)
    real(r8)    , pointer :: e2l_hksatc(:,:)
    real(r8)    , pointer :: e2l_bswc(:,:)
    real(r8)    , pointer :: e2l_sucsatc(:,:)

    real(r8)    , pointer :: dz(:,:)

    PetscScalar , pointer :: sat(:)
    PetscScalar , pointer :: watsat(:)
    PetscScalar , pointer :: hksat(:)
    PetscScalar , pointer :: bsw(:)
    PetscScalar , pointer :: sucsat(:)
    PetscErrorCode        :: ierr
    character(len= 128)   :: subname = 'UnpackInitializationData_TDycore2ELM'

    bounds_proc_begc = bounds_clump%begc
    bounds_proc_endc = bounds_clump%endc
    nlevmapped       = elm_tdycore_idata%nzelm_mapped

    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_active             , col_active   )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_landunit_index     , col_landunit )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_gridcell_index     , col_gridcell )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_col_dz                , dz           )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_landunit_type          , lun_type     )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_h2osoi_ice            , l2e_h2osoi_ice       )

    call e2l_init_list%GetPointerToReal1D(this%index_e2l_init_state_wtd             , e2l_zwt              )
    call e2l_init_list%GetPointerToReal1D(this%index_e2l_init_flux_mflx_snowlyr_col , e2l_mflx_snowlyr_col )

    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_state_h2osoi_liq      , e2l_h2osoi_liq       )
    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_state_h2osoi_ice      , e2l_h2osoi_ice       )
    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_state_h2osoi_vol      , e2l_h2osoi_vol       )

    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_parameter_watsatc     , e2l_watsatc )
    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_parameter_hksatc      , e2l_hksatc  )
    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_parameter_bswc        , e2l_bswc    )
    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_parameter_sucsatc     , e2l_sucsatc )

    ! Set initial value of for ELM
    e2l_mflx_snowlyr_col(:) = 0._r8
    e2l_zwt(:)              = 0._r8

    ! Initialize soil moisture
    call VecGetArrayF90(elm_tdycore_idata%sat_elm_svec,sat,ierr)
    CHKERRA(ierr);
    call VecGetArrayF90(elm_tdycore_idata%watsat_elm_svec, watsat , ierr)
    CHKERRA(ierr);
    call VecGetArrayF90(elm_tdycore_idata%sucsat_elm_svec, sucsat , ierr)
    CHKERRA(ierr);
    call VecGetArrayF90(elm_tdycore_idata%bsw_elm_svec, bsw , ierr)
    CHKERRA(ierr);
    call VecGetArrayF90(elm_tdycore_idata%watsat_elm_svec, bsw , ierr)
    CHKERRA(ierr);
    call VecGetArrayF90(elm_tdycore_idata%hksat_x_elm_svec, hksat , ierr)
    CHKERRA(ierr);

    watsat(:) = 0.5d0

    do c = bounds_proc_begc, bounds_proc_endc
       if (col_active(c) == 1) then
          l = col_landunit(c)
          if (lun_type(l) == istsoil) then
             g = col_gridcell(c)
             gcount = g - bounds_clump%begg
             do j = 1, nlevgrnd
                pf_j = gcount*nlevmapped + j

                if (j <= nlevmapped) then
                   e2l_h2osoi_liq(c,j) = sat(pf_j)*watsat(pf_j)*dz(c,j)*1.e3_r8

                   e2l_h2osoi_vol(c,j) = e2l_h2osoi_liq(c,j)/dz(c,j)/denh2o + &
                        l2e_h2osoi_ice(c,j)/dz(c,j)/denice
                   e2l_h2osoi_vol(c,j) = min(e2l_h2osoi_vol(c,j),watsat(pf_j))
                   e2l_h2osoi_ice(c,j) = 0._r8

                   e2l_watsatc(c,j) = watsat(pf_j)
                   e2l_hksatc(c,j)  = hksat(pf_j)
                   e2l_bswc(c,j)    = bsw(pf_j)
                   e2l_sucsatc(c,j) = sucsat(pf_j)
                else
                   e2l_h2osoi_liq(c,j) = e2l_h2osoi_liq(c,nlevmapped)
                   e2l_h2osoi_vol(c,j) = e2l_h2osoi_vol(c,nlevmapped)
                   e2l_h2osoi_ice(c,j) = 0._r8
                   e2l_watsatc(c,j)    = e2l_watsatc(c,nlevmapped)
                   e2l_hksatc(c,j)     = e2l_hksatc(c,nlevmapped)
                   e2l_bswc(c,j)       = e2l_bswc(c,nlevmapped)
                   e2l_sucsatc(c,j)    = e2l_sucsatc(c,nlevmapped)
                end if

             enddo
          else
             write(iulog,*)'WARNING: Land Unit type other than soil type is present within the domain'
             call endrun( trim(subname)//' ERROR: Land Unit type not supported' )
          endif
       endif
    enddo

    call VecRestoreArrayF90(elm_tdycore_idata%sat_elm_svec,sat,ierr)
    CHKERRA(ierr);
    call VecRestoreArrayF90(elm_tdycore_idata%watsat_elm_svec, watsat , ierr)
    CHKERRA(ierr);
    call VecRestoreArrayF90(elm_tdycore_idata%sucsat_elm_svec, sucsat , ierr)
    CHKERRA(ierr);
    call VecRestoreArrayF90(elm_tdycore_idata%bsw_elm_svec, bsw , ierr)
    CHKERRA(ierr);
    call VecRestoreArrayF90(elm_tdycore_idata%watsat_elm_svec, bsw , ierr)
    CHKERRA(ierr);
    call VecRestoreArrayF90(elm_tdycore_idata%hksat_x_elm_svec, hksat , ierr)
    CHKERRA(ierr);

  end subroutine UnpackInitializationData_TDycore2ELM

  !-----------------------------------------------------------------------

  subroutine InitializationData_ELM2TDycore(this)

    implicit none

    class(em_tdycore_type)  , intent(inout) :: this

    call PackMaterialProperties_ELM2TDycore(this)
    call MapMaterialProperties_ELM2TDycore(this)
    call UnpackMaterialProperties_ELM2TDycore(this)

  end subroutine InitializationData_ELM2TDycore

  !-----------------------------------------------------------------------
  subroutine PackMaterialProperties_ELM2TDycore(this)

    implicit none
    class(em_tdycore_type)  , intent(inout) :: this

    PetscReal, pointer :: hksat_x(:), hksat_y(:), hksat_z(:), watsat(:), thetares(:)
    PetscReal :: den, vis, grav, perm, porosity, Sr
    PetscInt :: count, g, j
    PetscErrorCode :: ierr

    call VecGetArrayF90(elm_tdycore_idata%hksat_x_elm_mvec, hksat_x, ierr)
    call VecGetArrayF90(elm_tdycore_idata%hksat_y_elm_mvec, hksat_y, ierr)
    call VecGetArrayF90(elm_tdycore_idata%hksat_z_elm_mvec, hksat_z, ierr)
    call VecGetArrayF90(elm_tdycore_idata%watsat_elm_mvec, watsat, ierr)
    call VecGetArrayF90(elm_tdycore_idata%thetares_elm_mvec, thetares, ierr)

    den = 998.2d0       ! [kg/m^3]  @ 20 degC
    vis = 0.001002d0    ! [N s/m^2] @ 20 degC
    grav = 9.81d0       ! [m/s^2]
    perm = 1.d-10       ! [m^2]
    porosity = 0.5d0
    Sr   = 0.115d0

    count = 0
    do g = 1, elm_tdycore_idata%nlelm_srf
       do j = 1, elm_tdycore_idata%nzelm_mapped
          count = count + 1;
          hksat_x(count) = perm/vis*(den*grav)*1000.d0
          hksat_y(count) = perm/vis*(den*grav)*1000.d0
          hksat_z(count) = perm/vis*(den*grav)*1000.d0

          watsat(count)   = porosity
          thetares(count) = porosity*Sr
       enddo
    enddo

    call VecRestoreArrayF90(elm_tdycore_idata%hksat_x_elm_mvec, hksat_x, ierr)
    CHKERRA(ierr);
    call VecRestoreArrayF90(elm_tdycore_idata%hksat_y_elm_mvec, hksat_y, ierr)
    CHKERRA(ierr);
    call VecRestoreArrayF90(elm_tdycore_idata%hksat_z_elm_mvec, hksat_z, ierr)
    CHKERRA(ierr);
    call VecRestoreArrayF90(elm_tdycore_idata%watsat_elm_mvec, watsat, ierr)
    CHKERRA(ierr);
    call VecRestoreArrayF90(elm_tdycore_idata%thetares_elm_mvec, thetares, ierr)
    CHKERRA(ierr);

  end subroutine PackMaterialProperties_ELM2TDycore

  !-----------------------------------------------------------------------
  subroutine MapMaterialProperties_ELM2TDycore(this)

    implicit none
    class(em_tdycore_type)  , intent(inout) :: this

    call MappingSourceToDestination(this%map_elm_sub_to_tdy_sub, &
        elm_tdycore_idata%hksat_x_elm_mvec, &
        elm_tdycore_idata%hksat_x_tdycore_svec)

    call MappingSourceToDestination(this%map_elm_sub_to_tdy_sub, &
        elm_tdycore_idata%hksat_y_elm_mvec, &
        elm_tdycore_idata%hksat_y_tdycore_svec)

    call MappingSourceToDestination(this%map_elm_sub_to_tdy_sub, &
        elm_tdycore_idata%hksat_z_elm_mvec, &
        elm_tdycore_idata%hksat_z_tdycore_svec)

    call MappingSourceToDestination(this%map_elm_sub_to_tdy_sub, &
        elm_tdycore_idata%watsat_elm_mvec, &
        elm_tdycore_idata%watsat_tdycore_svec)

    call MappingSourceToDestination(this%map_elm_sub_to_tdy_sub, &
        elm_tdycore_idata%thetares_elm_mvec, &
        elm_tdycore_idata%thetares_tdycore_svec)

  end subroutine MapMaterialProperties_ELM2TDycore

  !-----------------------------------------------------------------------
  subroutine UnpackMaterialProperties_ELM2TDycore(this)

    implicit none
    class(em_tdycore_type)  , intent(inout) :: this

    PetscReal :: den, vis, grav, perm
    PetscReal, pointer :: hksat_x(:), hksat_y(:), hksat_z(:), watsat(:), thetares(:)
    PetscReal, pointer :: blockPerm(:), porosity(:), residualSat(:)
    PetscInt, pointer :: index(:)
    PetscInt :: ngridcells, g
    PetscErrorCode :: ierr

    den = 998.2d0       ! [kg/m^3]  @ 20 degC
    vis = 0.001002d0    ! [N s/m^2] @ 20 degC
    grav = 9.81d0       ! [m/s^2]

    call VecGetSize(elm_tdycore_idata%hksat_x_tdycore_svec,ngridcells,ierr);
    CHKERRA(ierr);

    allocate(blockPerm(ngridcells*9))
    allocate(porosity(ngridcells))
    allocate(residualSat(ngridcells))
    allocate(index(ngridcells))

    call VecGetArrayF90(elm_tdycore_idata%hksat_x_tdycore_svec, hksat_x, ierr)
    call VecGetArrayF90(elm_tdycore_idata%hksat_y_tdycore_svec, hksat_y, ierr)
    call VecGetArrayF90(elm_tdycore_idata%hksat_z_tdycore_svec, hksat_z, ierr)
    call VecGetArrayF90(elm_tdycore_idata%watsat_tdycore_svec, watsat, ierr)
    call VecGetArrayF90(elm_tdycore_idata%thetares_tdycore_svec, thetares, ierr)

    blockPerm(:) = 0.d0

    do g = 1, ngridcells
       blockPerm((g-1)*9+1) = hksat_x(g)*vis/(den*grav)/1000.d0
       blockPerm((g-1)*9+4) = hksat_y(g)*vis/(den*grav)/1000.d0
       blockPerm((g-1)*9+7) = hksat_z(g)*vis/(den*grav)/1000.d0
       porosity(g) = watsat(g)
       residualSat(g) = thetares(g)/porosity(g)
       index(g) = g-1
    enddo

    call VecRestoreArrayF90(elm_tdycore_idata%hksat_x_tdycore_svec, hksat_x, ierr)
    call VecRestoreArrayF90(elm_tdycore_idata%hksat_y_tdycore_svec, hksat_y, ierr)
    call VecRestoreArrayF90(elm_tdycore_idata%hksat_z_tdycore_svec, hksat_z, ierr)
    call VecRestoreArrayF90(elm_tdycore_idata%watsat_tdycore_svec, watsat, ierr)
    call VecRestoreArrayF90(elm_tdycore_idata%thetares_tdycore_svec, thetares, ierr)

    call TDySetPorosityValuesLocal(this%tdy,ngridcells,index,porosity,ierr);
    CHKERRA(ierr);

    call TDySetBlockPermeabilityValuesLocal(this%tdy,ngridcells,index,blockPerm,ierr);
    CHKERRA(ierr);

    call TDySetResidualSaturationValuesLocal(this%tdy,ngridcells,index,residualSat,ierr);
    CHKERRA(ierr);

    deallocate(blockPerm)
    deallocate(porosity)
    deallocate(residualSat)
    deallocate(index)

  end subroutine UnpackMaterialProperties_ELM2TDycore


  !-----------------------------------------------------------------------
  subroutine InitializeMap(this, map)

    use spmdMod, only : mpicom, iam
    use clm_varpar                    , only : nlevsoi, nlevgrnd
    !
    implicit none

    ! !ARGUMENTS:
    class(em_tdycore_type)  , intent(inout) :: this
    type(mapping_type), pointer :: map

    PetscInt :: g, s_npts_local, d_npts_local, d_npts_ghost
    PetscInt, pointer :: s_cell_ids_nindex(:), d_cell_ids_nindex(:), d_is_local(:)

    map%elm_nlevsoi     = nlevsoi
    map%elm_nlevgrnd    = nlevgrnd
    map%elm_nlev_mapped = elm_tdycore_idata%nzelm_mapped
    map%tdy_nlev        = elm_tdycore_idata%nlelm_sub
    map%tdy_nlev_mapped = elm_tdycore_idata%nlelm_sub

    ! call MappingReadTxtFile(map, map%filename)
    call MappingSetupIdentityMap(map, elm_tdycore_idata%nlelm_srf)

    s_npts_local = elm_tdycore_idata%nzelm_mapped * elm_tdycore_idata%nlelm_srf

    allocate(s_cell_ids_nindex(s_npts_local))
    do g = 1, s_npts_local
      s_cell_ids_nindex(g) = g - 1
    enddo

    d_npts_local = s_npts_local
    d_npts_ghost = 0
    allocate(d_cell_ids_nindex(d_npts_local))
    allocate(d_is_local(d_npts_local))
    do g = 1, d_npts_local
      d_cell_ids_nindex(g) = g - 1
      d_is_local(g) = 1
    enddo

    call MappingSetSourceMeshCellIds(map, s_npts_local, s_cell_ids_nindex)
    call MappingSetDestinationMeshCellIds(map, d_npts_local, d_npts_ghost, &
          d_cell_ids_nindex, d_is_local)

    call MappingDecompose(map, mpicom)
    call MappingFindDistinctSourceMeshCellIds(map)
    call MappingCreateWeightMatrix(map, iam)
    call MappingCreateScatterOfSourceMesh(map, mpicom)

    deallocate(s_cell_ids_nindex)
    deallocate(d_cell_ids_nindex)
    deallocate(d_is_local)

  end subroutine InitializeMap

  !-----------------------------------------------------------------------
  subroutine CreateELMTDycoreInterfaceData(this, bounds)
    !
    ! !DESCRIPTION:
    ! Allocates memory for CLM-PFLOTRAN interface
    !
    ! !USES:
    use decompMod                     , only : bounds_type
    use clm_varpar                    , only : nlevsoi, nlevgrnd
    use shr_kind_mod                  , only: r8 => shr_kind_r8
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_tdycore_type)  , intent(inout) :: this
    type(bounds_type)       , intent(in)    :: bounds

    ! Initialize PETSc vector for data transfer between ELM and TDycore
    call ELMTDycore_InterfaceData_Init()

    !
    elm_tdycore_idata%nzelm_mapped = nlevgrnd
    elm_tdycore_idata%nlelm_srf    = bounds%endg-bounds%endg+1

    elm_tdycore_idata%nlelm_sub    = elm_tdycore_idata%nlelm_srf * elm_tdycore_idata%nzelm_mapped
    elm_tdycore_idata%ngelm_sub    = elm_tdycore_idata%nlelm_sub

    elm_tdycore_idata%nltdy_sub    = elm_tdycore_idata%nlelm_sub
    elm_tdycore_idata%ngtdy_sub    = elm_tdycore_idata%nlelm_sub

    ! Allocate vectors for data transfer between CLM and PFLOTRAN.
    call ELMTDycore_InterfaceData_CreateVectors(MPI_COMM_WORLD)

  end subroutine CreateELMTDycoreInterfaceData

  !------------------------------------------------------------------------
   subroutine EM_TDycore_Solve(this, em_stage, dt, nstep, clump_rank, l2e_list, e2l_list, &
        bounds_clump)
    !
    ! !DESCRIPTION:
    ! The TDycore dirver subroutine
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_tdycore_type)               :: this
    integer              , intent(in)    :: em_stage
    real(r8)             , intent(in)    :: dt
    integer              , intent(in)    :: nstep
    integer              , intent(in)    :: clump_rank
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list
    type(bounds_type)    , intent(in)    :: bounds_clump

  end subroutine EM_TDycore_Solve

  !-----------------------------------------------------------------------
  subroutine extract_data_for_elm(this, l2e_init_list, e2l_init_list, bounds_clump)
    !
    !DESCRIPTION
    !  Saves
    !
    use landunit_varcon           , only : istsoil
    use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants , only : VAR_MASS
    use MultiPhysicsProbConstants , only : VAR_SOIL_MATRIX_POT
    use clm_varcon                , only : denice, denh2o
    use clm_varpar                , only : nlevgrnd
    !
    implicit none
    !
    ! !ARGUMENTS
    class(em_tdycore_type)  , intent(inout) :: this
    class(emi_data_list) , intent(in)    :: l2e_init_list
    class(emi_data_list) , intent(inout) :: e2l_init_list
    type(bounds_type)    , intent(in)    :: bounds_clump

    ! !LOCAL VARIABLES:
    integer               :: c,j,g,l,pf_j  ! do loop indices
    integer               :: nlevmapped
    integer               :: gcount
    integer               :: bounds_proc_begc, bounds_proc_endc

    integer     , pointer :: col_active(:)
    integer     , pointer :: col_landunit(:)
    integer     , pointer :: col_gridcell(:)
    integer     , pointer :: lun_type(:)

    real(r8)    , pointer :: l2e_h2osoi_ice(:,:)

    real(r8)    , pointer :: e2l_h2osoi_liq(:,:)
    real(r8)    , pointer :: e2l_h2osoi_ice(:,:)
    real(r8)    , pointer :: e2l_h2osoi_vol(:,:)
    real(r8)    , pointer :: e2l_zwt(:)
    real(r8)    , pointer :: e2l_mflx_snowlyr_col(:)
    real(r8)    , pointer :: e2l_watsatc(:,:)
    real(r8)    , pointer :: e2l_hksatc(:,:)
    real(r8)    , pointer :: e2l_bswc(:,:)
    real(r8)    , pointer :: e2l_sucsatc(:,:)

    real(r8)    , pointer :: dz(:,:)

    PetscScalar , pointer :: sat_loc(:)
    PetscScalar , pointer :: watsat_loc(:)
    PetscScalar , pointer :: hksat_loc(:)
    PetscScalar , pointer :: bsw_loc(:)
    PetscScalar , pointer :: sucsat_loc(:)
    PetscErrorCode        :: ierr

    character(len= 128)   :: subname = 'extract_data_for_elm'
    !-----------------------------------------------------------------------

    bounds_proc_begc = bounds_clump%begc
    bounds_proc_endc = bounds_clump%endc
    nlevmapped       = elm_tdycore_idata%nzelm_mapped

    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_active             , col_active   )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_landunit_index     , col_landunit )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_col_gridcell_index     , col_gridcell )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_col_dz                , dz           )
    call l2e_init_list%GetPointerToInt1D(this%index_l2e_init_landunit_type          , lun_type     )
    call l2e_init_list%GetPointerToReal2D(this%index_l2e_init_h2osoi_ice            , l2e_h2osoi_ice       )

    call e2l_init_list%GetPointerToReal1D(this%index_e2l_init_state_wtd             , e2l_zwt              )
    call e2l_init_list%GetPointerToReal1D(this%index_e2l_init_flux_mflx_snowlyr_col , e2l_mflx_snowlyr_col )

    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_state_h2osoi_liq      , e2l_h2osoi_liq       )
    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_state_h2osoi_ice      , e2l_h2osoi_ice       )
    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_state_h2osoi_vol      , e2l_h2osoi_vol       )

    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_parameter_watsatc     , e2l_watsatc )
    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_parameter_hksatc      , e2l_hksatc  )
    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_parameter_bswc        , e2l_bswc    )
    call e2l_init_list%GetPointerToReal2D(this%index_e2l_init_parameter_sucsatc     , e2l_sucsatc )

    !call pflotranModelGetSoilProp(this%pflotran_m)

    ! Set initial value of for ELM
    e2l_mflx_snowlyr_col(:) = 0._r8
    e2l_zwt(:)              = 0._r8

    ! Initialize soil moisture
    call VecGetArrayF90(elm_tdycore_idata%sat_elm_svec     , sat_loc    , ierr)
    call VecGetArrayF90(elm_tdycore_idata%watsat_elm_svec  , watsat_loc , ierr)
    call VecGetArrayF90(elm_tdycore_idata%hksat_x_elm_svec , hksat_loc  , ierr)
    call VecGetArrayF90(elm_tdycore_idata%bsw_elm_svec     , bsw_loc    , ierr)
    call VecGetArrayF90(elm_tdycore_idata%sucsat_elm_svec  , sucsat_loc , ierr)

    do c = bounds_proc_begc, bounds_proc_endc
       if (col_active(c) == 1) then
          l = col_landunit(c)
          if (lun_type(l) == istsoil) then
             g = col_gridcell(c)
             gcount = g - bounds_clump%begg
             do j = 1, nlevgrnd
                pf_j = gcount*nlevmapped + j

                if (j <= nlevmapped) then
                   e2l_h2osoi_liq(c,j) = sat_loc(pf_j)*watsat_loc(pf_j)*dz(c,j)*1.e3_r8

                   e2l_h2osoi_vol(c,j) = e2l_h2osoi_liq(c,j)/dz(c,j)/denh2o + &
                        l2e_h2osoi_ice(c,j)/dz(c,j)/denice
                   e2l_h2osoi_vol(c,j) = min(e2l_h2osoi_vol(c,j),watsat_loc(pf_j))
                   e2l_h2osoi_ice(c,j) = 0._r8

                   e2l_watsatc(c,j) = watsat_loc(pf_j)
                   e2l_hksatc(c,j)  = hksat_loc(pf_j)
                   e2l_bswc(c,j)    = bsw_loc(pf_j)
                   e2l_sucsatc(c,j) = sucsat_loc(pf_j)
                else
                   e2l_h2osoi_liq(c,j) = e2l_h2osoi_liq(c,nlevmapped)
                   e2l_h2osoi_vol(c,j) = e2l_h2osoi_vol(c,nlevmapped)
                   e2l_h2osoi_ice(c,j) = 0._r8
                   e2l_watsatc(c,j)    = e2l_watsatc(c,nlevmapped)
                   e2l_hksatc(c,j)     = e2l_hksatc(c,nlevmapped)
                   e2l_bswc(c,j)       = e2l_bswc(c,nlevmapped)
                   e2l_sucsatc(c,j)    = e2l_sucsatc(c,nlevmapped)
                end if

             enddo
          else
             write(iulog,*)'WARNING: Land Unit type other than soil type is present within the domain'
             call endrun( trim(subname)//' ERROR: Land Unit type not supported' )
          endif
       endif
    enddo

    call VecRestoreArrayF90(elm_tdycore_idata%sat_elm_svec     , sat_loc    , ierr)
    call VecRestoreArrayF90(elm_tdycore_idata%watsat_elm_svec  , watsat_loc , ierr)
    call VecRestoreArrayF90(elm_tdycore_idata%hksat_x_elm_svec , hksat_loc  , ierr)
    call VecRestoreArrayF90(elm_tdycore_idata%bsw_elm_svec     , bsw_loc    , ierr)
    call VecRestoreArrayF90(elm_tdycore_idata%sucsat_elm_svec  , sucsat_loc , ierr)

   end subroutine extract_data_for_elm

end module ExternalModelTDycore

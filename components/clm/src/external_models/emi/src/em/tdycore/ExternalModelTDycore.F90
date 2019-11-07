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

  use petscsnes
  use petscvec
  use petscdm
  use petscts
  use Mapping_module

#ifdef USE_TDYCORE
  use tdycore
#endif

  implicit none

  PetscInt, parameter :: tdycore_solver_ts = 1
  PetscInt, parameter :: tdycore_solver_snes = 2

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
     SNES :: snes
     Vec :: U
     PetscInt :: tdycore_solver_type
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
    character(len=256)     :: solver_type

    allocate(this%tdy)

    this%tdycore_solver_type = tdycore_solver_ts

    call PetscOptionsGetString (PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
         '-tdycore_solver_type', solver_type, flg, ierr)
    if (flg) then
      select case(trim(solver_type))
      case('ts')
        this%tdycore_solver_type = tdycore_solver_ts
      case('snes')
        this%tdycore_solver_type = tdycore_solver_snes
      case default
        write(iulog,*)'Unknown -tdycore_solver_type ',trim(solver_type)
        call endrun('ERROR: Only ts or snes is supported -tdycore_solver_type' )
      end select
    endif

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

    call DMCreateGlobalVector(this%dm,this%U,ierr);
    CHKERRA(ierr);
    call VecSet(this%U,71325.d0,ierr);
    CHKERRA(ierr);

    select case(this%tdycore_solver_type)
    case (tdycore_solver_ts)
      call TSCreate(PETSC_COMM_WORLD,this%ts,ierr);
      CHKERRA(ierr);

      call TSSetEquationType(this%ts,TS_EQ_IMPLICIT,ierr);
      CHKERRA(ierr);

      call TSSetType(this%ts,TSBEULER,ierr);
      CHKERRA(ierr);

      call TDySetIFunction(this%ts,this%tdy,ierr);
      CHKERRA(ierr);

      call TDySetIJacobian(this%ts,this%tdy,ierr);
      CHKERRA(ierr);

      call TSSetDM(this%ts,this%dm,ierr);
      CHKERRA(ierr);

      call TSSetSolution(this%ts,this%U,ierr);
      CHKERRA(ierr);

      call TSSetMaxSteps(this%ts,10,ierr);
      CHKERRA(ierr);

      call TSSetMaxTime(this%ts,1800.d0,ierr);
      CHKERRA(ierr);

      call TSSetExactFinalTime(this%ts,TS_EXACTFINALTIME_STEPOVER,ierr);
      CHKERRA(ierr);

      call TSSetFromOptions(this%ts,ierr);
      CHKERRA(ierr);

      call TSSetUp(this%ts,ierr);
      CHKERRA(ierr);

    case (tdycore_solver_snes)
      call SNESCreate(PETSC_COMM_WORLD,this%snes,ierr);
      CHKERRA(ierr);

      call TDySetSNESFunction(this%snes,this%tdy,ierr);
      CHKERRA(ierr);

      call TDySetSNESJacobian(this%snes,this%tdy,ierr);
      CHKERRA(ierr);

      call SNESSetFromOptions(this%snes,ierr);
      CHKERRA(ierr);

      call TDySetInitialSolutionForSNESSolver(this%tdy,this%U,ierr);
      CHKERRA(ierr);

    case default
      call endrun('ERROR: TDycore solver setup only ts or snes' )

    end select

    ! ==========================================================================

    ! ==========================================================================
    ! Set initial condition

    call VecGetArrayF90(this%U,p_loc,ierr);
    CHKERRA(ierr);

    call TDyUpdateState(this%tdy,p_loc,ierr);
    CHKERRA(ierr);

    call VecRestoreArrayF90(this%U,p_loc,ierr);
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
    porosity = 0.115d0
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
       blockPerm((g-1)*9+5) = hksat_y(g)*vis/(den*grav)/1000.d0
       blockPerm((g-1)*9+9) = hksat_z(g)*vis/(den*grav)/1000.d0
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

    select case (em_stage)
    case (EM_TDYCORE_SOIL_HYDRO_STAGE)
       call EM_TDYCORE_Solve_Soil_Hydro(this, em_stage, dt, nstep, l2e_list, e2l_list, &
            bounds_clump)

    case default
       write(iulog,*)'EM_TDycore_Solve: Unknown em_stage = ' ,em_stage
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

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

    !------------------------------------------------------------------------
  subroutine EM_TDycore_Solve_Soil_Hydro(this, em_stage, dt, nstep, l2e_list, e2l_list, &
       bounds_clump)
    !
    ! !DESCRIPTION:
    ! Solve the Variably Saturated Flow Model (VSFM) in soil columns.
    !
#include <petsc/finclude/petsc.h>
    !
    ! !USES:
    use shr_kind_mod              , only : r8 => shr_kind_r8
    use abortutils                , only : endrun
    use shr_log_mod               , only : errMsg => shr_log_errMsg
    use MultiPhysicsProbConstants , only : VAR_BC_SS_CONDITION
    use MultiPhysicsProbConstants , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants , only : VAR_PRESSURE
    use MultiPhysicsProbConstants , only : VAR_LIQ_SAT
    use MultiPhysicsProbConstants , only : VAR_FRAC_LIQ_SAT
    use MultiPhysicsProbConstants , only : VAR_MASS
    use MultiPhysicsProbConstants , only : VAR_SOIL_MATRIX_POT
    use MultiPhysicsProbConstants , only : VAR_LATERAL_MASS_EXCHANGED
    use MultiPhysicsProbConstants , only : VAR_BC_MASS_EXCHANGED
    use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants , only : AUXVAR_BC
    use MultiPhysicsProbConstants , only : AUXVAR_SS
    use mpp_varpar                , only : nlevgrnd
    use clm_varcon                , only : denice, denh2o
    use petscsnes
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_tdycore_type)              :: this
    integer              , intent(in)    :: em_stage
    real(r8)             , intent(in)    :: dt
    integer              , intent(in)    :: nstep
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list
    type(bounds_type)    , intent(in)    :: bounds_clump
    !
    ! !LOCAL VARIABLES:
    integer                              :: p,c,fc,j,g                                                       ! do loop indices
    integer                              :: pi                                                               ! pft index
    real(r8)                             :: dzsum                                                            ! summation of dzmm of layers below water table (mm)
    real(r8)                             :: dtime

    real(r8)  , pointer                  :: mflx_et_col_1d         (:)
    real(r8)  , pointer                  :: mflx_infl_col_1d       (:)
    real(r8)  , pointer                  :: mflx_dew_col_1d        (:)
    real(r8)  , pointer                  :: mflx_drain_col_1d      (:)
    real(r8)  , pointer                  :: mflx_sub_snow_col_1d   (:)
    real(r8)  , pointer                  :: mflx_snowlyr_col_1d    (:)
    real(r8)  , pointer                  :: t_soil_col_1d          (:)
    integer    , pointer                 :: col_active             (:)
    integer    , pointer                 :: col_gridcell           (:)

    real(r8)  , pointer                  :: fliq_col_1d       (:)
    real(r8)  , pointer                  :: mass_col_1d       (:)
    real(r8)  , pointer                  :: smpl_col_1d       (:)
    real(r8)  , pointer                  :: soilp_col_1d      (:)
    real(r8)  , pointer                  :: sat_col_1d        (:)

    real(r8)  , pointer                  :: frac_ice                    (:,:) ! fraction of ice
    real(r8)  , pointer                  :: total_mass_flux_col         (:)            ! Sum of all source-sinks conditions for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_et_col      (:)            ! ET sink for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_infl_col    (:)            ! Infiltration source for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_dew_col     (:)            ! Dew source for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_drain_col   (:)            ! Drainage sink for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_snowlyr_col (:)            ! Flux due to disappearance of snow for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_sub_col     (:)            ! Sublimation sink for VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_lateral_col (:)            ! Lateral flux computed by VSFM solver at column level
    real(r8)  , pointer                  :: total_mass_flux_seepage_col (:)            ! Seepage flux computed by VSFM solver at column level
    real(r8)  , pointer                  :: qflx_seepage                (:)            ! Seepage flux computed by VSFM solver at column level
    real(r8)  , pointer                  :: mass_prev_col          (:,:) ! Mass of water before a VSFM solve
    real(r8)  , pointer                  :: dmass_col              (:)            ! Change in mass of water after a VSFM solve
    real(r8)  , pointer                  :: mass_beg_col                (:)            ! Total mass before a VSFM solve
    real(r8)  , pointer                  :: mass_end_col                (:)            ! Total mass after a VSFM solve
    integer                              :: ier                                                              ! error status

    integer                              :: begc, endc
    integer                              :: g_idx, c_idx

    PetscInt                             :: soe_auxvar_id                                                    ! Index of system-of-equation's (SoE's) auxvar
    PetscErrorCode                       :: ierr                                                             ! PETSc return error code

    PetscBool                            :: converged                                                        ! Did VSFM solver converge to a solution with given PETSc SNES tolerances
    PetscInt                             :: converged_reason                                                 ! SNES converged due to which criteria
    PetscReal                            :: atol_default                                                     ! Default SNES absolute convergance tolerance
    PetscReal                            :: rtol_default                                                     ! Default SNES relative convergance tolerance
    PetscReal                            :: stol_default                                                     ! Default SNES solution convergance tolerance
    PetscInt                             :: max_it_default                                                   ! Default SNES maximum number of iteration
    PetscInt                             :: max_f_default                                                    ! Default SNES maximum number of function evaluation
    PetscReal                            :: stol                                                             ! solution convergance tolerance
    PetscReal                            :: rtol                                                             ! relative convergance tolerance
    PetscReal,parameter                  :: stol_alternate = 1.d-10                                          ! Alternate solution convergance tolerance

    PetscReal                            :: mass_beg                                                         ! Sum of mass of water for all active soil columns before VSFM is called
    PetscReal                            :: mass_end                                                         ! Sum of mass of water for all active soil columns after VSFM is called
    PetscReal                            :: total_mass_flux_et                                               ! Sum of mass ET mass flux of water for all active soil columns
    PetscReal                            :: total_mass_flux_infl                                             ! Sum of mass infiltration mass flux of water for all active soil columns
    PetscReal                            :: total_mass_flux_dew                                              ! Sum of mass dew mass flux of water for all active soil columns
    PetscReal                            :: total_mass_flux_drain                                            ! Sum of mass drainage mass flux of water for all active soil columns
    PetscReal                            :: total_mass_flux_snowlyr                                          ! Sum of mass snow layer disappearance mass flux of water for all active soil columns
    PetscReal                            :: total_mass_flux_sub                                              ! Sum of mass sublimation mass flux of water for all active soil columns
    PetscReal                            :: total_mass_flux_lateral                                          ! Sum of lateral mass flux for all active soil columns
    PetscReal                            :: total_mass_flux                                                  ! Sum of mass ALL mass flux of water for all active soil columns
    PetscInt                             :: iter_count                                                       ! How many times VSFM solver is called

    PetscInt, parameter                  :: max_iter_count = 10                                              ! Maximum number of times VSFM can be called
    PetscInt                             :: diverged_count                                                   ! Number of time VSFM solver diverged
    PetscInt                             :: mass_bal_err_count                                               ! Number of time VSFM solver returns a solution that isn't within acceptable mass balance error threshold
    PetscReal                            :: abs_mass_error_col                                               ! Maximum absolute error for any active soil column
    PetscReal, parameter                 :: max_abs_mass_error_col  = 1.e-5                                  ! Acceptable mass balance error
    PetscBool                            :: successful_step                                                  ! Is the solution return by VSFM acceptable
    PetscReal , pointer                  :: soilp_col_ghosted_1d(:)
    PetscReal , pointer                  :: fliq_col_ghosted_1d(:)
    PetscReal , pointer                  :: mflx_lateral_col_1d(:)
    PetscReal , pointer                  :: lat_mass_exc_col_1d(:)
    PetscReal , pointer                  :: seepage_mass_exc_col_1d(:)
    PetscReal , pointer                  :: seepage_press_1d(:)

    integer                              :: jwt
    real(r8)                             :: z_dn, z_up

    real(r8)  , pointer                  :: l2e_mflux_infil(:)
    real(r8)  , pointer                  :: l2e_mflux_dew(:)
    real(r8)  , pointer                  :: l2e_mflux_sub_snow(:)
    real(r8)  , pointer                  :: l2e_mflux_snowlyr(:)
    real(r8)  , pointer                  :: l2e_mflux_et(:,:)
    real(r8)  , pointer                  :: l2e_mflux_drain(:,:)
    real(r8)  , pointer                  :: l2e_h2osoi_liq(:,:)
    real(r8)  , pointer                  :: l2e_h2osoi_ice(:,:)
    real(r8)  , pointer                  :: l2e_zi(:,:)
    real(r8)  , pointer                  :: col_dz(:,:)
    integer   , pointer                  :: l2e_filter_hydrologyc(:)
    integer                              :: l2e_num_hydrologyc

    real(r8)  , pointer                  :: e2l_h2osoi_liq(:,:)
    real(r8)  , pointer                  :: e2l_h2osoi_ice(:,:)
    real(r8)  , pointer                  :: e2l_smp(:,:)
    real(r8)  , pointer                  :: e2l_wtd(:)
    real(r8)  , pointer                  :: e2l_soilp(:,:)
    real(r8)  , pointer                  :: e2l_qrecharge(:)

    PetscScalar, pointer :: qflx_elm_loc(:)
    PetscScalar, pointer :: area_clm_loc(:)
    PetscScalar, pointer :: e2l_drain_perched(:)
    PetscScalar, pointer :: e2l_drain(:)
    PetscScalar, pointer :: e2l_qrgwl(:)
    PetscScalar, pointer :: e2l_rsub_sat(:)

    integer :: bounds_proc_begc, bounds_proc_endc
    integer :: nlevmapped
    real(r8) :: col_wtgcell
    real(r8) :: area
    integer :: ngridcells, nvalues
    integer, pointer :: index(:)
    PetscReal, pointer :: liquid_mass(:), mass_p(:), p_loc(:), qflx_p(:)
    SNESConvergedReason :: snes_reason

    !-----------------------------------------------------------------------

    bounds_proc_begc     = bounds_clump%begc
    bounds_proc_endc     = bounds_clump%endc

    ! Get time step

    dtime = dt

    call l2e_list%GetPointerToReal1D(this%index_l2e_flux_infil       , l2e_mflux_infil       )
    call l2e_list%GetPointerToReal1D(this%index_l2e_flux_dew         , l2e_mflux_dew         )
    call l2e_list%GetPointerToReal1D(this%index_l2e_flux_snow_sub    , l2e_mflux_sub_snow    )
    call l2e_list%GetPointerToReal1D(this%index_l2e_flux_snowlyr     , l2e_mflux_snowlyr     )

    call l2e_list%GetPointerToReal2D(this%index_l2e_flux_et          , l2e_mflux_et          )
    call l2e_list%GetPointerToReal2D(this%index_l2e_flux_drainage    , l2e_mflux_drain       )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_liq , l2e_h2osoi_liq        )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_ice , l2e_h2osoi_ice        )

    call l2e_list%GetPointerToInt1D(this%index_l2e_filter_hydrologyc , l2e_filter_hydrologyc )
    call l2e_list%GetIntValue(this%index_l2e_filter_num_hydrologyc   , l2e_num_hydrologyc    )

    call l2e_list%GetPointerToReal2D(this%index_l2e_column_zi        , l2e_zi                )
    call l2e_list%GetPointerToInt1D(this%index_l2e_col_active        , col_active            )
    call l2e_list%GetPointerToInt1D(this%index_l2e_col_gridcell      , col_gridcell          )
    call l2e_list%GetPointerToReal2D(this%index_l2e_init_col_dz      , col_dz                )

    call e2l_list%GetPointerToReal1D(this%index_e2l_state_wtd        , e2l_wtd               )
    call e2l_list%GetPointerToReal2D(this%index_e2l_state_h2osoi_liq , e2l_h2osoi_liq        )
    call e2l_list%GetPointerToReal2D(this%index_e2l_state_h2osoi_ice , e2l_h2osoi_ice        )
    !call e2l_list%GetPointerToReal2D(this%index_e2l_state_smp        , e2l_smp               )
    call e2l_list%GetPointerToReal2D(this%index_e2l_state_soilp      , e2l_soilp             )

    call e2l_list%GetPointerToReal1D(this%index_e2l_flux_qrecharge    , e2l_qrecharge        )
    call e2l_list%GetPointerToReal1D(this%index_e2l_flux_drain_perched, e2l_drain_perched    )
    call e2l_list%GetPointerToReal1D(this%index_e2l_flux_drain        , e2l_drain            )
    call e2l_list%GetPointerToReal1D(this%index_e2l_flux_qrgwl        , e2l_qrgwl            )
    call e2l_list%GetPointerToReal1D(this%index_e2l_flux_rsub_sat     , e2l_rsub_sat         )

    begc = bounds_proc_begc
    endc = bounds_proc_endc

    allocate(frac_ice                    (begc:endc,1:nlevgrnd))
    allocate(total_mass_flux_col         (begc:endc))
    allocate(total_mass_flux_et_col      (begc:endc))
    allocate(total_mass_flux_infl_col    (begc:endc))
    allocate(total_mass_flux_dew_col     (begc:endc))
    allocate(total_mass_flux_drain_col   (begc:endc))
    allocate(total_mass_flux_snowlyr_col (begc:endc))
    allocate(total_mass_flux_sub_col     (begc:endc))
    allocate(total_mass_flux_lateral_col (begc:endc))
    allocate(total_mass_flux_seepage_col (begc:endc))
    allocate(qflx_seepage                (begc:endc))
    allocate(mass_prev_col          (begc:endc,1:nlevgrnd))
    allocate(dmass_col              (begc:endc))
    allocate(mass_beg_col                (begc:endc))
    allocate(mass_end_col                (begc:endc))

    allocate(mflx_et_col_1d              ((endc-begc+1)*nlevgrnd))
    allocate(mflx_drain_col_1d           ((endc-begc+1)*nlevgrnd))
    allocate(mflx_infl_col_1d            (endc-begc+1))
    allocate(mflx_dew_col_1d             (endc-begc+1))
    allocate(mflx_sub_snow_col_1d        (endc-begc+1))
    allocate(mflx_snowlyr_col_1d         (endc-begc+1))
    allocate(t_soil_col_1d               ((endc-begc+1)*nlevgrnd))

    allocate(mass_col_1d            ((endc-begc+1)*nlevgrnd))
    allocate(fliq_col_1d            ((endc-begc+1)*nlevgrnd))
    allocate(smpl_col_1d            ((endc-begc+1)*nlevgrnd))
    allocate(soilp_col_1d           ((endc-begc+1)*nlevgrnd))
    allocate(sat_col_1d             ((endc-begc+1)*nlevgrnd))

    ! initialize
    mflx_et_col_1d(:)                = 0.d0
    mflx_infl_col_1d(:)              = 0.d0
    mflx_dew_col_1d(:)               = 0.d0
    mflx_drain_col_1d(:)             = 0.d0
    mflx_sub_snow_col_1d(:)          = 0.d0
    mflx_snowlyr_col_1d(:)           = 0.d0
    t_soil_col_1d(:)                 = 298.15d0

    mass_beg                         = 0.d0
    mass_end                         = 0.d0
    total_mass_flux                  = 0.d0
    total_mass_flux_et               = 0.d0
    total_mass_flux_infl             = 0.d0
    total_mass_flux_dew              = 0.d0
    total_mass_flux_drain            = 0.d0
    total_mass_flux_snowlyr          = 0.d0
    total_mass_flux_sub              = 0.d0
    total_mass_flux_lateral          = 0.d0

    mass_beg_col(:)                  = 0.d0
    mass_end_col(:)                  = 0.d0
    total_mass_flux_col(:)           = 0.d0
    total_mass_flux_et_col(:)        = 0.d0
    total_mass_flux_infl_col(:)      = 0.d0
    total_mass_flux_dew_col(:)       = 0.d0
    total_mass_flux_drain_col(:)     = 0.d0
    total_mass_flux_snowlyr_col(:)   = 0.d0
    total_mass_flux_sub_col(:)       = 0.d0
    total_mass_flux_lateral_col(:)   = 0.d0

    mass_prev_col(:,:)          = 0.d0
    dmass_col(:)                = 0.d0

    nlevmapped = elm_tdycore_idata%nzelm_mapped

    ! Get total mass
    call VecGetSize(elm_tdycore_idata%mass_elm_svec,ngridcells,ierr);
    CHKERRA(ierr);

    call VecGetArrayF90(this%U,p_loc,ierr);
    CHKERRA(ierr);
    call TDyUpdateState(this%tdy,p_loc,ierr);
    CHKERRA(ierr);
    call VecRestoreArrayF90(this%U,p_loc,ierr);
    CHKERRA(ierr);

    allocate(liquid_mass(ngridcells))

    call TDyGetLiquidMassValuesLocal(this%tdy,nvalues,liquid_mass,ierr)
    CHKERRA(ierr);

    call VecGetArrayF90(elm_tdycore_idata%mass_tdycore_mvec, mass_p, ierr)
    CHKERRA(ierr);
    do g = 1, ngridcells
      mass_p(g) = liquid_mass(g)
    enddo
    deallocate(liquid_mass)
    call VecRestoreArrayF90(elm_tdycore_idata%mass_tdycore_mvec, mass_p, ierr)
    CHKERRA(ierr);

    call MappingSourceToDestination(this%map_tdy_sub_to_elm_sub, &
        elm_tdycore_idata%mass_tdycore_mvec, &
        elm_tdycore_idata%mass_elm_svec)

    call VecGetArrayF90(elm_tdycore_idata%mass_elm_svec, mass_p, ierr); CHKERRQ(ierr)

    area = 1._r8
    do fc = 1, l2e_num_hydrologyc
       c = l2e_filter_hydrologyc(fc)
       g = col_gridcell(c)

       do j = 1, nlevmapped

          c_idx = (c - begc)*nlevgrnd + j
          g_idx = (g - bounds_clump%begg)*nlevmapped + j

          mflx_et_col_1d(c_idx)          = l2e_mflux_et(c,j)
          mflx_drain_col_1d(c_idx)       = l2e_mflux_drain(c,j)

          total_mass_flux_et           = total_mass_flux_et           + mflx_et_col_1d(c_idx)
          total_mass_flux_et_col(c)    = total_mass_flux_et_col(c)    + mflx_et_col_1d(c_idx)

          total_mass_flux_drain        = total_mass_flux_drain        + mflx_drain_col_1d(c_idx)
          total_mass_flux_drain_col(c) = total_mass_flux_drain_col(c) + mflx_drain_col_1d(c_idx)

          mass_beg                     = mass_beg                     + mass_p(g_idx)!/area_clm_loc(g_idx)
          mass_beg_col(c)              = mass_beg_col(c)              + mass_p(g_idx)!/area_clm_loc(g_idx)
          mass_prev_col(c,j)      = mass_col_1d(c_idx)
       end do

       c_idx = c - begc+1

       mflx_dew_col_1d(c_idx)        = l2e_mflux_dew(c)
       mflx_infl_col_1d(c_idx)       = l2e_mflux_infil(c)
       mflx_snowlyr_col_1d(c_idx)    = l2e_mflux_snowlyr(c)
       mflx_sub_snow_col_1d(c_idx)   = l2e_mflux_sub_snow(c)

       total_mass_flux_dew            = total_mass_flux_dew            + mflx_dew_col_1d(c_idx)
       total_mass_flux_dew_col(c)     = total_mass_flux_dew_col(c)     + mflx_dew_col_1d(c_idx)

       total_mass_flux_infl           = total_mass_flux_infl           + mflx_infl_col_1d(c_idx)
       total_mass_flux_infl_col(c)    = total_mass_flux_infl_col(c)    + mflx_infl_col_1d(c_idx)

       total_mass_flux_snowlyr        = total_mass_flux_snowlyr        + mflx_snowlyr_col_1d(c_idx)
       total_mass_flux_snowlyr_col(c) = total_mass_flux_snowlyr_col(c) + mflx_snowlyr_col_1d(c_idx)

       total_mass_flux_sub            = total_mass_flux_sub            + mflx_sub_snow_col_1d(c_idx)
       total_mass_flux_sub_col(c)     = total_mass_flux_sub_col(c)     + mflx_sub_snow_col_1d(c_idx)

       total_mass_flux_col(c) = total_mass_flux_et_col(c)      + &
            total_mass_flux_infl_col(c)    + &
            total_mass_flux_dew_col(c)     + &
            total_mass_flux_drain_col(c)   + &
            total_mass_flux_snowlyr_col(c) + &
            total_mass_flux_sub_col(c)     + &
            total_mass_flux_lateral_col(c)
    end do
    total_mass_flux        = &
         total_mass_flux_et        + &
         total_mass_flux_infl      + &
         total_mass_flux_dew       + &
         total_mass_flux_drain     + &
         total_mass_flux_snowlyr   + &
         total_mass_flux_sub       + &
         total_mass_flux_lateral
    call VecRestoreArrayF90(elm_tdycore_idata%mass_elm_svec, mass_p, ierr); CHKERRQ(ierr)

    call VecGetArrayF90(elm_tdycore_idata%qflx_elm_mvec, qflx_elm_loc, ierr); CHKERRQ(ierr)

    frac_ice(:,:)       = 0.d0
    do fc = 1, l2e_num_hydrologyc
       c = l2e_filter_hydrologyc(fc)
       do j = 1, nlevmapped
          frac_ice(c,j) = l2e_h2osoi_ice(c,j)/(l2e_h2osoi_liq(c,j) + l2e_h2osoi_ice(c,j))
       end do
    end do

    ! Initialize ET sink to ZERO
    do g = bounds_clump%begg, bounds_clump%endg
       do j = 1,nlevmapped
          g_idx = (g - bounds_clump%begg)*nlevmapped + j
          qflx_elm_loc(g_idx) = 0.0_r8
       end do
    end do

    ! Account for following fluxes in the top soil layer:
    ! - infiltration
    ! - dew
    ! - disapperance of snow layer
    ! - sublimation of snow
    j = 1
    col_wtgcell = 1._r8
    do c = bounds_proc_begc, bounds_proc_endc
       if (col_active(c) == 1) then
          ! Set gridcell indices
          g = col_gridcell(c)
          g_idx = (g - bounds_clump%begg)*nlevmapped + j
          c_idx = c - begc+1
          qflx_elm_loc(g_idx) = qflx_elm_loc(g_idx) + &
               (&
               mflx_infl_col_1d(c_idx)    + &
               mflx_dew_col_1d(c_idx)     + &
               mflx_snowlyr_col_1d(c_idx) + &
               mflx_sub_snow_col_1d(c_idx)  &
               )*col_wtgcell!*area_clm_loc(g_idx)
       end if
       total_mass_flux_col(c) = 0.d0
    enddo

    ! Account for evapotranspiration flux
    do c = bounds_proc_begc, bounds_proc_endc
       do j = 1,nlevmapped
          g = col_gridcell(c)
          g_idx = (g - bounds_clump%begg)*nlevmapped + j
          c_idx = (c - bounds_proc_begc)*nlevgrnd+j
          if (col_active(c) == 1) then
             qflx_elm_loc(g_idx) = qflx_elm_loc(g_idx) + &
                  mflx_et_col_1d(c_idx)!*area_clm_loc(g_idx)
             total_mass_flux_col(c) = total_mass_flux_col(c) + qflx_elm_loc(g_idx)!/area_clm_loc(g_idx)
          end if
       end do
    end do

    call VecRestoreArrayF90(elm_tdycore_idata%qflx_elm_mvec, qflx_elm_loc, ierr); CHKERRQ(ierr)

    call MappingSourceToDestination(this%map_elm_sub_to_tdy_sub, &
        elm_tdycore_idata%qflx_elm_mvec, &
        elm_tdycore_idata%qflx_tdycore_svec)

    call VecGetSize(elm_tdycore_idata%qflx_tdycore_svec,ngridcells,ierr);
    call VecGetArrayF90(elm_tdycore_idata%qflx_tdycore_svec,qflx_p,ierr); CHKERRQ(ierr)
    allocate(index(ngridcells))
    do g = 1, ngridcells
       index(g) = g-1
    enddo
    !qflx_p(:) = 0.d0
    call TDySetSourceSinkValuesLocal(this%tdy,ngridcells,index,qflx_p,ierr)
    call VecRestoreArrayF90(elm_tdycore_idata%qflx_tdycore_svec,qflx_p,ierr); CHKERRQ(ierr)
    deallocate(index)

    select case(this%tdycore_solver_type)
    case (tdycore_solver_ts)
      call TSSetTime(this%ts,0.d0,ierr);CHKERRQ(ierr)
      call TSSetTimeStep(this%ts,dt,ierr);CHKERRQ(ierr)
      call TSSetMaxTime(this%ts,dt,ierr);
      call TSSetStepNumber(this%ts,0,ierr);CHKERRQ(ierr)
      call TSSetExactFinalTime(this%ts,TS_EXACTFINALTIME_MATCHSTEP,ierr);CHKERRQ(ierr)
      call TSSolve(this%ts,this%U,ierr);

    case (tdycore_solver_snes)
      dtime = 1800.0
      call TDySetDtimeForSNESSolver(this%tdy,dtime,ierr);
      CHKERRA(ierr);

      call TDyPreSolveSNESSolver(this%tdy,ierr);
      CHKERRA(ierr);

      call SNESSolve(this%snes,PETSC_NULL_VEC,this%U,ierr);
      CHKERRA(ierr);

      call SNESGetConvergedReason(this%snes, snes_reason, ierr);CHKERRQ(ierr)
      CHKERRA(ierr);
      if (snes_reason<0) then
         call endrun('ERROR: TDycore snes failed to converged' )
      endif

      call TDyPostSolveSNESSolver(this%tdy,this%U,ierr);
      CHKERRA(ierr);


    case default
      call endrun('ERROR: TDycore solver solve only ts or snes' )
    end select

    allocate(liquid_mass(ngridcells))

    call TDyGetLiquidMassValuesLocal(this%tdy,nvalues,liquid_mass,ierr)
    CHKERRA(ierr);

    call VecGetArrayF90(elm_tdycore_idata%mass_tdycore_mvec, mass_p, ierr)
    CHKERRA(ierr);
    do g = 1, ngridcells
      mass_p(g) = liquid_mass(g)
    enddo
    deallocate(liquid_mass)
    call VecRestoreArrayF90(elm_tdycore_idata%mass_tdycore_mvec, mass_p, ierr)
    CHKERRA(ierr);

    call MappingSourceToDestination(this%map_tdy_sub_to_elm_sub, &
        elm_tdycore_idata%mass_tdycore_mvec, &
        elm_tdycore_idata%mass_elm_svec)

    call VecGetArrayF90(elm_tdycore_idata%mass_elm_svec, mass_p, ierr); CHKERRQ(ierr)

    abs_mass_error_col = 0.d0
    do fc = 1, l2e_num_hydrologyc
       c = l2e_filter_hydrologyc(fc)
       g = col_gridcell(c)

       ! initialization
       jwt = -1

       ! Loops in decreasing j so WTD can be computed in the same loop
       e2l_h2osoi_liq(c,:) = 0._r8
       e2l_h2osoi_ice(c,:) = 0._r8
       do j = nlevmapped, 1, -1
          g_idx = (g - bounds_clump%begg)*nlevmapped + j

          e2l_h2osoi_liq(c,j) =  (1.d0 - frac_ice(c,j))*mass_p(g_idx)!/area_clm_loc(g_idx)
          e2l_h2osoi_ice(c,j) =  frac_ice(c,j)         *mass_p(g_idx)!/area_clm_loc(g_idx)

          mass_end        = mass_end        + mass_p(g_idx)!/area_clm_loc(g_idx)
          mass_end_col(c) = mass_end_col(c) + mass_p(g_idx)!/area_clm_loc(g_idx)

       end do

       ! Find maximum water balance error over the column
       abs_mass_error_col = max(abs_mass_error_col,                     &
            abs(mass_beg_col(c) - mass_end_col(c) + &
            total_mass_flux_col(c)*dt))
       e2l_qrecharge     (c) = 0._r8

       e2l_wtd(c) = l2e_zi(c,nlevmapped)
    end do

    ! Save soil liquid pressure from VSFM for all (active+nonactive) cells.
    ! soilp_col is used for restarting VSFM.
    do c = begc, endc
       do j = 1, nlevgrnd
          c_idx = (c - begc)*nlevgrnd + j
          e2l_soilp(c,j) = 0._r8!soilp_col_1d(c_idx)
       end do
    end do

    e2l_drain_perched (:) = 0._r8
    e2l_drain         (:) = 0._r8
    e2l_qrgwl         (:) = 0._r8
    e2l_rsub_sat      (:) = 0._r8

    call VecRestoreArrayF90(elm_tdycore_idata%mass_elm_svec, mass_p, ierr); CHKERRQ(ierr)

    write(*,*)'abs_mass_error_col = ',abs_mass_error_col


    deallocate(frac_ice                    )
    deallocate(total_mass_flux_col         )
    deallocate(total_mass_flux_et_col      )
    deallocate(total_mass_flux_infl_col    )
    deallocate(total_mass_flux_dew_col     )
    deallocate(total_mass_flux_drain_col   )
    deallocate(total_mass_flux_snowlyr_col )
    deallocate(total_mass_flux_sub_col     )
    deallocate(total_mass_flux_lateral_col )
    deallocate(total_mass_flux_seepage_col )
    deallocate(qflx_seepage                )
    deallocate(mass_prev_col          )
    deallocate(dmass_col              )
    deallocate(mass_beg_col                )
    deallocate(mass_end_col                )

    deallocate(mflx_et_col_1d              )
    deallocate(mflx_drain_col_1d           )
    deallocate(mflx_infl_col_1d            )
    deallocate(mflx_dew_col_1d             )
    deallocate(mflx_sub_snow_col_1d        )
    deallocate(mflx_snowlyr_col_1d         )
    deallocate(t_soil_col_1d               )

    deallocate(mass_col_1d            )
    deallocate(fliq_col_1d            )
    deallocate(smpl_col_1d            )
    deallocate(soilp_col_1d           )
    deallocate(sat_col_1d             )

  end subroutine EM_TDycore_Solve_Soil_Hydro

end module ExternalModelTDycore

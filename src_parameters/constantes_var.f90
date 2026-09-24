!
!=================================================================================================================================
! MODULE       : constantes_var
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        constantes_var module contains most constantes like pi, Earth radius, etc...
!!              and all externalized parameters except pft-dependent constants.
!!
!!\n DESCRIPTION: This module contains most constantes and the externalized parameters of ORCHIDEE which 
!!                are not pft-dependent.\n
!!                In this module, you can set the flag diag_qsat in order to detect the pixel where the
!!                temperature is out of range (look qsatcalc and dev_qsatcalc in qsat_moisture.f90).\n
!!                The Earth radius is approximated by the Equatorial radius.The Earth's equatorial radius a,
!!                or semi-major axis, is the distance from its center to the equator and equals 6,378.1370 km.
!!                The equatorial radius is often used to compare Earth with other planets.\n
!!                The meridional mean is well approximated by the semicubic mean of the two axe yielding 
!!                6367.4491 km or less accurately by the quadratic mean of the two axes about 6,367.454 km
!!                or even just the mean of the two axes about 6,367.445 km.\n
!!                This module is already USE in module constantes. Therefor no need to USE it seperatly except
!!                if the subroutines in module constantes are not needed.\n
!!                
!! RECENT CHANGE(S):
!!
!! REFERENCE(S)	: 
!! - Louis, Jean-Francois (1979), A parametric model of vertical eddy fluxes in the atmosphere. 
!! Boundary Layer Meteorology, 187-202.\n
!!
!! SVN          :
!! $HeadURL: $
!! $Date: 2026-07-21 10:28:15 +0200 (Tue, 21 Jul 2026) $
!! $Revision: 9660 $
!! \n
!_ ================================================================================================================================

MODULE constantes_var

  USE defprec

  IMPLICIT NONE
!-

                         !-----------------------!
                         !  ORCHIDEE CONSTANTS   !
                         !-----------------------!

  !
  ! FLAGS 
  !
  LOGICAL, SAVE :: nc_restart_compression !! activate netcdf restart compression
!$OMP THREADPRIVATE(nc_restart_compression)
  LOGICAL :: river_routing      !! activate river routing
!$OMP THREADPRIVATE(river_routing)
  LOGICAL, SAVE :: ok_nudge_mc  !! Activate nudging of soil moisture 
!$OMP THREADPRIVATE(ok_nudge_mc)
  LOGICAL, SAVE :: ok_nudge_snow!! Activate nudging of snow variables
!$OMP THREADPRIVATE(ok_nudge_snow)
  LOGICAL, SAVE :: nudge_interpol_with_xios  !! Activate reading and interpolation with XIOS for nudging fields
!$OMP THREADPRIVATE(nudge_interpol_with_xios)
  LOGICAL :: do_floodplains     !! activate flood plains
!$OMP THREADPRIVATE(do_floodplains)
  LOGICAL :: do_irrigation      !! activate computation of irrigation flux
!$OMP THREADPRIVATE(do_irrigation)
  LOGICAL :: ok_sechiba         !! activate physic of the model
!$OMP THREADPRIVATE(ok_sechiba)
  LOGICAL :: ok_stomate         !! activate carbon cycle
!$OMP THREADPRIVATE(ok_stomate)
  LOGICAL :: ok_ncycle          !! activate nitrogen cycle
!$OMP THREADPRIVATE(ok_ncycle)
  LOGICAL :: impose_cn          !! impose the CN ratio of leaves
!$OMP THREADPRIVATE(impose_cn)
  LOGICAL :: reset_impose_cn    !! reset the CN ratio of leaves
!$OMP THREADPRIVATE(reset_impose_cn)
  LOGICAL :: read_cn          !! read the CN ratio of leaves
!$OMP THREADPRIVATE(read_cn)
  LOGICAL :: manure_orga          !! 1 if manure goes to litter pools
!$OMP THREADPRIVATE(manure_orga)
  LOGICAL :: ok_dgvm            !! activate dynamic vegetation
!$OMP THREADPRIVATE(ok_dgvm)
  LOGICAL :: do_wood_harvest    !! activate wood harvest
!$OMP THREADPRIVATE(do_wood_harvest)
  LOGICAL :: ok_pheno           !! activate the calculation of lai using stomate rather than a prescription
!$OMP THREADPRIVATE(ok_pheno)
  LOGICAL :: ok_bvoc            !! activate biogenic volatile organic coumpounds
!$OMP THREADPRIVATE(ok_bvoc)
  LOGICAL :: ok_leafage         !! activate leafage
!$OMP THREADPRIVATE(ok_leafage)
  LOGICAL :: ok_radcanopy       !! use canopy radiative transfer model
!$OMP THREADPRIVATE(ok_radcanopy)
  LOGICAL :: ok_multilayer      !! use canopy radiative transfer model with multi-layers
!$OMP THREADPRIVATE(ok_multilayer)
  LOGICAL :: ok_pulse_NOx       !! calculate NOx emissions with pulse
!$OMP THREADPRIVATE(ok_pulse_NOx)
  LOGICAL :: ok_bbgfertil_NOx   !! calculate NOx emissions with bbg fertilizing effect
!$OMP THREADPRIVATE(ok_bbgfertil_NOx)
  LOGICAL :: ok_cropsfertil_NOx !! calculate NOx emissions with fertilizers use
!$OMP THREADPRIVATE(ok_cropsfertil_NOx)

  LOGICAL :: ok_co2bvoc_poss    !! CO2 inhibition on isoprene activated following Possell et al. (2005) model
!$OMP THREADPRIVATE(ok_co2bvoc_poss)

  LOGICAL :: ok_co2bvoc_wilk    !! CO2 inhibition on isoprene activated following Wilkinson et al. (2006) model
!$OMP THREADPRIVATE(ok_co2bvoc_wilk)

  LOGICAL :: read_ONfert= .FALSE.    !! ONfert data for stics data fertilziation for crops
!$OMP THREADPRIVATE(read_ONfert)  
  
  LOGICAL, SAVE :: OFF_LINE_MODE = .FALSE.  !! ORCHIDEE detects if it is coupled with a GCM or 
                                            !! just use with one driver in OFF-LINE. (true/false)
!$OMP THREADPRIVATE(OFF_LINE_MODE)  
  LOGICAL, SAVE :: impose_param = .TRUE.    !! Flag impos_param : read all the parameters in the run.def file
!$OMP THREADPRIVATE(impose_param)
  CHARACTER(LEN=80), SAVE     :: restname_in       = 'NONE'                 !! Input Restart files name for Sechiba component  
!$OMP THREADPRIVATE(restname_in)
  CHARACTER(LEN=80), SAVE     :: restname_out      = 'sechiba_rest_out.nc'  !! Output Restart files name for Sechiba component
!$OMP THREADPRIVATE(restname_out)
  CHARACTER(LEN=80), SAVE     :: stom_restname_in  = 'NONE'                 !! Input Restart files name for Stomate component
!$OMP THREADPRIVATE(stom_restname_in)
  CHARACTER(LEN=80), SAVE     :: stom_restname_out = 'stomate_rest_out.nc'  !! Output Restart files name for Stomate component
!$OMP THREADPRIVATE(stom_restname_out)
  INTEGER, SAVE :: printlev=2       !! Standard level for text output [0, 1, 2, 3]
!$OMP THREADPRIVATE(printlev)

  !
  ! TIME
  !
  INTEGER(i_std), PARAMETER  :: spring_days_max = 40  !! Maximum number of days during which we watch for possible spring frost damage
  !
  ! SPECIAL VALUES 
  !
  INTEGER(i_std), PARAMETER :: undef_int = 999999999     !! undef integer for integer arrays (unitless)
  !-
  REAL(r_std), SAVE :: val_exp = 999999.                 !! Specific value if no restart value  (unitless)
!$OMP THREADPRIVATE(val_exp)
  REAL(r_std), PARAMETER :: undef = -9999.               !! Special value for stomate (unitless)
  
  REAL(r_std), PARAMETER :: min_sechiba = 1.E-8_r_std    !! Epsilon to detect a near zero floating point (unitless)
  REAL(r_std), PARAMETER :: undef_sechiba = 1.E+20_r_std !! The undef value used in SECHIBA (unitless)
  
  REAL(r_std), PARAMETER :: min_stomate = 1.E-8_r_std    !! Epsilon to detect a near zero floating point (unitless)
  REAL(r_std), PARAMETER :: large_value = 1.E33_r_std    !! some large value (for stomate) (unitless)

  REAL(r_std), SAVE :: alpha_nudge_mc                    !! Nudging constant for soil moisture 
!$OMP THREADPRIVATE(alpha_nudge_mc)
  REAL(r_std), SAVE :: alpha_nudge_snow                  !! Nudging constant for snow variables
!$OMP THREADPRIVATE(alpha_nudge_snow)

  !
  !  DIMENSIONING AND INDICES PARAMETERS  
  !
  INTEGER(i_std), PARAMETER :: ibare_sechiba = 1 !! Index for bare soil in Sechiba (unitless)
  INTEGER(i_std), PARAMETER :: ivis = 1          !! index for albedo in visible range (unitless)
  INTEGER(i_std), PARAMETER :: inir = 2          !! index for albeod i near-infrared range (unitless) 
  INTEGER(i_std), PARAMETER :: nnobio = 1        !! Number of other surface types: land ice (lakes,cities, ...) (unitless)
  INTEGER(i_std), PARAMETER :: iice = 1          !! Index for land ice (see nnobio) (unitless)
  !-
  !! Soil
  INTEGER(i_std), PARAMETER :: classnb = 9       !! Levels of soil colour classification (unitless)
  !-
  INTEGER(i_std), PARAMETER :: nleafages = 4     !! leaf age discretisation ( 1 = no discretisation )(unitless)
  !-
  !! litter fractions: indices (unitless)
  INTEGER(i_std), PARAMETER :: ileaf = 1         !! Index for leaf compartment (unitless)
  INTEGER(i_std), PARAMETER :: isapabove = 2     !! Index for sapwood above compartment (unitless)
  INTEGER(i_std), PARAMETER :: isapbelow = 3     !! Index for sapwood below compartment (unitless)
  INTEGER(i_std), PARAMETER :: iheartabove = 4   !! Index for heartwood above compartment (unitless)
  INTEGER(i_std), PARAMETER :: iheartbelow = 5   !! Index for heartwood below compartment (unitless)
  INTEGER(i_std), PARAMETER :: iroot = 6         !! Index for roots compartment (unitless)
  INTEGER(i_std), PARAMETER :: ifruit = 7        !! Index for fruits compartment (unitless)
  INTEGER(i_std), PARAMETER :: icarbres = 8      !! Index for reserve compartment (unitless)
  INTEGER(i_std), PARAMETER :: ilabile = 9       !! Index for reserve compartment (unitless) 
  INTEGER(i_std), PARAMETER :: nparts = 9        !! Number of biomass compartments (unitless)
  !-
  !! indices for assimilation parameters 
  INTEGER(i_std), PARAMETER :: ivcmax = 1        !! Index for vcmax (assimilation parameters) (unitless)
  INTEGER(i_std), PARAMETER :: inue = 2          !! Index for nue (assimilationbn parameters) (unitless)
  INTEGER(i_std), PARAMETER :: ileafN = 3        !! Index for leaf N (assimilationbn parameters) (unitless)
  INTEGER(i_std), PARAMETER :: npco2 = 3         !! Number of assimilation parameters (unitless)
  !-
  !! trees and litter: indices for the parts of heart-
  !! and sapwood above and below the ground 
  INTEGER(i_std), PARAMETER :: iabove = 1       !! Index for above part (unitless)
  INTEGER(i_std), PARAMETER :: ibelow = 2       !! Index for below part (unitless)
  INTEGER(i_std), PARAMETER :: nlevs = 2        !! Number of levels for trees and litter (unitless)
  !-
  !! litter: indices for metabolic and structural part
  INTEGER(i_std), PARAMETER :: imetabolic = 1   !! Index for metabolic litter (unitless)
  INTEGER(i_std), PARAMETER :: istructural = 2  !! Index for structural litter (unitless)
  INTEGER(i_std), PARAMETER :: iwoody = 3       !! Index for woody litter (unitless)
  INTEGER(i_std), PARAMETER :: nlitt = 3        !! Number of levels for litter compartments (unitless)
  !-
  !! carbon pools: indices
  INTEGER(i_std), PARAMETER :: iactive = 1      !! Index for active carbon pool (unitless)
  INTEGER(i_std), PARAMETER :: islow = 2        !! Index for slow carbon pool (unitless)
  INTEGER(i_std), PARAMETER :: ipassive = 3     !! Index for passive carbon pool (unitless)
  INTEGER(i_std), PARAMETER :: isurface = 4     !! Index for passive carbon pool (unitless)
  INTEGER(i_std), PARAMETER :: ncarb = 4        !! Number of soil carbon pools (unitless)
  !-
  !! For isotopes and nitrogen
  INTEGER(i_std), PARAMETER :: nelements = 2    !! Number of isotopes considered
  INTEGER(i_std), PARAMETER :: icarbon = 1      !! Index for carbon 
  INTEGER(i_std), PARAMETER :: initrogen = 2    !! Index for nitrogen 
  !! N-cycle : indices
  INTEGER(i_std), PARAMETER :: iammonium = 1    !! Index for Ammonium 
  INTEGER(i_std), PARAMETER :: initrate  = 2    !! Index for Nitrate
  INTEGER(i_std), PARAMETER :: inox      = 3    !! Index for NOX
  INTEGER(i_std), PARAMETER :: initrous  = 4    !! Index for N2O
  INTEGER(i_std), PARAMETER :: idinitro  = 5    !! Index for N2
  INTEGER(i_std), PARAMETER :: nionspec  = 2    !! Number of ionics form considered (ammonium, nitrate)
  INTEGER(i_std), PARAMETER :: nnspec    = 5    !! Number of N-species considered

  INTEGER(i_std), PARAMETER :: iatm_ammo = 1    !! Index for N input from Ammonium N atmospheric deposition
  INTEGER(i_std), PARAMETER :: iatm_nitr = 2    !! Index for N input from Nitrate N atmospheric deposition
  INTEGER(i_std), PARAMETER :: ibnf      = 3    !! Index for N input from BNF
  INTEGER(i_std), PARAMETER :: imanure   = 4    !! Index for N input from Manure
  INTEGER(i_std), PARAMETER :: ifert_ammo = 5   !! Index for N input from Fertilisation
  INTEGER(i_std), PARAMETER :: ifert_nitr = 6   !! Index for N input from Fertilisation
  INTEGER(i_std), PARAMETER :: ninput    = 6    !! Number of N-input considered  
 

  INTEGER(i_std), PARAMETER :: i_nh4_to_no3 = 1 !! Index for NO3 production
  INTEGER(i_std), PARAMETER :: i_nh4_to_no  = 2 !! Index for NO production
  INTEGER(i_std), PARAMETER :: i_nh4_to_n2o = 3 !! Index for N2O production

  INTEGER(i_std), PARAMETER :: i_no3_to_nox = 1  !! Index for NO3 consumption
  INTEGER(i_std), PARAMETER :: i_nox_to_n2o  = 2 !! Index for NO/Nox consumption
  INTEGER(i_std), PARAMETER :: i_n2o_to_n2 = 3   !! Index for N2O consumption


  !
  !! Indices used for analytical spin-up
  INTEGER(i_std), PARAMETER :: nbpools = 10              !! Total number of carbon pools (unitless)
  INTEGER(i_std), PARAMETER :: istructural_above = 1    !! Index for structural litter above (unitless)
  INTEGER(i_std), PARAMETER :: istructural_below = 2    !! Index for structural litter below (unitless)
  INTEGER(i_std), PARAMETER :: imetabolic_above = 3     !! Index for metabolic litter above (unitless)
  INTEGER(i_std), PARAMETER :: imetabolic_below = 4     !! Index for metabolic litter below (unitless)
  INTEGER(i_std), PARAMETER :: iwoody_above = 5         !! Index for woody litter above (unitless)
  INTEGER(i_std), PARAMETER :: iwoody_below = 6         !! Index for woody litter below (unitless)
  INTEGER(i_std), PARAMETER :: iactive_pool = 7         !! Index for active carbon pool (unitless)
  INTEGER(i_std), PARAMETER :: islow_pool   = 8         !! Index for slow carbon pool (unitless)
  INTEGER(i_std), PARAMETER :: ipassive_pool = 9        !! Index for passive carbon pool (unitless)
  INTEGER(i_std), PARAMETER :: isurface_pool = 10       !! Index for surface carbon pool (unitless)
  !
  !! Indicies used for output variables on Landuse tiles defined according to LUMIP project
  !! Note that ORCHIDEE do not represent pasture and urban land. Therefor the variables will have 
  !! val_exp as missing value for these tiles. 
  INTEGER(i_std), PARAMETER :: nlut=4                   !! Total number of landuse tiles according to LUMIP
  INTEGER(i_std), PARAMETER :: id_psl=1                 !! Index for primary and secondary land
  INTEGER(i_std), PARAMETER :: id_crp=2                 !! Index for crop land
  INTEGER(i_std), PARAMETER :: id_pst=3                 !! Index for pasture land
  INTEGER(i_std), PARAMETER :: id_urb=4                 !! Index for urban land


  !
  ! NUMERICAL AND PHYSICS CONSTANTS
  !
  !

  !-
  ! 1. Mathematical and numerical constants
  !-
  REAL(r_std), PARAMETER :: pi = 3.141592653589793238   !! pi souce : http://mathworld.wolfram.com/Pi.html (unitless)
  REAL(r_std), PARAMETER :: euler = 2.71828182845904523 !! e source : http://mathworld.wolfram.com/e.html (unitless)
  REAL(r_std), PARAMETER :: zero = 0._r_std             !! Numerical constant set to 0 (unitless)
  REAL(r_std), PARAMETER :: undemi = 0.5_r_std          !! Numerical constant set to 1/2 (unitless)
  REAL(r_std), PARAMETER :: un = 1._r_std               !! Numerical constant set to 1 (unitless)
  REAL(r_std), PARAMETER :: moins_un = -1._r_std        !! Numerical constant set to -1 (unitless)
  REAL(r_std), PARAMETER :: deux = 2._r_std             !! Numerical constant set to 2 (unitless)
  REAL(r_std), PARAMETER :: trois = 3._r_std            !! Numerical constant set to 3 (unitless)
  REAL(r_std), PARAMETER :: quatre = 4._r_std           !! Numerical constant set to 4 (unitless)
  REAL(r_std), PARAMETER :: cinq = 5._r_std             !![DISPENSABLE] Numerical constant set to 5 (unitless)
  REAL(r_std), PARAMETER :: six = 6._r_std              !![DISPENSABLE] Numerical constant set to 6 (unitless)
  REAL(r_std), PARAMETER :: huit = 8._r_std             !! Numerical constant set to 8 (unitless)
  REAL(r_std), PARAMETER :: mille = 1000._r_std         !! Numerical constant set to 1000 (unitless)

  !-
  ! 2 . Physics
  !-
  REAL(r_std), PARAMETER :: R_Earth = 6378000.              !! radius of the Earth : Earth radius ~= Equatorial radius (m)
  REAL(r_std), PARAMETER :: mincos  = 0.0001                !! Minimum cosine value used for interpolation (unitless) 
  REAL(r_std), PARAMETER :: pb_std = 1013.                  !! standard pressure (hPa)
  REAL(r_std), PARAMETER :: ZeroCelsius = 273.15            !! 0 degre Celsius in degre Kelvin (K)
  REAL(r_std), PARAMETER :: tp_00 = 273.15                  !! 0 degre Celsius in degre Kelvin (K)
  REAL(r_std), PARAMETER :: chalsu0 = 2.8345E06             !! Latent heat of sublimation (J.kg^{-1})
  REAL(r_std), PARAMETER :: chalev0 = 2.5008E06             !! Latent heat of evaporation (J.kg^{-1}) 
  REAL(r_std), PARAMETER :: chalfu0 = chalsu0-chalev0       !! Latent heat of fusion (J.kg^{-1}) 
  REAL(r_std), PARAMETER :: c_stefan = 5.6697E-8            !! Stefan-Boltzman constant (W.m^{-2}.K^{-4})
  REAL(r_std), PARAMETER :: cp_air = 1004.675               !! Specific heat of dry air (J.kg^{-1}.K^{-1}) 
  REAL(r_std), PARAMETER :: cte_molr = 287.05               !! Specific constant of dry air (kg.mol^{-1}) 
  REAL(r_std), PARAMETER :: kappa = cte_molr/cp_air         !! Kappa : ratio between specific constant and specific heat 
                                                            !! of dry air (unitless)
  REAL(r_std), PARAMETER :: msmlr_air = 28.964E-03          !! Molecular weight of dry air (kg.mol^{-1})
  REAL(r_std), PARAMETER :: msmlr_h2o = 18.02E-03           !! Molecular weight of water vapor (kg.mol^{-1}) 
  REAL(r_std), PARAMETER :: cp_h2o = &                      !! Specific heat of water vapor (J.kg^{-1}.K^{-1}) 
       & cp_air*(quatre*msmlr_air)/( 3.5_r_std*msmlr_h2o) 
  REAL(r_std), PARAMETER :: cte_molr_h2o = cte_molr/quatre  !! Specific constant of water vapor (J.kg^{-1}.K^{-1}) 
  REAL(r_std), PARAMETER :: retv = msmlr_air/msmlr_h2o-un   !! Ratio between molecular weight of dry air and water 
                                                            !! vapor minus 1(unitless)  
  REAL(r_std), PARAMETER :: rvtmp2 = cp_h2o/cp_air-un       !! Ratio between specific heat of water vapor and dry air
                                                            !! minus 1 (unitless)
  REAL(r_std), PARAMETER :: cepdu2 = (0.1_r_std)**2         !! Squared wind shear (m^2.s^{-2}) 
  REAL(r_std), PARAMETER :: ct_karman = 0.41_r_std          !! Van Karmann Constant (unitless)
  REAL(r_std), PARAMETER :: cte_grav = 9.80665_r_std        !! Acceleration of the gravity (m.s^{-2})
  REAL(r_std), PARAMETER :: pa_par_hpa = 100._r_std         !! Transform pascal into hectopascal (unitless)
  REAL(r_std), PARAMETER :: RR = 8.314                      !! Ideal gas constant (J.mol^{-1}.K^{-1})
  REAL(r_std), PARAMETER :: Sct = 1370.                     !! Solar constant (W.m^{-2}) 


  INTEGER(i_std), SAVE :: testpft = 6
!$OMP THREADPRIVATE(testpft)

  !-
  ! 3. Climatic constants
  !-
  !! Constantes of the Louis scheme 
  REAL(r_std), SAVE :: cb = 5._r_std              !! Constant of the Louis scheme (unitless);
                                                  !! reference to Louis (1979)
!$OMP THREADPRIVATE(cb)
  REAL(r_std), SAVE :: cc = 5._r_std              !! Constant of the Louis scheme (unitless);
                                                  !! reference to Louis (1979)
!$OMP THREADPRIVATE(cc)
  REAL(r_std), SAVE :: cd = 5._r_std              !! Constant of the Louis scheme (unitless);
                                                  !! reference to Louis (1979)
!$OMP THREADPRIVATE(cd)
  REAL(r_std), SAVE :: rayt_cste = 125.           !! Constant in the computation of surface resistance (W.m^{-2})
!$OMP THREADPRIVATE(rayt_cste)
  REAL(r_std), SAVE :: defc_plus = 23.E-3         !! Constant in the computation of surface resistance (K.W^{-1})
!$OMP THREADPRIVATE(defc_plus)
  REAL(r_std), SAVE :: defc_mult = 1.5            !! Constant in the computation of surface resistance (K.W^{-1})
!$OMP THREADPRIVATE(defc_mult)

  !-
  ! 4. Soil thermodynamics constants
  !-
  ! Look at constantes_soil.f90


  !
  ! OPTIONAL PARTS OF THE MODEL
  !
    LOGICAL, SAVE     :: long_print = .false.       !! To set for more printing
    !$OMP THREADPRIVATE(long_print)

  LOGICAL,PARAMETER :: diag_qsat = .TRUE.         !! One of the most frequent problems is a temperature out of range
                                                  !! we provide here a way to catch that in the calling procedure. 
                                                  !! (from Jan Polcher)(true/false) 
  LOGICAL, SAVE     :: almaoutput =.FALSE.        !! Selects the type of output for the model.(true/false)
                                                  !! Value is read from run.def in intersurf_history
!$OMP THREADPRIVATE(almaoutput)

  !
  ! DIVERSE
  !
  CHARACTER(LEN=100), SAVE :: stomate_forcing_name='NONE'  !! NV080800 Name of STOMATE forcing file (unitless)
                                                           ! Compatibility with Nicolas Viovy driver.
!$OMP THREADPRIVATE(stomate_forcing_name)
  CHARACTER(LEN=100), SAVE :: stomate_Cforcing_name='NONE' !! NV080800 Name of soil forcing file (unitless)
                                                           ! Compatibility with Nicolas Viovy driver.
!$OMP THREADPRIVATE(stomate_Cforcing_name)
  CHARACTER(LEN=100), SAVE :: stomate_Cforcing_discretization_name='NONE' !! Name of soil carbon discretization forcing file (unitless)
!$OMP THREADPRIVATE(stomate_Cforcing_discretization_name)
  INTEGER(i_std), SAVE :: forcing_id                 !! Index of the forcing file (unitless)
!$OMP THREADPRIVATE(forcing_id)
  LOGICAL, SAVE :: allow_forcing_write=.TRUE.        !! Allow writing of stomate_forcing file. 
                                                     !! This variable will be set to false for teststomate.
!$OMP THREADPRIVATE(allow_forcing_write)


                         !------------------------!
                         !  SECHIBA PARAMETERS    !
                         !------------------------!
 

  !
  ! GLOBAL PARAMETERS   
  !
  REAL(r_std), SAVE :: min_wind = 0.1      !! The minimum wind (m.s^{-1})
!$OMP THREADPRIVATE(min_wind)
  REAL(r_std), PARAMETER :: min_qc = 1.e-4 !! The minimum value for qc (qc=drag*wind) used in coupled(enerbil) and forced mode (enerbil and diffuco) 
  REAL(r_std), SAVE :: snowcri = 1.5       !! Sets the amount above which only sublimation occures (kg.m^{-2})
!$OMP THREADPRIVATE(snowcri)


  !
  ! FLAGS ACTIVATING SUB-MODELS
  !
  LOGICAL, SAVE :: treat_expansion = .FALSE.   !! Do we treat PFT expansion across a grid point after introduction? (true/false)
!$OMP THREADPRIVATE(treat_expansion)
  LOGICAL, SAVE :: ok_herbivores = .FALSE.     !! flag to activate herbivores (true/false)
!$OMP THREADPRIVATE(ok_herbivores)
  LOGICAL, SAVE :: harvest_agri = .TRUE.       !! flag to harvest aboveground biomass from agricultural PFTs)(true/false)
!$OMP THREADPRIVATE(harvest_agri)
  LOGICAL, SAVE :: lpj_gap_const_mort          !! constant moratlity (true/false). Default value depend on OK_DGVM.
!$OMP THREADPRIVATE(lpj_gap_const_mort)
  LOGICAL, SAVE :: disable_fire = .TRUE.       !! flag that disable fire (true/false)
!$OMP THREADPRIVATE(disable_fire)
  LOGICAL, SAVE :: spinup_analytic = .FALSE.   !! Flag to activate analytical resolution for spinup (true/false)
!$OMP THREADPRIVATE(spinup_analytic)
  LOGICAL, SAVE :: ok_soil_carbon_discretization !! Flag to activate soil carbon discretization (vertical carbon and soil carbon thermal insulation)
!$OMP THREADPRIVATE(ok_soil_carbon_discretization)
  LOGICAL, SAVE :: ok_forcesoil_write            !! Flag to activate soil carbon discretization output file to be used for forcesoil program
!$OMP THREADPRIVATE(ok_forcesoil_write)

  !
  ! CONFIGURATION VEGETATION
  !
  LOGICAL, SAVE :: agriculture = .TRUE.    !! allow agricultural PFTs (true/false)
!$OMP THREADPRIVATE(agriculture)
  LOGICAL, SAVE :: impveg = .FALSE.        !! Impose vegetation ? (true/false)
!$OMP THREADPRIVATE(impveg)
  LOGICAL, SAVE :: impsoilt = .FALSE.      !! Impose soil ? (true/false)
!$OMP THREADPRIVATE(impsoilt)
  LOGICAL, SAVE :: impose_ninput_dep = .FALSE. !! Impose N input values from the atmosphere ? (true/false)
!$OMP THREADPRIVATE(impose_ninput_dep)
  LOGICAL, SAVE :: impose_ninput_fert = .FALSE. !! Impose N input values from fertilizer ? (true/false)        
!$OMP THREADPRIVATE(impose_ninput_fert)    
 LOGICAL, SAVE :: impose_ninput_manure = .FALSE. !! Impose N input values from manure? (true/false)        
!$OMP THREADPRIVATE(impose_ninput_manure)                     
  LOGICAL, SAVE :: impose_ninput_bnf = .FALSE. !! Impose N input values from biological nitrogen fixation (BNF) ? (true/false)        
!$OMP THREADPRIVATE(impose_ninput_bnf) 

  LOGICAL, SAVE :: do_now_stomate_lcchange = .FALSE.  !! Time to call lcchange in stomate_lpj
!$OMP THREADPRIVATE(do_now_stomate_lcchange)
  LOGICAL, SAVE :: do_now_stomate_woodharvest = .FALSE.  !! Time to call woodharvest in stomate_lpj
!$OMP THREADPRIVATE(do_now_stomate_woodharvest)
  LOGICAL, SAVE :: do_now_stomate_wood_products = .FALSE. !! Time to call wood product use and decomposition
!$OMP THREADPRIVATE(do_now_stomate_wood_products)
  LOGICAL, SAVE :: done_stomate_lcchange = .FALSE.    !! If true, call lcchange in stomate_lpj has just been done. 
!$OMP THREADPRIVATE(done_stomate_lcchange)
  LOGICAL, SAVE :: read_lai = .FALSE.      !! Flag to read a map of LAI if STOMATE is not activated (true/false)
!$OMP THREADPRIVATE(read_lai)
  LOGICAL, SAVE :: vegetmap_reset = .FALSE.!! Reset the vegetation map and reset carbon related variables
!$OMP THREADPRIVATE(vegetmap_reset)
  INTEGER(i_std) , SAVE :: veget_update    !! Update frequency in years for landuse (nb of years)
!$OMP THREADPRIVATE(veget_update)
  LOGICAL, SAVE :: ninput_reinit = .TRUE.  !! To change N INPUT file in a run. (true/false)
!$OMP THREADPRIVATE(ninput_reinit)


  !
  ! PARAMETERS USED BY BOTH HYDROLOGY MODELS
  !
  REAL(r_std), SAVE :: max_snow_age = 50._r_std !! Maximum period of snow aging (days)
!$OMP THREADPRIVATE(max_snow_age)
  REAL(r_std), SAVE :: snow_trans = 0.2_r_std   !! Transformation time constant for snow (m), reduced from the value 0.3 (04/07/2016)
!$OMP THREADPRIVATE(snow_trans)
  REAL(r_std), SAVE :: sneige                   !! Lower limit of snow amount (kg.m^{-2})
!$OMP THREADPRIVATE(sneige)
  REAL(r_std), SAVE :: maxmass_snow = 3000.     !! The maximum mass of snow (kg.m^{-2})
!$OMP THREADPRIVATE(maxmass_snow)

  !! Heat capacity
  REAL(r_std), PARAMETER :: rho_water = 1000.           !! Density of water (kg/m3)
  REAL(r_std), PARAMETER :: rho_ice = 920.              !! Density of ice (kg/m3)
  REAL(r_std), PARAMETER :: rho_soil = 2700.            !! Density of soil particles (kg/m3), value from Peters-Lidard et al. 1998

  !! Thermal conductivities
  REAL(r_std), PARAMETER :: cond_water = 0.6            !! Thermal conductivity of liquid water (W/m/K)
  REAL(r_std), PARAMETER :: cond_ice = 2.2              !! Thermal conductivity of ice (W/m/K)
  REAL(r_std), PARAMETER :: cond_solid = 2.32           !! Thermal conductivity of mineral soil particles (W/m/K)

  !! Time constant of long-term soil humidity (s) 
  REAL(r_std), PARAMETER :: lhf = 0.3336*1.E6           !! Latent heat of fusion (J/kg)

  INTEGER(i_std), PARAMETER :: nsnow=3                  !! Number of levels in the snow for explicit snow scheme   
  REAL(r_std), PARAMETER    :: XMD    = 28.9644E-3 
  REAL(r_std), PARAMETER    :: XBOLTZ      = 1.380658E-23 
  REAL(r_std), PARAMETER    :: XAVOGADRO   = 6.0221367E+23 
  REAL(r_std), PARAMETER    :: XRD    = XAVOGADRO * XBOLTZ / XMD 
  REAL(r_std), PARAMETER    :: XCPD   = 7.* XRD /2. 
  REAL(r_std), PARAMETER    :: phigeoth = 0.057 ! 0. DKtest 
  REAL(r_std), PARAMETER    :: thick_min_snow = .01 

  !! The maximum snow density and water holding characterisicts 
  REAL(r_std), SAVE         :: xrhosmax = 750.  ! (kg m-3) 
!$OMP THREADPRIVATE(xrhosmax)
  REAL(r_std), SAVE         :: xwsnowholdmax1   = 0.03  ! (-) 
!$OMP THREADPRIVATE(xwsnowholdmax1)
  REAL(r_std), SAVE         :: xwsnowholdmax2   = 0.10  ! (-) 
!$OMP THREADPRIVATE(xwsnowholdmax2)
  REAL(r_std), SAVE         :: xsnowrhohold     = 200.0 ! (kg/m3) 
!$OMP THREADPRIVATE(xsnowrhohold)
  REAL(r_std), SAVE         :: xrhosmin = 50. 
!$OMP THREADPRIVATE(xrhosmin)
  REAL(r_std), PARAMETER    :: xci = 2.106e+3 
  REAL(r_std), PARAMETER    :: xrv = 6.0221367e+23 * 1.380658e-23 /18.0153e-3 

  !! ISBA-ES Critical snow depth at which snow grid thicknesses constant 
  REAL(r_std), PARAMETER    :: xsnowcritd = 0.03  ! (m) 

  !! The threshold of snow depth used for preventing numerical problem in thermal calculations
  REAL(r_std), PARAMETER    :: snowcritd_thermal = 0.01  ! (m)  
  
  !! ISBA-ES CROCUS (Pahaut 1976): snowfall density coefficients: 
  REAL(r_std), PARAMETER       :: snowfall_a_sn = 109.0  !! (kg/m3) 
  REAL(r_std), PARAMETER       :: snowfall_b_sn =   6.0  !! (kg/m3/K) 
  REAL(r_std), PARAMETER       :: snowfall_c_sn =  26.0  !! [kg/(m7/2 s1/2)] 

  REAL(r_std), PARAMETER       :: dgrain_new_max=  2.0e-4!! (m) : Maximum grain size of new snowfall 
  
  !! Used in explicitsnow to prevent numerical problems as snow becomes vanishingly thin. 
  REAL(r_std), PARAMETER                :: psnowdzmin = .0001   ! m 
  REAL(r_std), PARAMETER                :: xsnowdmin = .000001  ! m 

  REAL(r_std), PARAMETER                :: ph2o = 1000.         !! Water density [kg/m3] 
  
  ! ISBA-ES Thermal conductivity coefficients from Anderson (1976): 
  ! see Boone, Meteo-France/CNRM Note de Centre No. 70 (2002) 
  REAL(r_std), SAVE                     :: ZSNOWTHRMCOND1 = 0.02    ! [W/m/K] 
!$OMP THREADPRIVATE(ZSNOWTHRMCOND1)
  REAL(r_std), SAVE                     :: ZSNOWTHRMCOND2 = 2.5E-6  ! [W m5/(kg2 K)] 
!$OMP THREADPRIVATE(ZSNOWTHRMCOND2)
  
  ! ISBA-ES Thermal conductivity: Implicit vapor diffn effects 
  ! (sig only for new snow OR high altitudes) 
  ! from Sun et al. (1999): based on data from Jordan (1991) 
  ! see Boone, Meteo-France/CNRM Note de Centre No. 70 (2002) 
  ! 
  REAL(r_std), SAVE                       :: ZSNOWTHRMCOND_AVAP  = -0.06023 ! (W/m/K) 
!$OMP THREADPRIVATE(ZSNOWTHRMCOND_AVAP)
  REAL(r_std), SAVE                       :: ZSNOWTHRMCOND_BVAP  = -2.5425  ! (W/m) 
!$OMP THREADPRIVATE(ZSNOWTHRMCOND_BVAP)
  REAL(r_std), SAVE                       :: ZSNOWTHRMCOND_CVAP  = -289.99  ! (K) 
!$OMP THREADPRIVATE(ZSNOWTHRMCOND_CVAP)
  REAL(r_std),SAVE :: xansmax = 0.85      !! Maxmimum snow albedo
!$OMP THREADPRIVATE(xansmax)
  REAL(r_std),SAVE :: xansmin = 0.50      !! Miniumum snow albedo
!$OMP THREADPRIVATE(xansmin)
  REAL(r_std),SAVE :: xans_todry = 0.008  !! Albedo decay rate for dry snow
!$OMP THREADPRIVATE(xans_todry)
  REAL(r_std),SAVE :: xans_t = 0.240      !! Albedo decay rate for wet snow
!$OMP THREADPRIVATE(xans_t)

  ! ISBA-ES Thermal conductivity coefficients from Anderson (1976):
  ! see Boone, Meteo-France/CNRM Note de Centre No. 70 (2002)
  REAL(r_std), PARAMETER                  :: XP00 = 1.E5

  ! ISBA-ES Thermal conductivity: Implicit vapor diffn effects
  ! (sig only for new snow OR high altitudes)
  ! from Sun et al. (1999): based on data from Jordan (1991)
  ! see Boone, Meteo-France/CNRM Note de Centre No. 70 (2002)
  !
  REAL(r_std), SAVE          :: ZSNOWCMPCT_RHOD  = 150.0        !! (kg/m3)
!$OMP THREADPRIVATE(ZSNOWCMPCT_RHOD)
  REAL(r_std), SAVE          :: ZSNOWCMPCT_ACM   = 2.8e-6       !! (1/s)
!$OMP THREADPRIVATE(ZSNOWCMPCT_ACM)
  REAL(r_std), SAVE          :: ZSNOWCMPCT_BCM   = 0.04         !! (1/K) 
!$OMP THREADPRIVATE(ZSNOWCMPCT_BCM)
  REAL(r_std), SAVE          :: ZSNOWCMPCT_CCM   = 460.         !! (m3/kg)
!$OMP THREADPRIVATE(ZSNOWCMPCT_CCM)
  REAL(r_std), SAVE          :: ZSNOWCMPCT_V0    = 3.7e7        !! (Pa/s) 
!$OMP THREADPRIVATE(ZSNOWCMPCT_V0)
  REAL(r_std), SAVE          :: ZSNOWCMPCT_VT    = 0.081        !! (1/K)
!$OMP THREADPRIVATE(ZSNOWCMPCT_VT)
  REAL(r_std), SAVE          :: ZSNOWCMPCT_VR    = 0.018        !! (m3/kg)
!$OMP THREADPRIVATE(ZSNOWCMPCT_VR)

  !
  ! BVOC : Biogenic activity  for each age class
  !
  REAL(r_std), SAVE, DIMENSION(nleafages) :: iso_activity = (/0.5, 1.5, 1.5, 0.5/)     !! Biogenic activity for each 
                                                                                       !! age class : isoprene (unitless)
!$OMP THREADPRIVATE(iso_activity)
  REAL(r_std), SAVE, DIMENSION(nleafages) :: methanol_activity = (/1., 1., 0.5, 0.5/)  !! Biogenic activity for each
                                                                                       !! age class : methanol (unnitless)
!$OMP THREADPRIVATE(methanol_activity)

  !
  ! condveg.f90
  !

  ! 1. Scalar

  ! 1.1 Flags used inside the module

  LOGICAL, SAVE :: alb_bare_model = .FALSE. !! Switch for choosing values of bare soil 
                                            !! albedo (see header of subroutine)
                                            !! (true/false)
!$OMP THREADPRIVATE(alb_bare_model)
  LOGICAL, SAVE :: alb_bg_modis = .TRUE.    !! Switch for choosing values of bare soil 
                                            !! albedo read from file
                                            !! (true/false)
!$OMP THREADPRIVATE(alb_bg_modis)
  LOGICAL, SAVE :: impaze = .FALSE.         !! Switch for choosing surface parameters
                                            !! (see header of subroutine).  
                                            !! (true/false)
!$OMP THREADPRIVATE(impaze)
  LOGICAL, SAVE :: rough_dyn = .TRUE.       !! Chooses between two methods to calculate the 
                                            !! the roughness height : static or dynamic (varying with LAI)
                                            !! (true/false)
!$OMP THREADPRIVATE(rough_dyn)

  LOGICAL, SAVE :: sla_dyn = .TRUE.        !! Chooses between two methods to calculate the 
                                            !! specific leaf area: static or dynamic (varying with LAI or biomass)
                                            !! (true/false)
!$OMP THREADPRIVATE(sla_dyn)

  LOGICAL, SAVE :: new_watstress = .FALSE.
!$OMP THREADPRIVATE(new_watstress)

  REAL(r_std), SAVE :: alpha_watstress = 1.
!$OMP THREADPRIVATE(alpha_watstress)

  ! 1.2 Others 

  REAL(r_std), SAVE :: height_displacement = 0.66        !! Factor to calculate the zero-plane displacement
                                                         !! height from vegetation height (m)
!$OMP THREADPRIVATE(height_displacement)
  REAL(r_std), SAVE :: z0_bare = 0.01                    !! bare soil roughness length (m)
!$OMP THREADPRIVATE(z0_bare)
  REAL(r_std), SAVE :: z0_ice = 0.001                    !! ice roughness length (m)
!$OMP THREADPRIVATE(z0_ice)
  REAL(r_std), SAVE :: tcst_snowa = 10.0                 !! Time constant of the albedo decay of snow (days), increased from the value 5.0 (04/07/2016)
!$OMP THREADPRIVATE(tcst_snowa)
  REAL(r_std), SAVE :: snowcri_alb = 10.                 !! Critical value for computation of snow albedo (cm)
!$OMP THREADPRIVATE(snowcri_alb)
  REAL(r_std), SAVE :: fixed_snow_albedo = undef_sechiba !! To choose a fixed snow albedo value (unitless)
!$OMP THREADPRIVATE(fixed_snow_albedo)
  REAL(r_std), SAVE :: z0_scal = 0.15                    !! Surface roughness height imposed (m)
!$OMP THREADPRIVATE(z0_scal)
  REAL(r_std), SAVE :: roughheight_scal = zero           !! Effective roughness Height depending on zero-plane 
                                                         !! displacement height (m) (imposed)
!$OMP THREADPRIVATE(roughheight_scal)
  REAL(r_std), SAVE :: emis_scal = 1.0                   !! Surface emissivity imposed (unitless)
!$OMP THREADPRIVATE(emis_scal)

  REAL(r_std), SAVE :: c1 = 0.32                         !! Constant used in the formulation of the ratio of 
!$OMP THREADPRIVATE(c1)                                  !! friction velocity to the wind speed at the canopy top
                                                         !! see Ershadi et al. (2015) for more info
  REAL(r_std), SAVE :: c2 = 0.264                        !! Constant used in the formulation of the ratio of 
!$OMP THREADPRIVATE(c2)                                  !! friction velocity to the wind speed at the canopy top
                                                         !! see Ershadi et al. (2015) for more info
  REAL(r_std), SAVE :: c3 = 15.1                         !! Constant used in the formulation of the ratio of 
!$OMP THREADPRIVATE(c3)                                  !! friction velocity to the wind speed at the canopy top
                                                         !! see Ershadi et al. (2015) for more info
  REAL(r_std), SAVE :: Cdrag_foliage = 0.2               !! Drag coefficient of the foliage
!$OMP THREADPRIVATE(Cdrag_foliage)                       !! See Ershadi et al. (2015) and Su et. al (2001) for more info
  REAL(r_std), SAVE :: Ct = 0.01                         !! Heat transfer coefficient of the leaf
!$OMP THREADPRIVATE(Ct)                                  !! See Ershadi et al. (2015) and Su et. al (2001) for more info
  REAL(r_std), SAVE :: Prandtl = 0.71                    !! Prandtl number used in the calculation of Ct_star
!$OMP THREADPRIVATE(Prandtl)                             !! See Su et. al (2001) for more info



  ! 2. Arrays

  REAL(r_std), SAVE, DIMENSION(2) :: alb_deadleaf = (/ .12, .35/)    !! albedo of dead leaves, VIS+NIR (unitless)
!$OMP THREADPRIVATE(alb_deadleaf)
  REAL(r_std), SAVE, DIMENSION(2) :: alb_ice = (/ .60, .20/)         !! albedo of ice, VIS+NIR (unitless)
!$OMP THREADPRIVATE(alb_ice)
  REAL(r_std), SAVE, DIMENSION(2) :: albedo_scal = (/ 0.25, 0.25 /)  !! Albedo values for visible and near-infrared 
                                                                     !! used imposed (unitless) 
!$OMP THREADPRIVATE(albedo_scal)
  REAL(r_std) , SAVE, DIMENSION(classnb) :: vis_dry = (/0.24,&
       &0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10, 0.27/)  !! Soil albedo values to soil colour classification:
                                                          !! dry soil albedo values in visible range
!$OMP THREADPRIVATE(vis_dry)
  REAL(r_std), SAVE, DIMENSION(classnb) :: nir_dry = (/0.48,&
       &0.44, 0.40, 0.36, 0.32, 0.28, 0.24, 0.20, 0.55/)  !! Soil albedo values to soil colour classification:
                                                          !! dry soil albedo values in near-infrared range 
!$OMP THREADPRIVATE(nir_dry)
  REAL(r_std), SAVE, DIMENSION(classnb) :: vis_wet = (/0.12,&
       &0.11, 0.10, 0.09, 0.08, 0.07, 0.06, 0.05, 0.15/)  !! Soil albedo values to soil colour classification:
                                                          !! wet soil albedo values in visible range 
!$OMP THREADPRIVATE(vis_wet)
  REAL(r_std), SAVE, DIMENSION(classnb) :: nir_wet = (/0.24,&
       &0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10, 0.31/)  !! Soil albedo values to soil colour classification:
                                                          !! wet soil albedo values in near-infrared range
!$OMP THREADPRIVATE(nir_wet)
  REAL(r_std), SAVE, DIMENSION(classnb) :: albsoil_vis = (/ &
       &0.18, 0.16, 0.16, 0.15, 0.12, 0.105, 0.09, 0.075, 0.25/)   !! Soil albedo values to soil colour classification:
                                                                   !! Averaged of wet and dry soil albedo values
                                                                   !! in visible and near-infrared range
!$OMP THREADPRIVATE(albsoil_vis) 
  REAL(r_std), SAVE, DIMENSION(classnb) :: albsoil_nir = (/ &
       &0.36, 0.34, 0.34, 0.33, 0.30, 0.25, 0.20, 0.15, 0.45/)  !! Soil albedo values to soil colour classification:
                                                                !! Averaged of wet and dry soil albedo values
                                                                !! in visible and near-infrared range
!$OMP THREADPRIVATE(albsoil_nir)

  !
  ! diffuco.f90
  !

  ! 0. Constants

  REAL(r_std), PARAMETER :: Tetens_1 = 0.622         !! Ratio between molecular weight of water vapor and molecular weight  
                                                     !! of dry air (unitless)
  REAL(r_std), PARAMETER :: Tetens_2 = 0.378         !!
  REAL(r_std), PARAMETER :: ratio_H2O_to_CO2 = 1.6   !! Ratio of water vapor diffusivity to the CO2 diffusivity (unitless)
  REAL(r_std), PARAMETER :: mol_to_m_1 = 0.0244      !!
  REAL(r_std), PARAMETER :: RG_to_PAR = 0.5          !!
  REAL(r_std), PARAMETER :: W_to_mol = 4.6          !! W_to_mmol * RG_to_PAR = 2.3

  ! 1. Scalar

  INTEGER(i_std), SAVE :: nlai = 20             !! Number of LAI levels (unitless)
!$OMP THREADPRIVATE(nlai)
  LOGICAL, SAVE :: ldq_cdrag_from_gcm = .FALSE. !! Set to .TRUE. if you want q_cdrag coming from GCM
!$OMP THREADPRIVATE(ldq_cdrag_from_gcm)
  REAL(r_std), SAVE :: laimax = 12.             !! Maximal LAI used for splitting LAI into N layers (m^2.m^{-2})
!$OMP THREADPRIVATE(laimax)
  LOGICAL, SAVE :: downregulation_co2 = .TRUE.             !! Set to .TRUE. if you want CO2 downregulation.
!$OMP THREADPRIVATE(downregulation_co2)
  REAL(r_std), SAVE :: downregulation_co2_baselevel = 380. !! CO2 base level (ppm)
!$OMP THREADPRIVATE(downregulation_co2_baselevel)

  REAL(r_std), SAVE :: gb_ref = 1./25.                     !! Leaf bulk boundary layer resistance (s m-1)
!$OMP THREADPRIVATE(gb_ref)

  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: lai_level_depth = 0.15  !!
!$OMP THREADPRIVATE(lai_level_depth)
!
  REAL(r_std), SAVE, DIMENSION(6) :: dew_veg_poly_coeff = &            !! coefficients of the 5 degree polynomomial used
  & (/ 0.887773, 0.205673, 0.110112, 0.014843,  0.000824,  0.000017 /) !! in the equation of coeff_dew_veg
!$OMP THREADPRIVATE(dew_veg_poly_coeff)
!
  REAL(r_std), SAVE               :: Oi=210000.    !! Intercellular oxygen partial pressure (ubar)
!$OMP THREADPRIVATE(Oi)
  !
  ! slowproc.f90 
  !

  ! 1. Scalar

  INTEGER(i_std), SAVE :: ninput_year_orig = 0       !!  first year for N inputs (number)
!$OMP THREADPRIVATE(ninput_year_orig)
  LOGICAL, SAVE :: ninput_suffix_year = .FALSE.      !! Do the Ninput datasets have a 'year' suffix ? (y/n)  
!$OMP THREADPRIVATE(ninput_suffix_year)
  REAL(r_std), SAVE :: bulk_default = 1000.0           !! Default value for bulk density of soil (kg/m3)
!$OMP THREADPRIVATE(bulk_default)
  REAL(r_std), SAVE :: ph_default = 5.5              !! Default value for pH of soil (-)
!$OMP THREADPRIVATE(ph_default)

! The default value for clay fraction is an heritage, with no documentation nor justification.   
  REAL(r_std), SAVE :: clayfraction_default = 0.2    !! Default value for clay fraction (0-1, unitless)
!$OMP THREADPRIVATE(clayfraction_default)
! We need to output sand and silt fractiosn for SP-MIP, and the following default values, corresponding to a Loamy soil, are selected.
  REAL(r_std), SAVE :: sandfraction_default = 0.4    !! Default value for sand fraction (0-1, unitless)
!$OMP THREADPRIVATE(sandfraction_default)
  REAL(r_std), SAVE :: siltfraction_default = 0.4    !! Default value for silt fraction (0-1, unitless)
!$OMP THREADPRIVATE(siltfraction_default)
  REAL(r_std), SAVE :: min_vegfrac = 0.001           !! Minimal fraction of mesh a vegetation type can occupy (0-1, unitless)
!$OMP THREADPRIVATE(min_vegfrac)
  REAL(r_std), SAVE :: frac_nobio_fixed_test_1 = 0.0 !! Value for frac_nobio for tests in 0-dim simulations (0-1, unitless)
!$OMP THREADPRIVATE(frac_nobio_fixed_test_1)
  
  REAL(r_std), SAVE :: stempdiag_bid = 280.          !! only needed for an initial LAI if there is no restart file
!$OMP THREADPRIVATE(stempdiag_bid)


                           !-----------------------------!
                           !  STOMATE AND LPJ PARAMETERS !
                           !-----------------------------!


  !
  ! lpj_constraints.f90
  !
  
  ! 1. Scalar

  REAL(r_std), SAVE  :: too_long = 5.      !! longest sustainable time without 
                                           !! regeneration (vernalization) (years)
!$OMP THREADPRIVATE(too_long)


  !
  ! lpj_establish.f90
  !

  ! 1. Scalar

  REAL(r_std), SAVE :: estab_max_tree = 0.12   !! Maximum tree establishment rate (ind/m2/dt_stomate)
!$OMP THREADPRIVATE(estab_max_tree)
  REAL(r_std), SAVE :: estab_max_grass = 0.12  !! Maximum grass establishment rate (ind/m2/dt_stomate)
!$OMP THREADPRIVATE(estab_max_grass)
  
  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: establish_scal_fact = 5.  !!
!$OMP THREADPRIVATE(establish_scal_fact)
  REAL(r_std), SAVE :: max_tree_coverage = 0.98  !! (0-1, unitless)
!$OMP THREADPRIVATE(max_tree_coverage)
  REAL(r_std), SAVE :: ind_0_estab = 0.2         !! = ind_0 * 10.
!$OMP THREADPRIVATE(ind_0_estab)


  !
  ! lpj_fire.f90
  !

  ! 1. Scalar

  REAL(r_std), SAVE :: tau_fire = 30.           !! Time scale for memory of the fire index (days).
!$OMP THREADPRIVATE(tau_fire)
  REAL(r_std), SAVE :: litter_crit = 200.       !! Critical litter quantity for fire
                                                !! below which iginitions extinguish 
                                                !! @tex $(gC m^{-2})$ @endtex
!$OMP THREADPRIVATE(litter_crit)
  REAL(r_std), SAVE :: fire_resist_lignin = 0.5 !!
!$OMP THREADPRIVATE(fire_resist_lignin)
  ! 2. Arrays

  REAL(r_std), SAVE, DIMENSION(nparts) :: co2frac = &    !! The fraction of the different biomass 
       & (/ .95, .95, 0., 0.3, 0., 0., .95, .95, .95 /)       !! compartments emitted to the atmosphere 
!$OMP THREADPRIVATE(co2frac)                                                         !! when burned (unitless, 0-1)  

  ! 3. Coefficients of equations

  REAL(r_std), SAVE, DIMENSION(3) :: bcfrac_coeff = (/ .3,  1.3,  88.2 /)         !! (unitless)
!$OMP THREADPRIVATE(bcfrac_coeff)
  REAL(r_std), SAVE, DIMENSION(4) :: firefrac_coeff = (/ 0.45, 0.8, 0.6, 0.13 /)  !! (unitless)
!$OMP THREADPRIVATE(firefrac_coeff)

  !
  ! lpj_gap.f90
  !

  ! 1. Scalar

  REAL(r_std), SAVE :: ref_greff = 0.035         !! Asymptotic maximum mortality rate
                                                 !! @tex $(year^{-1})$ @endtex
!$OMP THREADPRIVATE(ref_greff)

  !               
  ! lpj_light.f90 
  !              

  ! 1. Scalar
  
  LOGICAL, SAVE :: annual_increase = .TRUE. !! for diagnosis of fpc increase, compare today's fpc to last year's maximum (T) or
                                            !! to fpc of last time step (F)? (true/false)
!$OMP THREADPRIVATE(annual_increase)
  REAL(r_std), SAVE :: min_cover = 0.05     !! For trees, minimum fraction of crown area occupied
                                            !! (due to its branches etc.) (0-1, unitless)
                                            !! This means that only a small fraction of its crown area
                                            !! can be invaded by other trees.
!$OMP THREADPRIVATE(min_cover)
  !
  ! lpj_pftinout.f90 
  !

  ! 1. Scalar

  REAL(r_std), SAVE :: min_avail = 0.01         !! minimum availability
!$OMP THREADPRIVATE(min_avail)
  REAL(r_std), SAVE :: ind_0 = 0.02             !! initial density of individuals
!$OMP THREADPRIVATE(ind_0)
  ! 3. Coefficients of equations
  
  REAL(r_std), SAVE :: RIP_time_min = 1.25      !! test whether the PFT has been eliminated lately (years)
!$OMP THREADPRIVATE(RIP_time_min)
  REAL(r_std), SAVE :: npp_longterm_init = 10.  !! Initialisation value for npp_longterm (gC.m^{-2}.year^{-1})
!$OMP THREADPRIVATE(npp_longterm_init)
  REAL(r_std), SAVE :: everywhere_init = 0.05   !!
!$OMP THREADPRIVATE(everywhere_init)

REAL(r_std), SAVE :: z_nitrogen = 0.2            !! Curvature of the root profile (m)
!$OMP THREADPRIVATE(z_nitrogen)

!  ! 3. Coefficients of equations

!  REAL(r_std), SAVE :: lai_max_to_happy = 0.5      !!
!$OMP THREADPRIVATE(lai_max_to_happy)




  !
  ! stomate_data.f90 
  !

  ! 1. Scalar 

  ! 1.2 climatic parameters 

  REAL(r_std), SAVE :: precip_crit = 100.        !! minimum precip, in (mm/year)
!$OMP THREADPRIVATE(precip_crit)
  REAL(r_std), SAVE :: gdd_crit_estab = 150.     !! minimum gdd for establishment of saplings
!$OMP THREADPRIVATE(gdd_crit_estab)
  REAL(r_std), SAVE :: fpc_crit = 0.95           !! critical fpc, needed for light competition and establishment (0-1, unitless)
!$OMP THREADPRIVATE(fpc_crit)

  ! 1.3 sapling characteristics

  REAL(r_std), SAVE :: alpha_grass = 0.5         !! alpha coefficient for grasses (unitless)
!$OMP THREADPRIVATE(alpha_grass)
  REAL(r_std), SAVE :: alpha_tree = 1.           !! alpha coefficient for trees (unitless)
!$OMP THREADPRIVATE(alpha_tree)
  REAL(r_std), SAVE :: struct_to_leaves = 0.05  !! Fraction of structural carbon in grass and crops as a share of the leaf
                                                !! carbon pool. Only used for grasses and crops (thus NOT for trees)
	                                        !! (unitless)
!$OMP THREADPRIVATE(struct_to_leaves)

  REAL(r_std), SAVE :: labile_to_total = 0.01   !! Fraction of the labile pool in trees, grasses and crops as a share of the
                                                !! total carbon pool (accounting for the N-content of the different tissues).
                                                !! (unitless)
!$OMP THREADPRIVATE(labile_to_total)



  ! 1.4  time scales for phenology and other processes (in days)
  REAL(r_std), SAVE :: tau_hum_month = 20.        !! (days)       
!$OMP THREADPRIVATE(tau_hum_month)
  REAL(r_std), SAVE :: tau_hum_week = 7.          !! (days)  
!$OMP THREADPRIVATE(tau_hum_week)
  REAL(r_std), SAVE :: tau_t2m_month = 20.        !! (days)      
!$OMP THREADPRIVATE(tau_t2m_month)
  REAL(r_std), SAVE :: tau_t2m_week = 7.          !! (days)  
!$OMP THREADPRIVATE(tau_t2m_week)
  REAL(r_std), SAVE :: tau_tsoil_month = 20.      !! (days)     
!$OMP THREADPRIVATE(tau_tsoil_month)
  REAL(r_std), SAVE :: tau_soilhum_month = 20.    !! (days)     
!$OMP THREADPRIVATE(tau_soilhum_month)
  REAL(r_std), SAVE :: tau_gpp_week = 7.          !! (days)  
!$OMP THREADPRIVATE(tau_gpp_week)
  REAL(r_std), SAVE :: tau_sugarload_week = 7.    !! (days)
!$OMP THREADPRIVATE(tau_sugarload_week)
  REAL(r_std), SAVE :: tau_gdd = 40.              !! (days)  
!$OMP THREADPRIVATE(tau_gdd)
  REAL(r_std), SAVE :: tau_ngd = 50.              !! (days)  
!$OMP THREADPRIVATE(tau_ngd)
  REAL(r_std), SAVE :: coeff_tau_longterm = 3.    !! (unitless)
!$OMP THREADPRIVATE(coeff_tau_longterm)
  REAL(r_std), SAVE :: tau_longterm_max           !! (days)  
!$OMP THREADPRIVATE(tau_longterm_max)

  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: bm_sapl_carbres = 5.             !!
!$OMP THREADPRIVATE(bm_sapl_carbres)
  REAL(r_std), SAVE :: bm_sapl_labile = 5.             !!
!$OMP THREADPRIVATE(bm_sapl_labile)
  REAL(r_std), SAVE :: bm_sapl_sapabove = 0.5           !!
!$OMP THREADPRIVATE(bm_sapl_sapabove)
  REAL(r_std), SAVE :: bm_sapl_heartabove = 2.          !!
!$OMP THREADPRIVATE(bm_sapl_heartabove)
  REAL(r_std), SAVE :: bm_sapl_heartbelow = 2.          !!
!$OMP THREADPRIVATE(bm_sapl_heartbelow)
  REAL(r_std), SAVE :: init_sapl_mass_leaf_nat = 0.1    !!
!$OMP THREADPRIVATE(init_sapl_mass_leaf_nat)
  REAL(r_std), SAVE :: init_sapl_mass_leaf_agri = 1.    !!
!$OMP THREADPRIVATE(init_sapl_mass_leaf_agri)
  REAL(r_std), SAVE :: init_sapl_mass_carbres = 5.      !!
!$OMP THREADPRIVATE(init_sapl_mass_carbres)
  REAL(r_std), SAVE :: init_sapl_mass_labile = 5.      !!
!$OMP THREADPRIVATE(init_sapl_mass_labile)
  REAL(r_std), SAVE :: init_sapl_mass_root = 0.1        !!
!$OMP THREADPRIVATE(init_sapl_mass_root)
  REAL(r_std), SAVE :: init_sapl_mass_fruit = 0.3       !!  
!$OMP THREADPRIVATE(init_sapl_mass_fruit)
  REAL(r_std), SAVE :: cn_sapl_init = 0.5               !!
!$OMP THREADPRIVATE(cn_sapl_init)
  REAL(r_std), SAVE :: migrate_tree = 10.*1.E3          !!
!$OMP THREADPRIVATE(migrate_tree)
  REAL(r_std), SAVE :: migrate_grass = 10.*1.E3         !!
!$OMP THREADPRIVATE(migrate_grass)
  REAL(r_std), SAVE :: lai_initmin_tree = 0.3           !!
!$OMP THREADPRIVATE(lai_initmin_tree)
  REAL(r_std), SAVE :: lai_initmin_grass = 0.1          !!
!$OMP THREADPRIVATE(lai_initmin_grass)
  REAL(r_std), SAVE, DIMENSION(2) :: dia_coeff = (/ 4., 0.5 /)            !!
!$OMP THREADPRIVATE(dia_coeff)
  REAL(r_std), SAVE, DIMENSION(2) :: maxdia_coeff =(/ 100., 0.01/)        !!
!$OMP THREADPRIVATE(maxdia_coeff)
  REAL(r_std), SAVE, DIMENSION(4) :: bm_sapl_leaf = (/ 4., 4., 0.8, 5./)  !!
!$OMP THREADPRIVATE(bm_sapl_leaf)


  !
  ! stomate_litter.f90 
  !

  ! 0. Constants

  REAL(r_std), PARAMETER :: Q10 = 10.               !!

  ! 1. Scalar

  REAL(r_std), SAVE :: z_decomp = 0.2               !!  Maximum depth for soil decomposer's activity (m)
!$OMP THREADPRIVATE(z_decomp)

  ! 2. Arrays

  REAL(r_std), SAVE :: frac_soil_struct_sua = 0.4    !! corresponding to frac_soil(istructural,isurface,iabove) 
!$OMP THREADPRIVATE(frac_soil_struct_sua)
  REAL(r_std), SAVE :: frac_soil_struct_ab = 0.45   !! corresponding to frac_soil(istructural,iactive,ibelow)
!$OMP THREADPRIVATE(frac_soil_struct_ab)
  REAL(r_std), SAVE :: frac_soil_struct_sa = 0.7    !! corresponding to frac_soil(istructural,islow,iabove)
!$OMP THREADPRIVATE(frac_soil_struct_sa)
  REAL(r_std), SAVE :: frac_soil_struct_sb = 0.7    !! corresponding to frac_soil(istructural,islow,ibelow)
!$OMP THREADPRIVATE(frac_soil_struct_sb)
  REAL(r_std), SAVE :: frac_soil_metab_sua = 0.4    !! corresponding to frac_soil(imetabolic,iactive,iabove)
!$OMP THREADPRIVATE(frac_soil_metab_sua)
  REAL(r_std), SAVE :: frac_soil_metab_ab = 0.45    !! corresponding to frac_soil(imetabolic,iactive,ibelow)
!$OMP THREADPRIVATE(frac_soil_metab_ab)
  REAL(r_std), SAVE, DIMENSION(nparts) :: CN_fix = & !! C/N ratio of each plant pool (0-100, unitless)
       & (/ 40., 40., 40., 40., 40., 40., 40., 40., 40. /) 
!$OMP THREADPRIVATE(CN_fix)

  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: metabolic_ref_frac = 0.85    !! used by litter and soilcarbon (0-1, unitless)
!$OMP THREADPRIVATE(metabolic_ref_frac)
  REAL(r_std), SAVE :: metabolic_LN_ratio = 0.018   !! (0-1, unitless)   
!$OMP THREADPRIVATE(metabolic_LN_ratio)
  ! Turnover rate (yr-1) - From Parton et al., 1993
  REAL(r_std), SAVE :: turn_metabolic = 15           !!
!$OMP THREADPRIVATE(turn_metabolic)
  REAL(r_std), SAVE :: turn_struct = 4                !!
!$OMP THREADPRIVATE(turn_struct)
  REAL(r_std), SAVE :: turn_woody = 1.33              !! from DOFOCO
!$OMP THREADPRIVATE(turn_woody)
  REAL(r_std), SAVE :: soil_Q10 = 0.69              !!= ln 2
!$OMP THREADPRIVATE(soil_Q10)
  REAL(r_std), SAVE :: soil_Q10_uptake = 0.45         
!$OMP THREADPRIVATE(soil_Q10_uptake)
  REAL(r_std), SAVE :: tsoil_ref = 30.              !!
!$OMP THREADPRIVATE(tsoil_ref)
  REAL(r_std), SAVE :: litter_struct_coef = 3.      !! 
!$OMP THREADPRIVATE(litter_struct_coef)
  REAL(r_std), SAVE, DIMENSION(3) :: moist_coeff = (/ 1.1,  2.4,  0.29 /) !!
!$OMP THREADPRIVATE(moist_coeff)
  REAL(r_std), SAVE :: moistcont_min = 0.25  !! minimum soil wetness to limit the heterotrophic respiration
!$OMP THREADPRIVATE(moistcont_min)


  !
  ! stomate_lpj.f90
  !

  ! 1. Scalar

  REAL(r_std), SAVE :: frac_turnover_daily = 0.55  !! (0-1, unitless)
!$OMP THREADPRIVATE(frac_turnover_daily)


  !
  ! stomate_npp.f90 
  !

  ! 1. Scalar

  REAL(r_std), SAVE :: tax_max = 0.8 !! Maximum fraction of allocatable biomass used 
                                     !! for maintenance respiration (0-1, unitless)
!$OMP THREADPRIVATE(tax_max)


  !
  ! stomate_phenology.f90
  !

  ! 1. Scalar

  REAL(r_std), SAVE :: min_growthinit_time = 300.  !! minimum time since last beginning of a growing season (days)
!$OMP THREADPRIVATE(min_growthinit_time)
  REAL(r_std), SAVE :: moiavail_always_tree = 0.8  !! moisture monthly availability above which moisture tendency doesn't matter
                                                   !!  - for trees (0-1, unitless)
!$OMP THREADPRIVATE(moiavail_always_tree)
  REAL(r_std), SAVE :: moiavail_always_grass = 0.4 !! moisture monthly availability above which moisture tendency doesn't matter
                                                   !! - for grass (0-1, unitless)
!$OMP THREADPRIVATE(moiavail_always_grass)
  REAL(r_std), SAVE :: t_always                    !! monthly temp. above which temp. tendency doesn't matter
!$OMP THREADPRIVATE(t_always)
  REAL(r_std), SAVE :: t_always_add = 10.          !! monthly temp. above which temp. tendency doesn't matter (C)
!$OMP THREADPRIVATE(t_always_add)

  ! 3. Coefficients of equations
  
  REAL(r_std), SAVE :: gddncd_ref = 603.           !!
!$OMP THREADPRIVATE(gddncd_ref)
  REAL(r_std), SAVE :: gddncd_curve = 0.0091       !!
!$OMP THREADPRIVATE(gddncd_curve)
  REAL(r_std), SAVE :: gddncd_offset = 64.         !!
!$OMP THREADPRIVATE(gddncd_offset)


  !
  ! stomate_prescribe.f90
  !

  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: bm_sapl_rescale = 40.       !!
!$OMP THREADPRIVATE(bm_sapl_rescale)


  !
  ! stomate_resp.f90
  !

  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: maint_resp_min_vmax = 0.3   !!
!$OMP THREADPRIVATE(maint_resp_min_vmax)
  REAL(r_std), SAVE :: maint_resp_coeff = 1.4      !!
!$OMP THREADPRIVATE(maint_resp_coeff)


  !
  ! stomate_som_dynamics.f90 (in stomate_soilcarbon.f90 or
  ! stomate_soil_carbon_discretization.f90)   
  !

  ! 2. Arrays 

  ! 2.1 Fixed fraction from one pool to another (or to CO2 emission)

  REAL(r_std), SAVE :: active_to_pass_ref_frac = 0.003  !! from active pool: depends on clay content  (0-1, unitless)
                                                        !! corresponding to frac_carb(:,iactive,ipassive)
!$OMP THREADPRIVATE(active_to_pass_ref_frac)
  REAL(r_std), SAVE :: surf_to_slow_ref_frac = 0.4      !! from surface pool
                                                        !! corresponding to frac_carb(:,isurf,islow)
!$OMP THREADPRIVATE(surf_to_slow_ref_frac)

  REAL(r_std), SAVE :: active_to_CO2_ref_frac  = 0.85   !! from active pool: depends on clay content  (0-1, unitless)
                                                        !! corresponding to frac_resp(:,iactive)
!$OMP THREADPRIVATE(active_to_CO2_ref_frac)
  REAL(r_std), SAVE :: slow_to_pass_ref_frac   = 0.03  !! from slow pool: depends on clay content  (0-1, unitless) 
                                                        !! corresponding to frac_carb(:,islow,ipassive)
!$OMP THREADPRIVATE(slow_to_pass_ref_frac)
  REAL(r_std), SAVE :: slow_to_CO2_ref_frac    = 0.55   !! from slow pool (0-1, unitless) 
                                                        !! corresponding to frac_resp(:,islow)
!$OMP THREADPRIVATE(slow_to_CO2_ref_frac)
  REAL(r_std), SAVE :: pass_to_active_ref_frac = 0.45   !! from passive pool (0-1, unitless)
                                                        !! corresponding to frac_carb(:,ipassive,iactive)
!$OMP THREADPRIVATE(pass_to_active_ref_frac)
  REAL(r_std), SAVE :: pass_to_slow_ref_frac   = 0.0    !! from passive pool (0-1, unitless)
                                                        !! corresponding to frac_carb(:,ipassive,islow)
!$OMP THREADPRIVATE(pass_to_slow_ref_frac)

  ! 3. Define Variable fraction from one pool to another (function of silt and clay fraction)
  REAL(r_std), SAVE :: active_to_pass_clay_frac     = 0.032
!$OMP THREADPRIVATE(active_to_pass_clay_frac)
  REAL(r_std), SAVE :: active_to_CO2_clay_silt_frac = 0.68
!$OMP THREADPRIVATE(active_to_CO2_clay_silt_frac)
  REAL(r_std), SAVE :: slow_to_pass_clay_frac   = 0 
!$OMP THREADPRIVATE(slow_to_pass_clay_frac)

  ! C to N target ratios of differnt pools
  REAL(r_std), SAVE ::  CN_target_iactive_ref  = 15. !! CN target ratio of active pool for soil min N = 0
!$OMP THREADPRIVATE(CN_target_iactive_ref)
  REAL(r_std), SAVE ::  CN_target_islow_ref    = 20. !! CN target ratio of slow pool for soil min N = 0
!$OMP THREADPRIVATE(CN_target_islow_ref)
  REAL(r_std), SAVE ::  CN_target_ipassive_ref = 10. !! CN target ratio of passive pool for soil min N = 0
!$OMP THREADPRIVATE(CN_target_ipassive_ref)
  REAL(r_std), SAVE ::  CN_target_isurface_ref = 20. !! CN target ratio of surface pool for litter nitrogen content = 0
!$OMP THREADPRIVATE(CN_target_isurface_ref)
  REAL(r_std), SAVE ::  CN_target_iactive_Nmin  = -6.  !! CN target ratio change per mineral N unit (g m-2) for active pool
!$OMP THREADPRIVATE(CN_target_iactive_Nmin) 
  REAL(r_std), SAVE ::  CN_target_islow_Nmin    = -4.  !! CN target ratio change per mineral N unit (g m-2) for slow pool 
!$OMP THREADPRIVATE(CN_target_islow_Nmin)
  REAL(r_std), SAVE ::  CN_target_ipassive_Nmin = -1.5 !! CN target ratio change per mineral N unit (g m-2) for passive pool 
!$OMP THREADPRIVATE(CN_target_ipassive_Nmin)
  REAL(r_std), SAVE ::  CN_target_isurface_pnc  = -5.  !! CN target ratio change per plant nitrogen content unit (%) for surface pool
  !! Turnover in SOM pools (year-1)
!$OMP THREADPRIVATE(CN_target_isurface_pnc)
  REAL(r_std), SAVE :: som_turn_isurface = 6.0           !! turnover of surface pool (year-1)
!$OMP THREADPRIVATE(som_turn_isurface)
  REAL(r_std), SAVE :: som_turn_iactive  = 7.3           !! turnover of active pool (year-1)
!$OMP THREADPRIVATE(som_turn_iactive)
  REAL(r_std), SAVE :: som_turn_islow    = 0.1           !! turnover of slow pool (year-1)
!$OMP THREADPRIVATE(som_turn_islow)
  REAL(r_std), SAVE :: som_turn_ipassive = 0.00225        !! turnover of passive pool (year-1)
!$OMP THREADPRIVATE(som_turn_ipassive)
  REAL(r_std), SAVE :: fslow = 37.0                      !! convertiing factor to go from active pool turnover to slow pool turnover from Guimberteau et al 2018 GMD
!$OMP THREADPRIVATE(fslow)
  REAL(r_std), SAVE :: fpassive = 1617.45                !! convertiing factor to go from active pool turnover to passive pool turnover from Guimberteau et al 2018 GMD
!$OMP THREADPRIVATE(fpassive)
  REAL(r_std), SAVE :: stomate_tau = 4.699E6 
!$OMP THREADPRIVATE(stomate_tau)
  REAL(r_std), SAVE :: depth_modifier = 1.E6             !! e-folding depth of turnover rates,following Koven et al.,2013,Biogeosciences. A very large value means no depth modification
!$OMP THREADPRIVATE(depth_modifier)
  REAL(r_std), SAVE :: som_turn_iactive_clay_frac = 0.75 !! clay-dependant parameter impacting on turnover rate of active pool 
                                                         !! Tm parameter of Parton et al. 1993 (-)
!$OMP THREADPRIVATE(som_turn_iactive_clay_frac)

  !
  ! stomate_turnover.f90
  !

  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: new_turnover_time_ref = 20. !!(days)
!$OMP THREADPRIVATE(new_turnover_time_ref)
  REAL(r_std), SAVE :: leaf_age_crit_tref = 20.    !! (C)
!$OMP THREADPRIVATE(leaf_age_crit_tref)
  REAL(r_std), SAVE, DIMENSION(3) :: leaf_age_crit_coeff = (/ 1.5, 0.75, 10./) !! (unitless)
!$OMP THREADPRIVATE(leaf_age_crit_coeff)


  !
  ! stomate_vmax.f90
  !
 
  ! 1. Scalar

  REAL(r_std), SAVE :: vmax_offset = 0.3        !! minimum leaf efficiency (unitless)
!$OMP THREADPRIVATE(vmax_offset)
  REAL(r_std), SAVE :: leafage_firstmax = 0.03  !! relative leaf age at which efficiency
                                                !! reaches 1 (unitless)
!$OMP THREADPRIVATE(leafage_firstmax)
  REAL(r_std), SAVE :: leafage_lastmax = 0.5    !! relative leaf age at which efficiency
                                                !! falls below 1 (unitless)
!$OMP THREADPRIVATE(leafage_lastmax)
  REAL(r_std), SAVE :: leafage_old = 1.         !! relative leaf age at which efficiency
                                                !! reaches its minimum (vmax_offset) 
                                                !! (unitless)
!$OMP THREADPRIVATE(leafage_old)
  REAL(r_std), SAVE :: sugar_load_min = 0.0    !! Lower bound for sugar loading when used to regulate NUE
!$OMP THREADPRIVATE(sugar_load_min)
  REAL(r_std), SAVE :: sugar_load_max = 1.0    !! Upper bound for sugar loading when used to regulate NUE
!$OMP THREADPRIVATE(sugar_load_max)

  !
  ! nitrogen_dynamics (in stomate_soilcarbon.f90) 
  !

  ! 0. Constants
  REAL(r_std), PARAMETER :: D_air = 1.73664     !! Oxygen diffusion rate in the air = 0.07236 m2/h
                                               !! from Table 2 of Li et al, 2000
                                               !! (m**2/day)

  REAL(r_std), PARAMETER :: C_molar_mass = 12  !! Carbon Molar mass (gC mol-1) 

  REAL(r_std), PARAMETER :: Pa_to_hPa    = 0.01      !! Conversion factor from Pa to hPa (-)
  REAL(r_std), PARAMETER :: V_O2         = 0.209476  !! Volumetric fraction of O2 in air (-)

  REAL(r_std), PARAMETER :: pk_NH4 = 9.25      !! The negative logarithm of the acid dissociation constant K_NH4     
                                               !! See Table 4 of Li et al. 1992 and Appendix A of Zhang et al. 2002    

                                     
  ! 1. Scalar

  ! Coefficients for defining maximum porosity
  ! From Saxton, K.E., Rawls, W.J., Romberger, J.S., Papendick, R.I., 1986
  ! Estimationg generalized soil-water characteristics from texture. 
  ! Soil Sci. Soc. Am. J. 50, 1031-1036
  ! Cited in Table 5 (page 444) of
  ! Y. Pachepsky, W.J. Rawls
  ! Development of Pedotransfer Functions in Soil Hydrology
  ! Elsevier, 23 nov. 2004 - 542 pages
  ! http://books.google.fr/books?id=ar_lPXaJ8QkC&printsec=frontcover&hl=fr#v=onepage&q&f=false
  REAL(r_std), SAVE :: h_saxton = 0.332          !! h coefficient 
!$OMP THREADPRIVATE(h_saxton)
  REAL(r_std), SAVE :: j_saxton = -7.251*1e-4    !! j coefficient 
!$OMP THREADPRIVATE(j_saxton)
  REAL(r_std), SAVE :: k_saxton = 0.1276         !! k coefficient
!$OMP THREADPRIVATE(k_saxton)

  ! Values of the power used in the equation defining the diffusion of oxygen in soil
  ! from Table 2 of Li et al, 2000
  REAL(r_std), SAVE :: diffusionO2_power_1 = 3.33 !! (unitless)
!$OMP THREADPRIVATE(diffusionO2_power_1)
  REAL(r_std), SAVE :: diffusionO2_power_2 = 2.0  !! (unitless) 
!$OMP THREADPRIVATE(diffusionO2_power_2)

  ! Temperature-related Factors impacting on Oxygen diffusion rate
  ! From eq. 2 of Table 2 (Li et al, 2000)
  REAL(r_std), SAVE ::   F_nofrost = 1.2          !! (unitless)
!$OMP THREADPRIVATE(F_nofrost)
  REAL(r_std), SAVE ::   F_frost   = 0.8          !! (unitless)
!$OMP THREADPRIVATE(F_frost)

  ! Coefficients used in the calculation of Volumetric fraction of anaerobic microsites 
  ! a and b constants are not specified in Li et al., 2000
  ! S. Zaehle used a=0.85 and b=1 without mention to any publication
  REAL(r_std), SAVE ::   a_anvf    = 0.85 !! (-)
!$OMP THREADPRIVATE(a_anvf)
  REAL(r_std), SAVE ::   b_anvf    = 1.   !! (-)
!$OMP THREADPRIVATE(b_anvf)

  ! Coefficients used in the calculation of the Fraction of adsorbed NH4+
  ! Li et al. 1992, JGR, Table 4
  REAL(r_std), SAVE ::   a_FixNH4 = 0.41  !! (-)
!$OMP THREADPRIVATE(a_FixNH4)
  REAL(r_std), SAVE ::   b_FixNH4 = -0.47 !! (-)
!$OMP THREADPRIVATE(b_FixNH4)
  REAL(r_std), SAVE ::   clay_max = 0.63  !! (-)
!$OMP THREADPRIVATE(clay_max)

  ! Coefficients used in the calculation of the Response of Nitrification
  ! to soil moisture
  ! Zhang et al. 2002, Ecological Modelling, appendix A, page 101
  REAL(r_std), SAVE ::   fwnit_0 =  -0.0243  !! (-)
!$OMP THREADPRIVATE(fwnit_0)
  REAL(r_std), SAVE ::   fwnit_1 =   0.9975  !! (-)
!$OMP THREADPRIVATE(fwnit_1)
  REAL(r_std), SAVE ::   fwnit_2 =  -5.5368  !! (-)
!$OMP THREADPRIVATE(fwnit_2)
  REAL(r_std), SAVE ::   fwnit_3 =  17.651   !! (-)
!$OMP THREADPRIVATE(fwnit_3)
  REAL(r_std), SAVE ::   fwnit_4 = -12.904   !! (-)
!$OMP THREADPRIVATE(fwnit_4)

  ! Coefficients used in the calculation of the Response of Nitrification
  ! to Temperature
  ! Zhang et al. 2002, Ecological Modelling, appendix A, page 101
  REAL(r_std), SAVE ::   ft_nit_0 =  -0.0233 !! (-)
!$OMP THREADPRIVATE(ft_nit_0)
  REAL(r_std), SAVE ::   ft_nit_1 =   0.3094 !! (-)
!$OMP THREADPRIVATE(ft_nit_1)
  REAL(r_std), SAVE ::   ft_nit_2 =  -0.2234 !! (-)
!$OMP THREADPRIVATE(ft_nit_2)
  REAL(r_std), SAVE ::   ft_nit_3 =   0.1566 !! (-)
!$OMP THREADPRIVATE(ft_nit_3)
  REAL(r_std), SAVE ::   ft_nit_4 =  -0.0272 !! (-)
!$OMP THREADPRIVATE(ft_nit_4)

  ! Coefficients used in the calculation of the Response of Nitrification
  ! to pH
  ! Zhang et al. 2002, Ecological Modelling, appendix A, page 101
  REAL(r_std), SAVE ::   fph_0 = -1.2314  !! (-)
!$OMP THREADPRIVATE(fph_0)
  REAL(r_std), SAVE ::   fph_1 = 0.7347   !! (-)
!$OMP THREADPRIVATE(fph_1)
  REAL(r_std), SAVE ::   fph_2 = -0.0604  !! (-)
!$OMP THREADPRIVATE(fph_2)

  ! Coefficients used in the calculation of the response of NO2 or NO 
  ! production during nitrificationof to Temperature
  ! Zhang et al. 2002, Ecological Modelling, appendix A, page 102
  REAL(r_std), SAVE ::   ftv_0 = 2.72   !! (-)
!$OMP THREADPRIVATE(ftv_0)
  REAL(r_std), SAVE ::   ftv_1 = 34.6   !! (-)
!$OMP THREADPRIVATE(ftv_1)
  REAL(r_std), SAVE ::   ftv_2 = 9615.  !! (-) 
!$OMP THREADPRIVATE(ftv_2)

  REAL(r_std), SAVE ::   k_nitrif = 2.0         !! Nitrification rate at 20 ◦C and field capacity (day-1)
                                                !! Schmid et al., 2001 give a value of 0.2 per day.
                                                !! (https://doi.org/10.1023/A:1012694218748)
                                                !! OCN used a value of 2.0.  We keep the OCN value.
!$OMP THREADPRIVATE(k_nitrif)

  REAL(r_std), SAVE ::   n2o_nitrif_p = 0.0006  !! Reference n2o production per N-NO3 produced g N-N2O  (g N-NO3)-1
                                                !! From Zhang et al., 2002 - Appendix A p. 102
!$OMP THREADPRIVATE(n2o_nitrif_p)
  REAL(r_std), SAVE ::   no_nitrif_p = 0.0025   !! Reference NO production per N-NO3 produced g N-NO  (g N-NO3)-1
                                                !! From Zhang et al., 2002 - Appendix A p. 102
!$OMP THREADPRIVATE(no_nitrif_p)

  ! NO production from chemodenitrification
  ! based on Kesik et al., 2005, Biogeosciences
  ! Coefficients used in the calculation of the Response to Temperature
  REAL(r_std), SAVE ::   chemo_t0  = -31494. !! (-)
!$OMP THREADPRIVATE(chemo_t0)
  ! Coefficients use in the calculation of the Response to pH
  REAL(r_std), SAVE ::   chemo_ph0 = -1.62   !! (-)
!$OMP THREADPRIVATE(chemo_ph0)
  ! Coefficients used in the calculation of NO production from chemodenitrification
  REAL(r_std), SAVE ::   chemo_0   = 30.     !! (-)
!$OMP THREADPRIVATE(chemo_0)
  REAL(r_std), SAVE ::   chemo_1   = 16565.  !! (-)
!$OMP THREADPRIVATE(chemo_1)

  ! Denitrification processes
  ! Li et al, 2000, JGR Table 4 eq 1, 2 and 4
  !
  ! Coefficients used in the Temperature response of 
  ! relative growth rate of total denitrifiers - Eq. 2 Table 4 of Li et al., 2000
  REAL(r_std), SAVE ::   ft_denit_0 = 2.     !! (-)
!$OMP THREADPRIVATE(ft_denit_0)
  REAL(r_std), SAVE ::   ft_denit_1 = 22.5   !! (-)
!$OMP THREADPRIVATE(ft_denit_1)
  REAL(r_std), SAVE ::   ft_denit_2 = 10.    !! (-)
!$OMP THREADPRIVATE(ft_denit_2)
  !
  ! Coefficients used in the pH response of 
  ! relative growth rate of total denitrifiers - Eq. 2 Table 4 of Li et al., 2000
  REAL(r_std), SAVE ::   fph_no3_0  = 4.25    !! (-) 
!$OMP THREADPRIVATE(fph_no3_0)
  REAL(r_std), SAVE ::   fph_no3_1  = 0.5     !! (-)
!$OMP THREADPRIVATE(fph_no3_1)
  REAL(r_std), SAVE ::   fph_no_0  = 5.25     !! (-)
!$OMP THREADPRIVATE(fph_no_0)
  REAL(r_std), SAVE ::   fph_no_1  = 1.       !! (-)
!$OMP THREADPRIVATE(fph_no_1)
  REAL(r_std), SAVE ::   fph_n2o_0  = 6.25    !! (-)
!$OMP THREADPRIVATE(fph_n2o_0)
  REAL(r_std), SAVE ::   fph_n2o_1  = 1.5     !! (-)
!$OMP THREADPRIVATE(fph_n2o_1)

  REAL(r_std), SAVE ::   Kn = 0.083           !! Half Saturation of N oxydes (kgN/m3)
                                              !! Table 4 of Li et al., 2000
!$OMP THREADPRIVATE(Kn)

  
  REAL(r_std), SAVE ::   cte_bact = 0.00015
!$OMP THREADPRIVATE(cte_bact)
  
  ! Maximum Relative growth rate of Nox denitrifiers
  ! Eq.1 Table 4 Li et al., 2000 
  REAL(r_std), SAVE ::   mu_no3_max = 0.67   !! (hour-1)
!$OMP THREADPRIVATE(mu_no3_max)
  REAL(r_std), SAVE ::   mu_no_max  = 0.67   !! (hour-1)
!$OMP THREADPRIVATE(mu_no_max)
  REAL(r_std), SAVE ::   mu_n2o_max = 0.67   !! (hour-1)
!$OMP THREADPRIVATE(mu_n2o_max)

  ! Maximum growth yield of NOx denitrifiers on N oxydes
  ! Table 4 Li et al., 2000
  REAL(r_std), SAVE ::   Y_no3 = 0.401 !! (kgC / kgN)
!$OMP THREADPRIVATE(Y_no3)
  REAL(r_std), SAVE ::   Y_no  = 0.428 !! (kgC / kgN)
!$OMP THREADPRIVATE(Y_no)
  REAL(r_std), SAVE ::   Y_n2o = 0.151 !! (kgC / kgN)
!$OMP THREADPRIVATE(Y_n2o)

  ! Maintenance coefficient on N oxyde
  ! Table 4 Li et al., 2000
  REAL(r_std), SAVE ::   M_no3 = 0.09   !! (kgN / kgC / hour)
!$OMP THREADPRIVATE(M_no3)
  REAL(r_std), SAVE ::   M_no  = 0.035  !! (kgN / kgC / hour)
!$OMP THREADPRIVATE(M_no)
  REAL(r_std), SAVE ::   M_n2o = 0.079  !! (kgN / kgC / hour)
!$OMP THREADPRIVATE(M_n2o)

        
  REAL(r_std), SAVE ::   Maint_c = 0.0076    !! Maintenance coefficient of carbon (kgC/kgC/h)
                                        !! Table 4 Li et al., 2000
!$OMP THREADPRIVATE(Maint_c)
  REAL(r_std), SAVE ::   Yc = 0.503     !! Maximum growth yield on soluble carbon (kgC/kgC)
                                        !! Table 4 Li et al., 2000
!$OMP THREADPRIVATE(Yc)

  !! Coefficients used in the eq. defining the response of N-emission to clay fraction (-)
  !! from  Table 4, Li et al. 2000
  REAL(r_std), SAVE ::   F_clay_0 = 0.13   
!$OMP THREADPRIVATE(F_clay_0) 
  REAL(r_std), SAVE ::   F_clay_1 = -0.079
!$OMP THREADPRIVATE(F_clay_1)


  REAL(r_std), SAVE ::   ratio_nh4_fert = 0.5  !! Proportion of ammonium in the fertilizers (ammo-nitrate) 
                                           
!$OMP THREADPRIVATE(ratio_nh4_fert)

  REAL(r_std), SAVE ::   cn_ratio_manure = 13.7  !! C:N ratio of organic fertilizer (average over table 3-1-1 from Fuchs et al. 2014) 
 
!$OMP THREADPRIVATE(cn_ratio_manure)

  ! 2. Arrays
  REAL(r_std), SAVE, DIMENSION(2)  :: vmax_uptake = (/ 3. , 3. /) !! Vmax of nitrogen uptake by plants 
                                                                  !! for Ammonium (ind.1) and Nitrate (ind.2)
                                                                  !! (in umol (g DryWeight_root)-1 h-1)
                                                                  !! from  Kronzucker et al. (1995, 1996)
!$OMP THREADPRIVATE(vmax_uptake)

  REAL(r_std), SAVE, DIMENSION(2)  :: K_N_min = (/ 30., 30. /)    !! [NH4+] (resp. [NO3-]) for which the Nuptake 
                                                                  !! equals vmax/2.   (umol per litter)
                                                                  !! from Kronzucker, 1995
!$OMP THREADPRIVATE(K_N_min)

  REAL(r_std), SAVE, DIMENSION(2)  :: low_K_N_min = (/ 0.0002, 0.0002 /) !! Rate of N uptake not associated with 
                                                                         !! Michaelis- Menten Kinetics for Ammonium 
                                                                         !! (ind.1) and Nitrate (ind.2)
                                                                         !! from Kronzucker, 1995 ((umol)-1)
!$OMP THREADPRIVATE(low_K_N_min)

  REAL(r_std), SAVE :: timing_fert    = 183    !! Number of days on which fertilisation appears 
!$OMP THREADPRIVATE(timing_fert)

  REAL(r_std), SAVE      :: emm_fac = 0.2         !! Factor for reducing NH3 emission  
!$OMP THREADPRIVATE(emm_fac)
  REAL(r_std), SAVE      :: fact_kn_no = 0.012    !! Factor for adusting kn constant for NOx production
!$OMP THREADPRIVATE(fact_kn_no)
  REAL(r_std), SAVE      :: fact_kn_n2o = 0.04    !! Factor for adusting kn constant for N2O production
!$OMP THREADPRIVATE(fact_kn_n2o)
  REAL(r_std), SAVE      :: kfwdenit = -5.        !! Factor for adjusting sensitivity of denitrification to water content                     
!$OMP THREADPRIVATE(kfwdenit)  
 REAL(r_std), SAVE      :: fwdenitfc = 0.05       !! Value at field capacity of the sensitivity function of denitrification to water content                     
!$OMP THREADPRIVATE(fwdenitfc)          
  REAL(r_std), SAVE      :: fracn_drainage = 1.0  !! Fraction of NH3/NO3 loss by drainage 
!$OMP THREADPRIVATE(fracn_drainage)
  REAL(r_std), SAVE      :: fracn_runoff = 1.0  !! Fraction of NH3/NO3 loss by runoff
!$OMP THREADPRIVATE(fracn_runoff)


  !! Other N-related parameters
  REAL(r_std), SAVE                                  :: Dmax = 0.25      !! Parameter te be clarified (what it is, units, ...)
                                                                         !! used in stomate_growth_fun_all
!$OMP THREADPRIVATE(Dmax)
  REAL(r_std), SAVE :: reserve_time_tree = 30.     !! Maximum number of days during which
                                                   !! carbohydrate reserve may be used for 
                                                   !! trees (days)
!$OMP THREADPRIVATE(reserve_time_tree)
  REAL(r_std), SAVE :: reserve_time_grass = 20.    !! Maximum number of days during which
                                                   !! carbohydrate reserve may be used for 
                                                   !! grasses (days)
!$OMP THREADPRIVATE(reserve_time_grass)

  !
  ! stomate_season.f90 
  !

  ! 1. Scalar

  REAL(r_std), SAVE :: gppfrac_dormance = 0.2  !! report maximal GPP/GGP_max for dormance (0-1, unitless)
!$OMP THREADPRIVATE(gppfrac_dormance)
  REAL(r_std), SAVE :: tau_climatology = 20.   !! tau for "climatologic variables (years)
!$OMP THREADPRIVATE(tau_climatology)
  REAL(r_std), SAVE :: hvc1 = 0.019            !! parameters for herbivore activity (unitless)
!$OMP THREADPRIVATE(hvc1)
  REAL(r_std), SAVE :: hvc2 = 1.38             !! parameters for herbivore activity (unitless)
!$OMP THREADPRIVATE(hvc2)
  REAL(r_std), SAVE :: leaf_frac_hvc = 0.33    !! leaf fraction (0-1, unitless)
!$OMP THREADPRIVATE(leaf_frac_hvc)
  REAL(r_std), SAVE :: tlong_ref_max = 303.1   !! maximum reference long term temperature (K)
!$OMP THREADPRIVATE(tlong_ref_max)
  REAL(r_std), SAVE :: tlong_ref_min = 253.1   !! minimum reference long term temperature (K)
!$OMP THREADPRIVATE(tlong_ref_min)

  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: ncd_max_year = 3.
!$OMP THREADPRIVATE(ncd_max_year)
  REAL(r_std), SAVE :: gdd_threshold = 5.
!$OMP THREADPRIVATE(gdd_threshold)
  REAL(r_std), SAVE :: green_age_ever = 2.
!$OMP THREADPRIVATE(green_age_ever)
  REAL(r_std), SAVE :: green_age_dec = 0.5
!$OMP THREADPRIVATE(green_age_dec)

  INTEGER(i_std), SAVE :: ncirc = 1                  !! Number of circumference classes used to calculate C allocation
!$OMP THREADPRIVATE(ncirc)
  REAL(r_std), PARAMETER :: kilo_to_unit = 1.0E03   !! Convert Kilo to unit


  REAL(r_std), SAVE :: sync_threshold = 0.0001     !! The threshold above which a warning is generated when the
!$OMP THREADPRIVATE(sync_threshold)
  LOGICAL,PARAMETER :: ld_biomass=.FALSE.   ! a flag to turn on debug statements
  INTEGER(i_std), SAVE        :: test_pft = 4                              !! Number of PFT for which detailed output 
!$OMP THREADPRIVATE(test_pft)

  INTEGER(i_std), SAVE        :: test_grid = 1                                 !! Number of the grid square for which detailed output 
!$OMP THREADPRIVATE(test_grid)

  LOGICAL,PARAMETER                       :: ld_stop=.FALSE.      ! a flag to turn on some stop statements.
  LOGICAL,SAVE                            :: ld_alloc=.FALSE.     ! a flag to turn on debug statements
!$OMP THREADPRIVATE(ld_alloc)
  LOGICAL,PARAMETER                       :: ld_warn=.FALSE.      ! a flag to turn on various warnings
  LOGICAL,PARAMETER                       :: jr_nextstep = .FALSE.   ! set this to TRUE to activate the 
  LOGICAL,PARAMETER                       :: ld_massbal=.FALSE.   ! a flag to turn on debug statements
  INTEGER(i_std), PARAMETER :: ipoolchange = 5           !! change in pool size i.e. change in biomass
  INTEGER(i_std), PARAMETER :: ilat2in = 4               !! incoming lateral flux i.e. N deposition for the land
  INTEGER(i_std), PARAMETER :: ilat2out = 3              !! outgoing lateral flux i.e. DOC leaching for the litter routine
  INTEGER(i_std), PARAMETER :: iatm2land = 1             !! atmosphere to land fluxes such as GPP and co2_2_bm
  INTEGER(i_std), PARAMETER :: iland2atm = 2             !! land to atmosphere fluxes such as Rh, Ra and product decomposition
  INTEGER(i_std), PARAMETER :: nmbcomp = 5               !! The total nomber of components in our mass balance check
  

  REAL(r_std), SAVE :: max_delta_KF = 0.1          !! Maximum change in KF from one time step to another (m)
                                                   !! This is a bit arbitrary. 
!$OMP THREADPRIVATE(max_delta_KF)

  REAL(r_std), SAVE :: maint_from_gpp = 0.8        !! Some carbon needs to remain to support the growth, hence, 
                                                   !! respiration will be limited. In this case resp_maint 
                                                   !! (gC m-2 dt-1) should not be more than 80% (::maint_from_gpp) 
                                                   !! of the GPP (gC m-2 s-1)
!$OMP THREADPRIVATE(maint_from_gpp)
 
  REAL(r_std), PARAMETER :: m2_to_ha = 10000.       !! Conversion from m2 to hectares
  REAL(r_std), PARAMETER :: ha_to_m2 = 0.0001       !! Conversion from hectares (forestry) to m2 (rest of the code)

  ! Soil carbon related variables
  LOGICAL, SAVE                   :: soilc_isspinup = .FALSE.            !! Activate soil carbon spinup
!$OMP THREADPRIVATE(soilc_isspinup)
  REAL(r_std), PARAMETER          :: O2_init_conc = 298.3                !! gO2/m**3 mean for Cherskii
  REAL(r_std), PARAMETER          :: CH4_init_conc = 0.001267            !! gCH4/m**3 mean for Cherskii
  REAL(r_std), PARAMETER          :: z_root_max = 0.5                    !! Depth at which litter carbon input decays e-fold (root depth); 0.5 for compar w/WH
  REAL(r_std), PARAMETER          :: diffO2_air = 1.596E-5               !! oxygen diffusivity in air (m**2/s)
  REAL(r_std), PARAMETER          :: diffO2_w = 1.596E-9                 !! oxygen diffusivity in water (m**2/s)
  REAL(r_std), PARAMETER          :: O2_surf = 0.209                     !! oxygen concentration in surface air (molar fraction)
  REAL(r_std), PARAMETER          :: diffCH4_air = 1.702E-5              !! methane diffusivity in air (m**2/s)
  REAL(r_std), PARAMETER          :: diffCH4_w = 2.0E-9                  !! methane diffusivity in water (m**2/s)
  REAL(r_std), PARAMETER          :: CH4_surf = 1700.E-9                 !! methane concentration in surface air (molar fraction)
  REAL(r_std), SAVE               :: tetasat =  .5                       !! volumetric water content at saturation (porosity)
!$OMP THREADPRIVATE(tetasat)
  REAL(r_std), SAVE               :: tetamoss =  0.92                    !! porosity of moss
!$OMP THREADPRIVATE(tetamoss)
  REAL(r_std), SAVE               :: rho_moss                            !! density of moss
!$OMP THREADPRIVATE(rho_moss)
  REAL(r_std), PARAMETER          :: zmoss = 0.2                         !! thickness of moss layer (in permafrost regions,m) 0. ! 0.001 DKtest for compar w/WH
  REAL(r_std), PARAMETER          :: h_snowmoss = 0.2                    !! snow height above which we consider the moss layer to be to compressed to be effective (m)
  REAL(r_std), PARAMETER          :: BunsenO2 = 0.038                    !! Bunsen coefficient for O2 (10C, 1bar)
  REAL(r_std), PARAMETER          :: BunsenCH4 = 0.043                   !! Bunsen coefficient for CH4 (10C, Wiesenburg et Guinasso, Jr., 1979)
  REAL(r_std), PARAMETER          :: ebuthr = 0.9                        !! Soil humidity threshold for ebullition
  REAL(r_std), PARAMETER          :: wCH4 = 16.                          !! molar weight of CH4 (g/mol)
  REAL(r_std), PARAMETER          :: wO2 = 32.                           !! molar weight of O2 (g/mol)
  REAL(r_std), PARAMETER          :: wC = 12.                            !! molar weight of C (g/mol)
  REAL(r_std), PARAMETER          :: avm = .01                           !! minimum air volume (m**3 air/m**3 soil)
  REAL(r_std), PARAMETER          :: hmin_tcalc = .001                   !! minimum total snow layer thickness below which we ignore diffusion across the snow layer



                           !-----------------------------!
                           !  STICS PARAMETERS !
                           !-----------------------------!

  ! STICS NITROGEN
  ! IN THIS VERSION, WE DO NOT HAVE NITROGEN PROCESS, SO NO NITROGEN LIMITATION
  !REAL(r_std), SAVE :: innlai = 0.8 ! nitrogen stress for lai
  REAL(r_std), SAVE :: inns = 1.   ! nitrogen stress for grain production
!$OMP THREADPRIVATE(inns)
  REAL(r_std), SAVE :: innsenes = 1. ! nitrogen stress for senescence
!$OMP THREADPRIVATE(innsenes)
  REAL(r_std), SAVE :: exolai = 1.  !! water logging effects on lai [0 1], 1, no slowing effect
!$OMP THREADPRIVATE(exolai)
  REAL(r_std), SAVE :: resperenne = 0.0  !! C crop reserve, during the cropping season, or during the intercrop period (for perenial crops) // t ha-1
!$OMP THREADPRIVATE(resperenne)
  INTEGER(i_std), SAVE :: numcult = 1
!$OMP THREADPRIVATE(numcult)
  INTEGER(i_std), SAVE :: vlength = 70
!$OMP THREADPRIVATE(vlength)
  !INTEGER(i_std), SAVE :: maxgs = 300 ! This is an arbitrary setting of maximum
  !growing length for annual crop to avoid a too long growing season length
  INTEGER(i_std), SAVE :: matrec = 15 ! This is an arbitrary setting for a period between nmat and nrec for grain drying and nutrition dynamics
!$OMP THREADPRIVATE(matrec)

  REAL(r_std), SAVE :: tetstomate = 0.53  ! m3/m3 soil water contents
!$OMP THREADPRIVATE(tetstomate)
  REAL(r_std), SAVE :: teturg = 0.73  !
!$OMP THREADPRIVATE(teturg)
  REAL(r_std), SAVE :: tetsen = 0.4  !
!$OMP THREADPRIVATE(tetsen)

  LOGICAL, SAVE :: nitrogen_use = .FALSE.       !! IF we consider the nitrogen fertilization (true/false) 
!$OMP THREADPRIVATE(nitrogen_use)
  LOGICAL, SAVE :: fix_nfert = .TRUE.       !! IF we use the average nitrogen fertilization (true/false) 
!$OMP THREADPRIVATE(fix_nfert)
  LOGICAL, SAVE :: iplt_1d = .TRUE.  !! do we use the imposed planting date
!$OMP THREADPRIVATE(iplt_1d)
  REAL(r_std), SAVE :: N_effect = 0.65         !! this is the maximum additive effect of nitrogen fertilization on Vcmax
!$OMP THREADPRIVATE(N_effect)
  REAL(r_std), SAVE :: prc_residual = 0.5      !! percentage of straws leaving on the croplands
!$OMP THREADPRIVATE(prc_residual)
  INTEGER(i_std), SAVE :: codesla = 2     !! 1 for old sla calc, 2 for new sla calc



END MODULE constantes_var

! =================================================================================================================================
! MODULE 	: constantes_soil
!
! CONTACT       : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE       : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF         "constantes_soil" module contains subroutine to initialize the parameters related to soil and hydrology.
!!
!!\n DESCRIPTION : "constantes_soil" module contains subroutine to initialize the parameters related to soil and hydrology.
!!                 This module alos USE constates_soil and can therfor be used to acces the subroutines and the constantes.
!!                 The constantes declarations can also be used seperatly with "USE constantes_soil_var".
!!
!! RECENT CHANGE(S): 
!!
!! REFERENCE(S)	:
!!
!! SVN          :
!! $HeadURL: $
!! $Date: $
!! $Revision: $
!! \n
!_ ================================================================================================================================

MODULE constantes_soil

  USE constantes_soil_var
  USE constantes
  USE ioipsl_para 

  IMPLICIT NONE

CONTAINS


!! ================================================================================================================================
!! SUBROUTINE   : config_soil_parameters
!!
!>\BRIEF        This subroutine reads in the configuration file all the parameters related to soil and hydrology. 
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S) :
!!
!! FLOWCHART    :
!! \n
!_ ================================================================================================================================

  SUBROUTINE config_soil_parameters()

    USE ioipsl

    IMPLICIT NONE

    !! 0. Variables and parameters declaration

    !! 0.4 Local variables 

    INTEGER(i_std), PARAMETER      :: error_level = 3         !! Switch to 2 to turn fatal errors into warnings.(1-3, unitless)
    REAL(r_std)                    :: fr_center_c             !! Local value of fr_center in Celsius
    LOGICAL                        :: ok_freeze               !! Local variable used to set default values for all flags 
    !! controling the soil freezing scheme

    !_ ================================================================================================================================

    ! Following initializations are only done for option impose_param
    IF ( ok_sechiba .AND. impose_param ) THEN


       !Config Key   = SNOW_HEAT_COND
       !Config Desc  = Thermal Conductivity of snow
       !Config If    = OK_SECHIBA  
       !Config Def   = 0.3
       !Config Help  = 
       !Config Units = [W.m^{-2}.K^{-1}]
       CALL getin_p("SNOW_HEAT_COND",sn_cond)

       !! Check
       IF ( sn_cond <= zero ) THEN
          CALL ipslerr_p(error_level, "config_soil_parameters.", &
               &     "Wrong parameter value for SNOW_HEAT_COND.", &
               &     "This parameter should be positive. ", &
               &     "Please, check parameter value in run.def. ")
       END IF


       !Config Key   = SNOW_DENSITY
       !Config Desc  = Snow density for the soil thermodynamics 
       !Config If    = OK_SECHIBA 
       !Config Def   = 330.0
       !Config Help  = 
       !Config Units = [-] 
       CALL getin_p("SNOW_DENSITY",sn_dens)

       !! Check parameter value (correct range)
       IF ( sn_dens <= zero ) THEN
          CALL ipslerr_p(error_level, "config_soil_parameters.", &
               &     "Wrong parameter value for SNOW_DENSITY.", &
               &     "This parameter should be positive. ", &
               &     "Please, check parameter value in run.def. ")
       END IF


       !! Calculation of snow capacity
       !! If sn_dens is redefined by the user, sn_capa needs to be reset
       sn_capa = 2100.0_r_std*sn_dens


       !Config Key   = NOBIO_WATER_CAPAC_VOLUMETRI
       !Config Desc  = 
       !Config If    = 
       !Config Def   = 150.
       !Config Help  = 
       !Config Units = [s/m^2]
       CALL getin_p('NOBIO_WATER_CAPAC_VOLUMETRI',mx_eau_nobio)

       !! Check parameter value (correct range)
       IF ( mx_eau_nobio <= zero ) THEN
          CALL ipslerr_p(error_level, "config_soil_parameters.", &
               &     "Wrong parameter value for NOBIO_WATER_CAPAC_VOLUMETRI.", &
               &     "This parameter should be positive. ", &
               &     "Please, check parameter value in run.def. ")
       END IF


       !Config Key   = SECHIBA_QSINT 
       !Config Desc  = Interception reservoir coefficient
       !Config If    = OK_SECHIBA 
       !Config Def   = 0.02
       !Config Help  = Transforms leaf area index into size of interception reservoir
       !Config         for slowproc_derivvar or stomate
       !Config Units = [m]
       CALL getin_p('SECHIBA_QSINT',qsintcst)

       !! Check parameter value (correct range)
       IF ( qsintcst <= zero ) THEN
          CALL ipslerr_p(error_level, "config_soil_parameters.", &
               &     "Wrong parameter value for SECHIBA_QSINT.", &
               &     "This parameter should be positive. ", &
               &     "Please, check parameter value in run.def. ")
       END IF


    END IF ! IF ( ok_sechiba .AND. impose_param ) THEN



    !! Variables related to soil freezing in thermosoil module
    !
    !Config Key  = OK_FREEZE
    !Config Desc = Activate the complet soil freezing scheme
    !Config If   = OK_SECHIBA 
    !Config Def  = TRUE
    !Config Help = Activate soil freezing thermal effects. Activates soil freezing hydrological effects in CWRR scheme.
    !Config Units= [FLAG]

    ! ok_freeze is a flag that controls the default values for several flags controling 
    ! the different soil freezing processes
    ! Set ok_freeze=true for the complete soil freezing scheme
    ! ok_freeze is a local variable only used in this subroutine
    ok_freeze = .TRUE.
    CALL getin_p('OK_FREEZE',ok_freeze)


    !Config Key  = READ_REFTEMP
    !Config Desc = Initialize soil temperature using climatological temperature
    !Config If   = 
    !Config Def  = True/False depening on OK_FREEZE
    !Config Help = 
    !Config Units= [FLAG]

    IF (ok_freeze) THEN
       read_reftemp = .TRUE.
    ELSE
       read_reftemp = .FALSE.
    END IF
    CALL getin_p ('READ_REFTEMP',read_reftemp)

    !Config Key  = USE_INITSOM
    !Config Desc = Initialize SOM using SOM map
    !Config If   = OK_SOIL_CARBON_DISCRETIZATION
    !Config Def  = False
    !Config Help = 
    !Config Units= [FLAG]
    use_initsom = .FALSE.
    CALL getin_p ('USE_INITSOM',use_initsom)
    

    IF (use_initsom) THEN
       !Config Key   = INITSOM_FRAC_A
       !Config Desc  = Fraction of initSOM to put in active C pool
       !Config If    = USE_INITSOM
       !Config Def   = 0.1
       !Config Help  =
       !Config Units = [FILE]
       initSOM_frac_a = 0.1
       CALL getin_p('INITSOM_FRAC_A',initSOM_frac_a)
       
       !Config Key   = INITSOM_FRAC_S
       !Config Desc  = Fraction of initSOM to put in slow C pool
       !Config If    = USE_INITSOM
       !Config Def   = 0.2
       !Config Help  =
       !Config Units = [FILE]
       initSOM_frac_s = 0.2
       CALL getin_p('INITSOM_FRAC_S',initSOM_frac_s)
       
       !Config Key   = INITSOM_FRAC_P
       !Config Desc  = Fraction of initSOM to put in passive C pool
       !Config If    = USE_INITSOM
       !Config Def   = 0.7
       !Config Help  =
       !Config Units = [FILE]
       initSOM_frac_p = 0.7
       CALL getin_p('INITSOM_FRAC_P',initSOM_frac_p)

    END IF

       
    !Config Key  = OK_FREEZE_THERMIX
    !Config Desc = Activate thermal part of the soil freezing scheme
    !Config If   = 
    !Config Def  = True if OK_FREEZE else false
    !Config Help = 
    !Config Units= [FLAG]

    IF (ok_freeze) THEN
       ok_freeze_thermix = .TRUE.
    ELSE
       ok_freeze_thermix = .FALSE.
    END IF
    CALL getin_p ('OK_FREEZE_THERMIX',ok_freeze_thermix)


    !Config Key  = OK_ECORR
    !Config Desc = Energy correction for freezing
    !Config If   = OK_FREEZE_THERMIX
    !Config Def  = True if OK_FREEZE else false
    !Config Help = Energy conservation : Correction to make sure that the same latent heat is 
    !Config        released and consumed during freezing and thawing
    !Config Units= [FLAG]
    IF (ok_freeze) THEN
       ok_Ecorr = .TRUE.
    ELSE
       ok_Ecorr = .FALSE.
    END IF
    CALL getin_p ('OK_ECORR',ok_Ecorr)
    IF (ok_Ecorr .AND. .NOT. ok_freeze_thermix) THEN
       CALL ipslerr_p(3,'thermosoil_init','OK_ECORR cannot be activated without OK_FREEZE_THERMIX', &
            'Adapt run parameters with OK_FREEZE_THERMIX=y','')
    END IF

    !Config Key  = OK_FREEZE_THAW_LATENT_HEAT
    !Config Desc = Activate latent heat part of the soil freezing scheme
    !Config If   = 
    !Config Def  = FALSE 
    !Config Help = 
    !Config Units= [FLAG]

    ok_freeze_thaw_latent_heat = .FALSE.
    CALL getin_p ('OK_FREEZE_THAW_LATENT_HEAT',ok_freeze_thaw_latent_heat)


    !Config Key = fr_dT
    !Config Desc = Freezing window width
    !Config If = OK_SECHIBA
    !Config Def = 2.0
    !Config Help = 
    !Config Units = [K] 
    fr_dT=2.0
    CALL getin_p('FR_DT',fr_dT)

    !Config Key = FR_CENTER_C
    !Config Desc = Center value of the freezing window in degree Celsius
    !Config If = OK_SECHIBA
    !Config Def = 0.0
    !Config Help =
    !Config Units = [°C]
    fr_center_c=0.0
    CALL getin_p('FR_CENTER_C',fr_center_c)
    fr_center=fr_center_c+ZeroCelsius

    !Config Key   = SOILC_MAX
    !Config Desc  = Soil carbon above which soil thermal properties equals to organic soil properties 
    !Config If    = OK_SOIL_CARBON_DISCRETIZATION and USE_SOILC_TEMPDIFF
    !Config Def   = 130000
    !Config Help  = 
    !Config Units = [gC/m3] 
    CALL getin_p("SOILC_MAX", soilc_max)


    !! Variables related to soil Freezing in hydrol module

    !Config Key  = OK_FREEZE_CWRR
    !Config Desc = CWRR freezing scheme by I. Gouttevin
    !Config If   = 
    !Config Def  = True if OK_FREEZE else false
    !Config Help =
    !Config Units= [FLAG]

    IF (ok_freeze) THEN
       ok_freeze_cwrr = .TRUE.
    ELSE
       ok_freeze_cwrr = .FALSE.
    END IF
    CALL getin_p('OK_FREEZE_CWRR',ok_freeze_cwrr)


    IF (ok_freeze_cwrr) THEN
       !Config Key  = OK_THERMODYNAMICAL_FREEZING
       !Config Desc = Calculate frozen fraction thermodynamically 
       !Config If   = OK_FREEZE_CWRR
       !Config Def  = True
       !Config Help = Calculate frozen fraction thermodynamically if true,
       !Config      = else calculate frozen fraction linearly 
       !Config Units= [FLAG]
       ok_thermodynamical_freezing = .TRUE.
       CALL getin_p('OK_THERMODYNAMICAL_FREEZING',ok_thermodynamical_freezing)
    END IF


    !Config Key   = CHECK_CWRR
    !Config Desc  = Calculate diagnostics to check CWRR water balance
    !Config Def   = n
    !Config If    = 
    !Config Help  = Diagnostics from module hydrol. The verifictaions are done in post-treatement
    !Config Units = [FLAG]
    check_cwrr = .FALSE.
    CALL getin_p('CHECK_CWRR', check_cwrr)


    !Config Key  = OK_MOSS
    !Config Desc = Activate moss thermal effect
    !Config If   = OK_SOIL_CARBON_DISCRETIZATION 
    !Config Def  = FALSE
    !Config Help = Activate moss thermal effect (thermal capacity and conductivity).
    !Config Units= [FLAG]
    ok_moss = .FALSE.
    CALL getin_p('OK_MOSS',ok_moss)

    IF (ok_moss .AND. ok_soil_carbon_discretization==.FALSE.) THEN
       CALL ipslerr_p(3, "config_soil_parameters.", &
            "OK_MOSS can not be activated without OK_SOIL_CARBON_DISCRETIZATION=y",&
            "Please, change parameter value in run.def","")
    END IF

!    IF (ok_moss) THEN
       !Config Key = MOSS_LAYER_THICKNESS
       !Config Desc = Moss layer thickness
       !Config If = OK_SOIL_CARBON_DISCRETIZATION AND OK_MOSS
       !Config Def = 0.05
       !Config Help = 
       !Config Units = [m] 
       moss_layer_thickness=0.05
       CALL getin_p('MOSS_LAYER_THICKNESS',moss_layer_thickness)
       
       !Config Key = SOIL_MOSS_DEPTH
       !Config Desc = Maximum depth over which to compute moss thermal properties
       !Config If = OK_SOIL_CARBON_DISCRETIZATION AND OK_MOSS
       !Config Def = 0.05
       !Config Help =
       !Config Units = [m]
       soil_moss_depth=0.05
       CALL getin_p('SOIL_MOSS_DEPTH',soil_moss_depth)
!    END IF

  END SUBROUTINE config_soil_parameters


END MODULE constantes_soil

!***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular puR_companionose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************

! added little comment

      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use auto_diff
      use energy

      implicit none


    ! These variables can be saved in photos and restored at restarts
      real(dp) :: Orbital_separation, Deltar, Deltar_tides, R_bondi, R_influence
      real(dp) :: stop_age

    ! the routines that take care of doing the save/restore are the following:
    ! alloc_extra_info and unpack_extra_info << called by extras_startup
    ! store_extra_info << called by extras_finish_step
    ! these routines call move_extra_info.
    ! it must know about each of your variables to be saved/restored.
    ! so edit move_extra_info when you change the set of variables.

    ! these routines are called by the standard run_star check_model
      contains

      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

       ! this is the place to set any procedure pointers you want to change
       ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)

       ! Uncomment these lines if you wish to use the functions in this file,
       ! otherwise we use a null_ version which does nothing.
         s% other_energy => energy_routine

         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

       ! Once you have set the function pointers you want,
       ! then uncomment this (or set it in your star_job inlist)
       ! to disable the printed warning message,
       ! s% job% warn_run_star_extras =.false.


      end subroutine extras_controls

      subroutine energy_routine(id, ierr)
        integer, intent(in) :: id
        integer, intent(out) :: ierr
        logical :: restart, first
        type (star_info), pointer :: s
        integer :: k, nz
        integer :: krr_center, krr_bottom_bondi, krr_top_bondi, krr_bottom_companion, krr_top_companion
        real(dp) :: e_orbit, M_companion, R_companion, area, de, sound_speed, R_influence
        real(dp) :: rr, v_kepler, rho_bar_companion, rho_bar_bondi, rho_bar_drag
        real(dp) :: dmsum_companion,dmsum_bondi,dmsum_drag, de_heat
        real(dp) :: f_disruption
        real(dp) :: penetration_depth
        real(dp) :: t_tide
        real(dp) :: de_orbital_change, enclosed_mass
        ierr = 0

      ! Reads model infos from star structure s. Initialize variables.
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return


        nz = s% nz               ! Mesh size (primary)
         
        do k = 1, nz
            s% extra_heat(k) = 0d0    ! Initialize extra_heat vector
        end do
        

      ! Initialize injected energy and radial coordinate change
        de = 0d0
        Deltar = 0d0
        f_disruption = 0d0
        R_bondi = 0d0


      ! Mass and radius of injested companion from inlist. Also include a stop point for companion (x_ctrl(3) in Rsun) in case we want to stop before destruction.
        M_companion = s% x_ctrl(1) * Msun
        R_companion = s% x_ctrl(2) * Rsun


      ! Orbital_separation is the coordinate of the planet's center wrt the primary's core.
      ! If it's a restart, MESA will remember the radial location of the
      ! companion, Orbital_separation, from a photo. This is because we are moving Orbital_separation data in
      ! photos using 'move_extra_info' and this data is retrieved in 'extras_startup' using 'unpack_extra_info'

      ! Calculate orbital keplerian velocity of the companion (we assume circular orbits)
      ! If the companion's centre is outside the primary this is easyly done, but if it is inside, we need to use
      ! only the mass of the primary inside the orbit, so we need to locate the index where the companion core is.
      ! Calculate the bondi radius of the companion using a sound speed of 10 km/s outside the star (typical ISM)
        
        krr_center=1

        if (Orbital_separation > s% r(1)) then
            call orbital_velocity(s% m(1), Orbital_separation, v_kepler)
            sound_speed = 10. * 1.d5 ! set c_sound to be ISM in cgs
        else
            do while (krr_center >= 1 .and. krr_center < nz .and. s% r(krr_center) >= Orbital_separation)
                krr_center = krr_center + 1
            end do
            !write(*,*) 'Orb sep, r(krr_center) ',Orbital_separation, s% r(krr_center)
            !Orbital_separation = s% r(krr_center) ! This is to accommodate for errors in energy calculation introduced by gridding
            call orbital_velocity(s% m(krr_center), s% r(krr_center), v_kepler)
            sound_speed = s% csound(krr_center)
        endif

        call bondi_radius (M_companion, sound_speed, v_kepler, R_bondi)
        R_influence = max(R_bondi,R_companion)  ! Choose radius to be used for drag routine (R_bondi -> Gravodrag, R_companion -> Aerodynamic drag)


      ! Find gridpoint corresponding to the location of the engulfed companion bottom, center and top
        krr_bottom_companion=1
        do while (krr_bottom_companion >= 1 .and. &
                  krr_bottom_companion < nz .and. &
                  s% r(krr_bottom_companion) >= Orbital_separation-R_companion)
            krr_bottom_companion = krr_bottom_companion + 1
        end do
        krr_center = krr_bottom_companion
        do while (krr_center >=2 .and. s% r(krr_center) < Orbital_separation)
           krr_center = krr_center - 1
        end do
        krr_top_companion = krr_center
        do while (krr_top_companion >= 2 .and. s% r(krr_top_companion) < Orbital_separation+R_companion)
            krr_top_companion = krr_top_companion - 1
        end do


       ! Find gridpoint corresponding to the location of the engulfed companion's Bondi radius,  bottom, center and top
         krr_bottom_bondi=1
         do while (krr_bottom_bondi >= 1 .and. krr_bottom_bondi < nz .and. s% r(krr_bottom_bondi) >= Orbital_separation-R_bondi)
             krr_bottom_bondi = krr_bottom_bondi + 1
         end do
         krr_center = krr_bottom_bondi
         do while (krr_center >=2 .and. s% r(krr_center) < Orbital_separation)
            krr_center = krr_center - 1
         end do
         krr_top_bondi = krr_center
         do while (krr_top_bondi >= 2 .and. s% r(krr_top_bondi) < Orbital_separation+R_bondi)
             krr_top_bondi = krr_top_bondi - 1
         end do


      ! Calculate mass contained in the spherical shell occupied by the companion (shellular approximation)
      ! and the mass-weighted density of the region of impact for drag calculation
        dmsum_companion = sum(s% dm(krr_top_companion:krr_bottom_companion))
        dmsum_bondi = sum(s% dm(krr_top_bondi:krr_bottom_bondi))
        rho_bar_companion = dot_product &
                            (s% rho(krr_top_companion:krr_bottom_companion), &
                             s% dm(krr_top_companion:krr_bottom_companion))/dmsum_companion
        rho_bar_bondi = dot_product &
                            (s% rho(krr_top_bondi:krr_bottom_bondi), &
                             s% dm(krr_top_bondi:krr_bottom_bondi))/dmsum_bondi
        if (R_bondi >= R_companion) then
           dmsum_drag = dmsum_bondi
           rho_bar_drag = rho_bar_bondi
        else
           dmsum_drag = dmsum_companion
           rho_bar_drag = rho_bar_companion
        endif
      ! write(*,*)'indeces,bottom,centre,top',krr_bottom_bondi,krr_bottom_companion,krr_center,krr_top_companion,krr_top_bondi

      ! Check if the companion has been destroyed by ram pressure (f>1). This is important only for planets.
        f_disruption = check_disruption(M_companion,R_companion,v_kepler,rho_bar_drag)

      ! Calculate area used for drag calculation 
        penetration_depth = 0d0
        area = 0d0 ! Initialize cross section of companion (physical or Bondi) for calculating aerodynamic or gravitational drag

      ! Do the calculation only if this is a grazing collision and if the planet has not been destroyed yet
        if (Orbital_separation > s% r(1) + R_influence) then
            penetration_depth = 0.d0
        else
            penetration_depth = penetration_depth_function(R_influence,s% r(1), Orbital_separation)
        endif


        if (penetration_depth >= 0.0 .and. (Orbital_separation >= (s% r(1) - R_influence)) .and. (f_disruption <= 1d0)) then
            ! Calculate intersected area. Rstar-rr is x in sketch
              area = intercepted_area (penetration_depth, R_influence)
            !  write(*,*) 'Grazing Collision. Engulfed area fraction: ', s% model_number, area/(pi * pow(R_influence, 2.0))
        else
            ! Full engulfment. Cross section area = Planet area
              area = pi * pow(R_influence, 2d0)
            !  write(*,*) 'Full engulfment. R_influence, area',s% model_number,R_influence/Rsun,area
        end if

      !  write(*,'(a,i5,4f11.6,3e14.5)') &
      !          'Orbital_separation, R_bondi, R_influence, penetration depth, rho_bar_bondi, rho_bar_companion, area ',&
      !           s% model_number, Orbital_separation/Rsun, R_bondi/Rsun, R_influence/Rsun, &
      !           penetration_depth/Rsun, rho_bar_bondi, rho_bar_companion, area



        !########################### TIDES ######################################

        ! Calculate tidal timescale (according to Hansen et al. 2010, which uses Hut formalism)
        ! When Orbital_separation < R_star we assume that only the mass and radius of the star within the orbital separation play a role
        if (s% x_ctrl(8) > 0.0) then
            call tidal_timescale(s% m(krr_center), M_companion, s% r(krr_center), &
             R_companion, Orbital_separation, s% x_ctrl(8), t_tide)
            Deltar_tides = (s% dt/t_tide) * Orbital_separation
          else
            Deltar_tides = 0.0
            t_tide = 0.0 ! This should be +inf
        end if

        if (.not. s% x_logical_ctrl(1) .and. Orbital_separation <= s% r(1)) then ! No tides if a<R (if s% x_logical_ctrl(1) = .false.)
          Deltar_tides = 0.0
          t_tide = 0.0 ! This should be +inf*
        end if

        !########################### END TIDES ######################################

      ! If the companion has not been destroyed by ram pressure, deposit drag luminosity and heat the envelope
      ! Spread in the region occupied by the planet or by Bondi sphere, whichever is larger. Update radial coordinate of the engulfed planet too.
        if ( f_disruption <= 1d0 ) then

            ! Note we use s% r(krr_bottom)+R_companion instead of s% r(krr_center) because during grazing phase krr_center = krr_bottom
              call drag (s% m(krr_center), M_companion, area, rho_bar_drag, s% dt, Orbital_separation, de, Deltar)

              de_orbital_change = calculate_orbital_energy(s% m(krr_center),M_companion,Orbital_separation)
              ! Loop to find approximate enclosed mass at next step
              
             ! write(*,*) 'k center',krr_center
             ! do while (krr_center >= 1 .and. &
             !     krr_center < nz .and. &
             !     s% r(krr_center) > Orbital_separation-Deltar-Deltar_tides)
             !   krr_center = krr_center + 1
             ! end do
               
              de_orbital_change = de_orbital_change - calculate_orbital_energy(s% m(krr_center),M_companion,Orbital_separation-Deltar-Deltar_tides)
              !write(*,*) 'k center next step, de', krr_center, de_orbital_change 
              write(*,*) 'de_orbital_change, de', de_orbital_change, de ! Slight discrepancy between these two because De is approximate 
              
              !de = de_orbital_change ! Set this to be the exact change in orbital energy 

            !  write(*,'(A,i4,2f10.4,2e15.4,f12.4,e12.4,e12.4)')'after call drag', s% model_number, s% m(krr_center)/Msun, &
						!		       area, rho_bar_drag, s% dt, Orbital_separation/Rsun, de, Deltar/Rsun

            ! If the planet has not been destroyed by ram pressure, deposit drag luminosity and heat the envelope
            ! Spread in the region occupied by the planet or by the Bondi sphere.
            ! Update radial coordinate of the engulfed planet in 'extras_finish_step' using Deltar.
              do k = min(krr_top_bondi,krr_top_companion), max(krr_bottom_bondi,krr_bottom_companion)
                 s% extra_heat(k) = (de/dmsum_drag/s% dt) ! Uniform heating (erg/g/sec)
              end do

              ! Calculate orbital energy

              call orbital_energy(s% m(krr_center), M_companion,Orbital_separation, e_orbit)
            !  write(*,*) 'Injected Energy / Orbital Energy: ', abs(de/e_orbit)
        else
              write(*,*) '***************** Planet destroyed at R/Rsun = ', Orbital_separation/Rsun,'*********************'
              Deltar = 0d0
              s% use_other_energy = .false.
              stop_age = s% star_age + 1d1 * s% kh_timescale
              write(*,*)'stop_age and KH timescale are',stop_age, s% kh_timescale
        endif

                
        ! Save enclosed mass at timestep
        enclosed_mass = s% m(krr_center)

      ! write(*,*) 'Tidal Timescale (yrs): ',t_tide/secyer, 'Tidal Da (Rsun): ', Deltar_tides/Rsun, s% dt

        ! Save variables for history

          s% xtra(1) = v_kepler/1d5                ! Orbital velocity (km/s)
          s% xtra(2) = Deltar                      ! Infall distance due to drag (cm)
          s% xtra(3) = de                          ! Injected energy (erg)
          s% xtra(4) = f_disruption                ! Disruption factor
          s% xtra(5) = area/(pi * pow(max(R_companion,R_bondi), 2d0))  ! Engulfed fraction
          s% xtra(6) = dmsum_drag/Msun             ! Heated mass (Msun)
          s% xtra(7) = R_bondi/Rsun                ! Bondi radius (Rsun)
          s% xtra(8) = sound_speed/1.d5            ! Sound speed (km/s)
          s% xtra(9) = t_tide/secyer               ! Tidal timescale (yrs)
          s% xtra(10) = Deltar_tides               ! Infall distance due to tides (cm)
          s% xtra(11) = Deltar/s% dt/1e5           ! Infall velocity (drag) (km/s)
          s% xtra(12) = enclosed_mass              ! Enclosed stellar mass at companion location  
          s% xtra(13) = rho_bar_drag               ! Average Density at companion location

      end subroutine energy_routine




        subroutine extras_startup(id, restart, ierr)
           integer, intent(in) :: id
           logical, intent(in) :: restart
           integer, intent(out) :: ierr
           real(dp) :: v_kepler, sound_speed, R_influence !, Orbital_separation_prior
           type (star_info), pointer :: s
           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return

           if (.not. restart) then 
            Orbital_separation = s% x_ctrl(6)*Rsun ! Set initial separation to inlist value 
            sound_speed = 10. * 1.d5 ! set c_sound to be ISM in cgs
            call orbital_velocity(s% m(1), Orbital_separation, v_kepler)
            call bondi_radius (s% x_ctrl(1)*Msun, sound_speed, v_kepler, R_bondi)
            R_influence = max(R_bondi, s% x_ctrl(2)*Rsun)
            stop_age = -101d0
            call alloc_extra_info(s)
          else ! it is a restart -> Unpack value of Orbital_separation from photo
             call unpack_extra_info(s)
          end if

         ! We need to increase the resolution around the area where the extra heat is deposited
         ! We will do this at the startup and also in the extra_check model, since the position
         ! of the companion will be changing
          write(*,*) 'From Startup', s% R_function2_param1, s% R_function2_param2, Orbital_separation, s%use_other_energy
          if (Orbital_separation <= s% r(1) .and. s% use_other_energy ) then
            s% R_function2_param1 = Orbital_separation/(s%r(1)/Rsun) + 2.0 *  s% x_ctrl(2) * Rsun/s%r(1)
            s% R_function2_param2 = Orbital_separation/(s%r(1)/Rsun) - 2.0 *  s% x_ctrl(2) * Rsun/s%r(1)
            write(*,*) 'From Startup', s% R_function2_param1, s% R_function2_param2
          endif
        end subroutine extras_startup


     

      integer function extras_start_step(id)
           integer, intent(in) :: id
           integer :: ierr
           type (star_info), pointer :: s
           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return
           extras_start_step = 0
      end function extras_start_step

     


      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going

         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'
      
          write(*,*) 'From extras_check_model', Orbital_separation

         if (Orbital_separation <= s% r(1) + s% x_ctrl(2) * Rsun .and. s% use_other_energy ) then
            s% R_function2_param1 = Orbital_separation/(s%r(1)/Rsun) + 2.0 * 1.0 * s% x_ctrl(2) * Rsun/s%r(1)
            s% R_function2_param2 = Orbital_separation/(s%r(1)/Rsun) - 2.0 * 1.0 * s% x_ctrl(2) * Rsun/s%r(1)
            write(*,*) 'From extras_check_model', s% R_function2_param1/Rsun, s% R_function2_param2/Rsun
          endif

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 16
      end function how_many_extra_history_columns


      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         use num_lib
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         names(1) = 'Orbital_separation' ! Radial distance from stellar center of engulfed planet
         names(2) = 'Orbital_velocity' ! v_kepler
         names(3) = 'Log_Infall_distance'  ! dr
         names(4) = 'Log_Injected_energy'  ! de=dl*dt
         names(5) = 'Log_Destruction_factor'  ! Eq. 5 from Jia & Spruit 2018
         names(6) = 'Engulfed_fraction'  ! Cross section of the planet/Bondi area engulfed in the star (plane parallel approx, i.e. Max(Rplanet, Rbondi) << Rstar)
         names(7) = 'Total_mass_affected'
         names(8) = 'Planet_mass'
         names(9) = 'Planet_radius'
         names(10) = 'Bondi_radius'
         names(11) = 'Sound_speed'
         names(12) = 'Tidal_timescale' ! In years
         names(13) = 'Log_Infall_distance_tides' ! dr_tides
         names(14) = 'Infall_velocity' ! v_r [kms]
         names(15) = 'Enclosed_mass' ! msun 
         names(16) = 'Local_density' ! cgs 
         vals(1) = Orbital_separation / Rsun
         vals(2) = s% xtra(1)                 ! Orbital velocity
         vals(3) = safe_log10( s% xtra(2)) ! Infall distance
         vals(4) = safe_log10( s% xtra(3)) ! Injected energy
         vals(5) = safe_log10( s% xtra(4)) ! Disruption factor (ratio between ram pressure and binding energy density)
         vals(6) = s% xtra(5)                 ! Engulfed fraction
         vals(7) = s% xtra(6)
         vals(8) = s% x_ctrl(1)
         vals(9) = s% x_ctrl(2)
         vals(10) = s% xtra(7)
         vals(11) = s% xtra(8)
         vals(12) = s% xtra(9)
         vals(13) = safe_log10( s% xtra(10))
         vals(14) = s% xtra(11)
         vals(15) = s% xtra(12)
         vals(16) = s% xtra(13)
         ! note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in log column options.
         ! it must not include the new column names you are adding here.
      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id)
         use star_def, only: star_info
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 1
      end function how_many_extra_profile_columns


      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         use num_lib
         use math_lib
     !    use auto_diff
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! Adding extra heating so that we can plot on pgstar
         names(1) = 'engulfment_heating'
         do k = 1, nz
            vals(k,1) =  safe_log10( s% extra_heat(k)%val )  
         end do

         !note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.
      end subroutine data_for_extra_profile_columns


       integer function how_many_extra_history_header_items(id)
           integer, intent(in) :: id
           integer :: ierr
           type (star_info), pointer :: s
           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return
           how_many_extra_history_header_items = 0
        end function how_many_extra_history_header_items


        subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
           integer, intent(in) :: id, n
           character (len=maxlen_history_column_name) :: names(n)
           real(dp) :: vals(n)
           type(star_info), pointer :: s
           integer, intent(out) :: ierr
           ierr = 0
           call star_ptr(id,s,ierr)
           if(ierr/=0) return

           ! here is an example for adding an extra history header item
           ! also set how_many_extra_history_header_items
           ! names(1) = 'mixing_length_alpha'
           ! vals(1) = s% mixing_length_alpha

        end subroutine data_for_extra_history_header_items



  

       integer function how_many_extra_profile_header_items(id)
           integer, intent(in) :: id
           integer :: ierr
           type (star_info), pointer :: s
           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return
           how_many_extra_profile_header_items = 0
        end function how_many_extra_profile_header_items




      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
           integer, intent(in) :: id, n
           character (len=maxlen_profile_column_name) :: names(n)
           real(dp) :: vals(n)
           type(star_info), pointer :: s
           integer, intent(out) :: ierr
           ierr = 0
           call star_ptr(id,s,ierr)
           if(ierr/=0) return

           ! here is an example for adding an extra profile header item
           ! also set how_many_extra_profile_header_items
           ! names(1) = 'mixing_length_alpha'
           ! vals(1) = s% mixing_length_alpha

        end subroutine data_for_extra_profile_header_items



      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.

      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr, k
         logical :: grazing_phase
         real(dp) :: dr, delta_e, area, energy, R_influence, penetration_depth, dr_next
         real(dp) :: v_kepler
         character (len=90) :: fmt
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         call store_extra_info(s)

       ! Update the radial coordinate of the engulfed companion
       if ((.not. s% doing_first_model_of_run) .and. s% use_other_energy ) then
         Orbital_separation = max(Orbital_separation-Deltar-Deltar_tides, 0d0)
       end if


       ! Decrease timestep to fill da_tides/a tolerance
       ! write(*,*) 'Tff: ', s% xtra9, 'Tolerance: ',s% x_ctrl(7) , 'Orbital Sep: ',Orbital_separation
       !write(*,*) 'Deltar_tides/Orbital_separation',Deltar_tides/Orbital_separation,'Tolerance', s% x_ctrl(7)
       !write(*,*) 'Dt, Dt_next', s% dt, s% dt_next
       ! write(*,*) ' Deltar_tides/Orbital_separation > tolerance (Bad if > 1)', Deltar_tides/Orbital_separation, s% x_ctrl(7)
       if ((Deltar_tides/Orbital_separation >= s% x_ctrl(7)) .and. s% xtra(9) > 0.0 .and. s% use_other_energy) then
          s% dt_next = s% x_ctrl(7)* secyer * s% xtra(9)  ! Dt = tolerance * t_tide
          write(*,*) 'TIDES SETTING NEXT TIMESTEP: ', s% x_ctrl(7), s% xtra(9), s% dt_next
       end if
       ! CALCULATE R_influence AGAIN as it is not available to this part of the code
       ! But before we need bondi radius and orbital velocity.
       ! Calculate approx gridpoint location of planet center


         k=1
         do while (s% r(k) > Orbital_separation)
            k=k+1
         end do

         call orbital_velocity(s% m(1), Orbital_separation, v_kepler)
         call bondi_radius (s% x_ctrl(1)*Msun, s% csound(k), v_kepler, R_bondi)
         R_influence = max(R_bondi,s% x_ctrl(2)*Rsun)
        ! write(*,*)'From early',Orbital_separation/Rsun, Deltar/Rsun, R_influence/Rsun,R_bondi/Rsun

       ! Stop the run if: 1) we are at or past the stop age, but only if it has been set (default value is -101d0)
       !                  2) Orbital_separation is smaller than inlist-provided stop point x_ctrl(3) in Rsun or if
       !                  2B) Orbital_separation has become smaller than the companion radius or Bondi Radius, whichever is largest
         if (stop_age .GT. 0d0 .AND. s% star_age .GT. stop_age) then
           write(*,*) "Star should be thermally relaxed. Stopping."
           extras_finish_step = terminate
         endif
         if (Orbital_separation < R_influence .or. Orbital_separation < s% x_ctrl(3)*Rsun) then
          ! write(*,'(A,f10.4,A,f10.4,A,f10.4)') "Reached stop point from inlist. Orbital_separation=", Orbital_separation/Rsun, &
          !                                      'R_influence=',R_influence/Rsun, 'inslist stop:',s% x_ctrl(3)
           !extras_finish_step = terminate
           s% use_other_energy = .false. ! This also stops the tidal orbital evolution
           stop_age = s% star_age + 1d1 * s% kh_timescale
         endif

       ! Determine next timestep dt_next so that companion  infall distance dr is not too large compared to the influence radius
         delta_e = 0d0
         dr = 0d0
         energy = 0d0
         penetration_depth = 0d0
         grazing_phase = .false.

       ! Only do this if the planet is still around and falling into the star
         !write(*,*)'grazer outside if f_disruption, Deltar',s% model_number,s% xtra4,Deltar/Rsun
         if (s% xtra(4) < 1d0 .and. Deltar >= 0d0) then
         ! Calculate predicted dr in two cases: grazing phase and full engulfment
           penetration_depth = penetration_depth_function(R_influence,s% r(1), Orbital_separation)
           !write(*,*)'From outside grazer',s% model_number,R_influence/Rsun,(R_influence + s% r(1) - Orbital_separation)/Rsun,&
          !               penetration_depth/Rsun, R_influence/Rsun,s% r(1)/Rsun,Orbital_separation/Rsun
           if (penetration_depth <= 2d0*R_influence) then
              area = intercepted_area (penetration_depth, R_influence)
              grazing_phase = .true.
          !  write(*,*)'From grazer 1: r_infl, Orbital_separation-r_infl,Rstar,penetration,Deltar', &
          !              s% model_number,R_influence/Rsun,(Orbital_separation-R_influence)/Rsun,s% r(1)/Rsun, &
          !              penetration_depth/Rsun, Deltar/Rsun
           else
              area = pi * pow(R_influence, 2d0)
            !  write(*,*)'From grazer 2: r_infl, Orbital_separation-r_infl,Rstar,penetration,Deltar', &
            !            s% model_number,R_influence/Rsun,(Orbital_separation-R_influence)/Rsun,s% r(1)/Rsun, &
            !            penetration_depth/Rsun, Deltar/Rsun
           end if

         ! Estimate dr
           !write(*,*)'From inside finish step, just before drag call',s% model_number,s% m(k)/Msun
           call drag(s% m(k), s% x_ctrl(1)*Msun, area, s% rho(k), s% dt_next, s% r(k), delta_e, dr_next)
           !write(*,'(A,i4,f12.4,6e12.4)')'From end drag', &
          !      s% model_number,s% m(k)/Msun,area,s% rho(k),s% dt_next,s% r(k)/Rsun,delta_e,dr_next/Rsun

           if (grazing_phase) then                       ! Grazing Phase (requires small dr)
             do while (dr_next/R_influence > s% x_ctrl(4))
               s% dt_next = s% dt_next/2d0               ! There are better strategies, but this is simple enough
               call drag (s% m(k), s% x_ctrl(1)*Msun,area,s% rho(k),s% dt_next,s% r(k),delta_e,dr_next)
               write(*,*) s% model_number,&
                         'GRAZING ENGULFMENT SETTING DTNEXT: dr/R_influence too large: ' &
                         , dr_next/R_influence,'Decreasing dt to ', s% dt_next
             end do
           else
             do while (dr_next/R_influence > s% x_ctrl(5))   ! or Full Engulfment (allow for larger dr)
               s% dt_next = s% dt_next/2d0
               call drag (s% m(k), s% x_ctrl(1)*Msun,area,s% rho(k),s% dt_next,s% r(k),delta_e,dr_next)
               write(*,*) s% model_number,&
                         ! 'Engulfed dr/r_p too large: ', dr/(max(R_bondi,s% x_ctrl(2) * Rsun)),'Decreasing timestep to ', s% dt_next
                         'ENGULFMENT SETTING DTNEXT: dr/R_influence too large: ' &
                         , dr_next/R_influence,'Decreasing dt to ', s% dt_next
             end do
           end if
         end if
        ! #################### PRINT OUT #######################################
         fmt = '(f11.2,f11.3,f11.6,f11.6,f11.3,f11.3,f11.3,f11.5,f11.5,f11.6)'
         fmt=trim(fmt)

         write(*,'(a)') &
            '_______________________________________________________________________' // &
            '___________________________________________________________________________'
         write(*,*)
         write(*,'(a)') &
            '   a (rsun)      t_tide    Dr_drag   Dr_tides  &
             R_influence   R_Bondi   R_Star    Dr_next    Dt_next   v_infall [km/s]'

        write(*,fmt=fmt) Orbital_separation/Rsun, s% xtra(9), Deltar/Rsun,&
         Deltar_tides/Rsun, R_influence/Rsun,R_bondi/Rsun,s% r(1)/Rsun, dr_next/Rsun, s% dt_next/Rsun, dr_next/s% dt_next/1d5
         write(*,'(a)') &
            '_______________________________________________________________________' // &
            '___________________________________________________________________________'

         ! to save a profile,
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step



      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve


      ! routines for saving and restoring extra data so can do restarts

         ! put these defs at the top and delete from the following routines
         !integer, parameter :: extra_info_alloc = 1
         !integer, parameter :: extra_info_get = 2
         !integer, parameter :: extra_info_put = 3


      subroutine alloc_extra_info(s)
         integer, parameter :: extra_info_alloc = 1
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_alloc)
      end subroutine alloc_extra_info


      subroutine unpack_extra_info(s)
         integer, parameter :: extra_info_get = 2
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_get)
      end subroutine unpack_extra_info


      subroutine store_extra_info(s)
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_put)
      end subroutine store_extra_info


      subroutine move_extra_info(s,op)
         integer, parameter :: extra_info_alloc = 1
         integer, parameter :: extra_info_get = 2
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         integer, intent(in) :: op

         integer :: i, j, num_ints, num_dbls, ierr

         i = 0
         ! call move_int or move_flg
         num_ints = i

         i = 0
         ! call move_dbl
         ! (TAs) Important to understand what this is and what is done here. This is essential for MESA to remember Orbital_separation between timesteps
         ! and to allow for restarts from photos
         call move_dbl(Orbital_separation)
         call move_dbl(stop_age)
         num_dbls = i

         if (op /= extra_info_alloc) return
         if (num_ints == 0 .and. num_dbls == 0) return

         ierr = 0
         call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_alloc_extras'
            write(*,*) 'alloc_extras num_ints', num_ints
            write(*,*) 'alloc_extras num_dbls', num_dbls
            stop 1
         end if

         contains

         subroutine move_dbl(dbl)
            real(dp) :: dbl
            i = i+1
            select case (op)
            case (extra_info_get)
               dbl = s% extra_work(i)
            case (extra_info_put)
               s% extra_work(i) = dbl
            end select
         end subroutine move_dbl

         subroutine move_int(int)
            integer :: int
            i = i+1
            select case (op)
            case (extra_info_get)
               int = s% extra_iwork(i)
            case (extra_info_put)
               s% extra_iwork(i) = int
            end select
         end subroutine move_int

         subroutine move_flg(flg)
            logical :: flg
            i = i+1
            select case (op)
            case (extra_info_get)
               flg = (s% extra_iwork(i) /= 0)
            case (extra_info_put)
               if (flg) then
                  s% extra_iwork(i) = 1
               else
                  s% extra_iwork(i) = 0
               end if
            end select
         end subroutine move_flg

      end subroutine move_extra_info

      end module run_star_extras

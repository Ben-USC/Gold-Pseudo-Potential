

subroutine read_inpufile
   use common_variables
   integer :: iunit  = 10
   integer :: nclshl =  0
   integer :: ios, istat
   integer :: i
   character(20) :: filename = 'in7.dat'
   character(10) :: word, clshl

   call set_default

   open( unit=iunit, file=filename,  &
   &     iostat=ios, status='old', action='read' )
   if( ios /= 0 ) stop 'Error opening file input-file'
   print*, 'open :', filename

   preparedo: do
      read( unit=iunit, fmt="(a10)", iostat=istat ) word

      preblockif: if( istat == 0 ) then

         presectionif: if( word == '*pseudo-po' ) then
            presectiondo: do
            read( unit=iunit, fmt="(a10)", iostat=istat ) word

               presubif: if( word == '(pseudo-po' ) then

                  read( unit=iunit, fmt=*, iostat=istat ) methodpp
                  if( istat /= 0 ) stop "Read error-0000 in read_inputfile"

                  ! ----- NCPP
                  if(     methodpp == 1 ) then
                      lncpp = .true.
                  ! ----- USPP
                  elseif( methodpp == 2 ) then
                      luspp = .true.
                  ! ----- PAW
                  elseif( methodpp == 3 ) then
                      lpaw  = .true.
                  else
                      stop 'Pseudo-potential is not supported.'
                  end if

               elseif( word == '(closed sh' ) then presubif
                  read( unit=iunit, fmt=*, iostat=istat ) clshl
                  if( istat /= 0 ) stop "Read error-0001 in read_inputfile"

                  if( clshl/='[He]' .and. clshl/='[Ne]' .and. clshl/='[Ar]' .and.  &
                  &   clshl/='[Kr]' .and. clshl/='[Xe]' .and. clshl/='[Rn]' )  &
                  &   stop 'Read error-0002 in read_inputfile'

                  do
                     nljc(1)  = 100 ;  wnlj(1)  =  2.d0 !  1s
                     nclshl =  1
                     if( clshl == '[He]' ) exit
                     nljc(2)  = 200 ;  wnlj(2)  =  2.d0 !  2s
                     nljc(3)  = 210 ;  wnlj(3)  =  6.d0 !  2p
                     nclshl =  3
                     if( clshl == '[Ne]' ) exit
                     nljc(4)  = 300 ;  wnlj(4)  =  2.d0 !  3s
                     nljc(5)  = 310 ;  wnlj(5)  =  6.d0 !  3p
                     nclshl =  5
                     if( clshl == '[Ar]' ) exit
                     nljc(6)  = 320 ;  wnlj(6)  = 10.d0 !  3d
                     nljc(7)  = 400 ;  wnlj(7)  =  2.d0 !  4s
                     nljc(8)  = 410 ;  wnlj(8)  =  6.d0 !  4p
                     nclshl =  8
                     if( clshl == '[Kr]' ) exit
                     nljc(9)  = 420 ;  wnlj(9)  = 10.d0 !  4d
                     nljc(10) = 500 ;  wnlj(10) =  2.d0 !  5s
                     nljc(11) = 510 ;  wnlj(11) =  6.d0 !  5p
                     nclshl = 11
                     if( clshl == '[Xe]' ) exit
                     nljc(12) = 430 ;  wnlj(12) = 14.d0 !  4f
                     nljc(13) = 520 ;  wnlj(13) = 10.d0 !  5d
                     nljc(14) = 600 ;  wnlj(14) =  2.d0 !  6s
                     nljc(15) = 610 ;  wnlj(15) =  6.d0 !  6p
                     nclshl = 15
                     ! if( clshl == '[Rn]' ) exit
                     exit
                  end do

               elseif( word == '*(end)    ' ) then presubif
                  exit presectiondo
               end if presubif
            end do presectiondo
         end if presectionif

      elseif( istat < 0 ) then preblockif

            exit preparedo

      else preblockif

            stop 'error-1 in read_inputfile'

      end if preblockif
   end do preparedo

   if( methodpp == 0 ) stop 'Does not read section: (pseudo-potential)'

   rewind(iunit)

   readdo: do
      read( unit=iunit, fmt="(a10)", iostat=istat ) word

      blockif: if( istat == 0 ) then

         sectonif: if( word == '*pseudo-po' ) then
            sectiondo: do
               read( unit=iunit, fmt="(a10)", iostat=istat ) word

               subif: if( word == '(atomic nu' ) then

                  read( unit=iunit, fmt=*, iostat=istat ) zatm
                  if( istat /= 0 ) stop "Read error-0003 in read_inputfile"
                  read( unit=iunit, fmt=*, iostat=istat ) xion
                  if( istat /= 0 ) stop "Read error-0004 in read_inputfile"

               elseif( word == '(starting ' ) then subif

                  read( unit=iunit, fmt=*, iostat=istat ) lsg
                  if( istat /= 0 ) stop "Read error-0005 in read_inputfile"

               elseif( word == '(print AE ' ) then subif

                  read( unit=iunit, fmt=*, iostat=istat ) lprn
                  if( istat /= 0 ) stop "Read error-0006 in read_inputfile"
                  read( unit=iunit, fmt=*, iostat=istat ) lplot
                  if( istat /= 0 ) stop "Read error-0007 in read_inputfile"

               elseif( word == '(latter co' ) then subif

                  read( unit=iunit, fmt=*, iostat=istat ) latt
                  if( istat /= 0 ) stop "Read error-0008 in read_inputfile"

               elseif( word == '(configura' ) then subif

                  read( unit=iunit, fmt=*, iostat=istat ) nshl
                  if( istat /= 0 ) stop "Read error-0009 in read_inputfile"

                  nshl = nshl + nclshl
                  do i = nclshl+1, nshl, 1
                     read( unit=iunit, fmt=*, iostat=istat ) nljc(i), wnlj(i)
                     if( istat /= 0 ) stop "Read error-0010 in read_inputfile"
                  end do

               elseif( word == '(accuracy)' ) then subif

                  read( unit=iunit, fmt=*, iostat=istat ) edel
                  if( istat /= 0 ) stop "Read error-0011 in read_inputfile"
                  read( unit=iunit, fmt=*, iostat=istat ) vdel
                  if( istat /= 0 ) stop "Read error-0012 in read_inputfile"

               elseif( word == '(charge mi' ) then subif

                  read( unit=iunit, fmt=*, iostat=istat ) phi
                  if( istat /= 0 ) stop "Read error-0013 in read_inputfile"

               elseif( word == '(radial me' ) then subif

                  read( unit=iunit, fmt=*, iostat=istat ) mesh1
                  if( istat /= 0 ) stop "Read error-0014 in read_inputfile"
                  read( unit=iunit, fmt=*, iostat=istat ) xh
                  if( istat /= 0 ) stop "Read error-0015 in read_inputfile"
                  read( unit=iunit, fmt=*, iostat=istat ) rmax1
                  if( istat /= 0 ) stop "Read error-0016 in read_inputfile"

!                elseif( word == '(pseudo-po' ) then subif

!                   read( unit=iunit, fmt=*, iostat=istat ) methodpp
!                   if( istat /= 0 ) stop "Read error-0017 in read_inputfile"

               elseif( word == '(valence) ' ) then subif

                  read( unit=iunit, fmt=*, iostat=istat ) nval
                  if( istat /= 0 ) stop "Read error-0018 in read_inputfile"

                  do i = 1, nval, 1
                     read( unit=iunit, fmt=*, iostat=istat ) iref(i)
                     if( istat /= 0 ) stop "Read error-0019 in read_inputfile"

                     if( iref(i) == 1 ) then
                         read( unit=iunit, fmt=*, iostat=istat )         rus(i), rnc(i)
                         if( istat /= 0 ) stop "Read error-0020 in read_inputfile"
                     elseif( iref(i) == 2 ) then
                         read( unit=iunit, fmt=*, iostat=istat ) ref(i), rus(i), rnc(i)
                         if( istat /= 0 ) stop "Read error-0021 in read_inputfile"
                     end if
                  end do

               elseif( word == '(local pot' ) then subif

                  read( unit=iunit, fmt=*, iostat=istat ) v0
                  if( istat /= 0 ) stop "Read error-0023 in read_inputfile"
                  read( unit=iunit, fmt=*, iostat=istat ) rloc
                  if( istat /= 0 ) stop "Read error-0024 in read_inputfile"

               elseif( word == '(transfera' ) then subif

                  read( unit=iunit, fmt=*, iostat=istat ) rtra
                  if( istat /= 0 ) stop "Read error-0022 in read_inputfile"

               elseif( word == '(G-functio' ) then subif

                  if( lpaw ) then
                      read( unit=iunit, fmt=*, iostat=istat ) fcomp
                      if( istat /= 0 ) stop "Read error-0025 in read_inputfile"
                  end if

               elseif( word == '*(end)    ' ) then subif
                  exit
               end if subif

            end do sectiondo
         end if sectonif
      elseif( istat < 0 ) then blockif
         ! ----- EOF
         exit
      else blockif
         stop 'error-2 in read_inputfile'
      end if blockif
   end do readdo

   close( unit=iunit )

   ! ----- error trap
   if( zatm == 0.d0 ) stop 'Does not read section: (atomic number).'
   if( nshl == 0    ) stop 'Does not read section: (configuration).'
   if( nval == 0    ) stop 'Does not read section: (valence).'

end subroutine read_inpufile



subroutine set_default
   use common_variables

   ! No default value. These should be read
   zatm = 0.d0
   nshl = 0
   nval = 0
   methodpp = 0
   lncpp = .false.
   luspp = .false.
   lpaw  = .false.

   ! default value
   xion   = 0.d0
   lsg    = .false.
   lprn   = .true.
   lplot  = .true.
   latt   = .false.
   edel   = 5.0E-09
   vdel   = 1.0E-06
   phi    = 0.3d0
   mesh1  = 2001
   xh     = 100.d0
   rmax1  = 100.d0
   rtra   = 5.d0
   v0     = -2.5d0
   rloc   = 1.8d0
   fcomp  = 1.5d0

end subroutine set_default


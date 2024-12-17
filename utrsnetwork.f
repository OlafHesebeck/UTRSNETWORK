c     Olaf Hesebeck, Fraunhofer IFAM (2024):
c     This user subroutine implements a tabellaric shift function
c     similar to *TRS, tabular.
c     The shift is limited in order to avoid too small relaxation times.

c     Properties to be defined in the input file:
c     The first property given is the minimum admissible shift factor "a".
c     It is followed by an arbitrary number of pairs of numbers,
c     which consist of the base 10 logarithm of the shift function "log10(a(T))"
c     and the corresponding temperature "T", just like the table given to
c     *TRS, tabular.


      subroutine utrsnetwork (
C Must be updated
     *   outputData,
C Can be updated
     *   statev,
C Information (Read only)
     *   nOutput,
     *   nstatv,
     *   networkid,
     *   coords,
     *   temp,
     *   dtemp,
     *   nfield,
     *   predef,
     *   dpred,
     *   nprops,
     *   props,
     *   i_array,
     *   niarray,
     *   r_array,
     *   nrarray,
     *   c_array,
     *   ncarray)
C
      include 'aba_param.inc'
C
      parameter( io_trs_shift_begin = 1,
     *           io_trs_shift_end   = 2 )
C
      parameter( i_trs_kstep   = 1,
     *           i_trs_kinc    = 2,
     *           i_trs_noel    = 3,
     *           i_trs_npt     = 4,
     *           i_trs_layer   = 5,
     *           i_trs_kspt    = 6 )
C
      parameter( ir_trs_step_time  = 1,
     *           ir_trs_total_time = 2,
     *           ir_trs_creep_time = 3,
     *           ir_trs_timeinc    = 4 )

C
      parameter( ic_trs_material_name = 1 )
C
      dimension 
     *   statev(nstatv),
     *   predef(nfield),
     *   dpred(nfield),
     *   props(nprops),
     *   coords(*),
     *   outputData(nOutput),
     *   i_array(niarray),
     *   r_array(nrarray)

      character*80 c_array(ncarray)

      parameter( dln10=2.30258509299d0)
      
      shiftMin = props(1)
      numT = (nprops - 1) / 2

      temperatur = temp - dtemp

      if (temperatur .le. props(3)) then
         aLog = props(2)
      elseif (temperatur .ge. props(nprops)) then
          aLog = props(nprops-1)
      else
         do i=1, numT-1
            t1 = props(2*i+1)
            t2 = props(2*i+3)
            if (t1 .le. temperatur .and. t2 .gt. temperatur) then
               aLog1 = props(2*i)
               aLog2 = props(2*i+2)
               aLog = aLog1 + (aLog2-aLog1) * (temperatur-t1) / (t2-t1)
            end if
         end do
      endif
      a = exp(dln10 * aLog)

      outputData(io_trs_shift_begin) = max(shiftMin, a)

      temperatur = temp

      if (temperatur .le. props(3)) then
         aLog = props(2)
      elseif (temperatur .ge. props(nprops)) then
          aLog = props(nprops-1)
      else
         do i=1, numT-1
            t1 = props(2*i+1)
            t2 = props(2*i+3)
            if (t1 .le. temperatur .and. t2 .gt. temperatur) then
               aLog1 = props(2*i)
               aLog2 = props(2*i+2)
               aLog = aLog1 + (aLog2-aLog1) * (temperatur-t1) / (t2-t1)
            end if
         end do
      endif
      a = exp(dln10 * aLog)

      outputData(io_trs_shift_end) = max(shiftMin, a)

      return
      end

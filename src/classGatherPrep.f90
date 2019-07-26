!> \file
!! Assigns values to pointer vectors relating the location of elements on the "gathered" variable
!! vectors to elements on the original two-dimensional arrays (latitude circle x mosaic tiles) for land grid
!! cells.
!! @author D. Verseghy, M. Lazare, E. Chan
!
subroutine classGatherPrep(ILMOS,JLMOS,IWMOS,JWMOS, & ! Formerly GATPREP
                           NML,NMW,GCROW,FAREA,MOSID, &
                           NL,NM,ILG,IL1,IL2,IM)
  !
  !     * JAN 12/17 - D.VERSEGHY. NOTE: THIS VERSION OF classGatherPrep
  !     *                         IS DESIGNED SPECIFICALLY FOR LAND.
  !     *                         THE VERSION CURRENTLY USED IN THE AGCM
  !     *                         (SINCE JAN. 2014) HAS BEEN GENERALIZED
  !     *                         FOR LAND, LAKES AND WATER/ICE.
  !     * DEC 28/11 - D.VERSEGHY. CHANGE ILGM BACK TO ILG AND
  !     *                         ILG TO NL FOR CONSISTENCY WITH
  !     *                         BOTH STAND-ALONE AND GCM
  !     *                         CONVENTIONS.
  !     * OCT 22/11 - M.LAZARE. REMOVE OCEAN/ICE CODE (NOW DONE
  !     *                       IN COISS).
  !     * OCT 21/11 - M.LAZARE. COSMETIC: ILG->ILGM AND NLAT->ILG,
  !     *                       TO BE CONSISTENT WITH MODEL
  !     *                       CONVENTION. ALSO GCGRD->GCROW.
  !     * JUN 12/06 - E.CHAN.  DIMENSION IWAT AND IICE BY ILG.
  !     * NOV 03/04 - D.VERSEGHY. ADD "IMPLICIT NONE" COMMAND.
  !     * AUG 09/02 - D.VERSEGHY/M.LAZARE. DETERMINE INDICES FOR
  !     *                        GATHER-SCATTER OPERATIONS ON
  !     *                        CURRENT LATITUDE LOOP.
  !
  implicit none
  !
  !     * INTEGER CONSTANTS.
  !
  integer, intent(inout)  :: NML  !< Total number of mosaic tiles in land surface gather vectors
  integer, intent(inout)  :: NMW  !< Total number of mosaic tiles in inland water gather vectors
  integer, intent(in)  :: NL   !< Number of latitude grid cells in land surface scatter vectors
  integer, intent(in)  :: NM   !< Total number of CLASS mosaic tiles in land surface gather vectors
  integer, intent(in)  :: ILG  !< Total number of mosaic tiles per latitude grid cell in land surface scatter vector
  integer, intent(in)  :: IL1, IL2  !<
  integer, intent(in)  :: IM   !< Maximum number of mosaic tiles within the grid cells in the array under consideration
  integer  :: I, J
  !
  !     * OUTPUT FIELDS.
  !
  integer, intent(out)  :: ILMOS  (ILG) !< Index of grid cell corresponding to current element
  !< of gathered vector of land surface variables [ ]
  integer, intent(out)  :: JLMOS  (ILG) !< Index of mosaic tile corresponding to current element
  !< of gathered vector of land surface variables [ ]
  integer, intent(out)  :: IWMOS  (ILG) !< Index of grid cell corresponding to current element of gathered vector
  !< of inland water body variables [ ]
  integer, intent(out)  :: JWMOS  (ILG) !< Index of mosaic tile corresponding to current element of gathered vector
  !< of inland water body variables [ ]
  !
  !     * INPUT FIELDS.
  !
  real     :: GCROW (NL)    !< Real number identifier indicating whether the grid cell
  !< is land (-1.0), sea ice (+1.0), or ocean (0.0)
  real     :: FAREA (NL,NM) !< Fractional coverage of mosaic tile on grid cell [ ]
  !
  integer  :: MOSID (NL,NM) !< Mosaic tile type identifier (1 for land, 0 for inland water)
  !---------------------------------------------------------------------
  NML = 0
  NMW = 0
  !>
  !! A looping operation is performed over the latitude circle, or array of grid cells, under consideration. If
  !! the grid cell is a land one (GCROW = -1.0), an additional internal loop is performed over all the mosaic
  !! tiles present. For each mosaic tile, if its fractional coverage is greater than zero, then if the mosaic type
  !! identifier MOSID is equal to 1 (indicating land), the counter of total mosaic tiles in the land surface gather
  !! vectors, NML, is incremented by one, and the elements of the vectors ILMOS and JLMOS corresponding
  !! to NML are set to the indices of the current grid cell and mosaic tile respectively. If MOSID is equal to
  !! zero (indicating inland water), the counter of total mosaic tiles in the inland water gather vectors, NMW,
  !! is incremented by one, and the elements of the vectors IWMOS and JWMOS corresponding to NMW are
  !! set to the indices of the current grid cell and mosaic tile respectively.
  !!
  do I = IL1,IL2 ! loop 200
    if (GCROW(I) <= - 0.5) then
      do J = 1,IM ! loop 100
        if (FAREA(I,J) > 0.0) then
          if (MOSID(I,J) > 0) then
            NML = NML + 1
            ILMOS(NML) = I
            JLMOS(NML) = J
          else
            NMW = NMW + 1
            IWMOS(NMW) = I
            JWMOS(NMW) = J
          end if
        end if
      end do ! loop 100
    end if
  end do ! loop 200

  return
end

c
c Loriano Storchi this file is part of libdstcltmtx
c 
c libdstcltmtx is free software: you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation, either version 3 of the License, or
c (at your option) any later version.
c
c Foobar is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c
c You should have received a copy of the GNU General Public License
c along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
c    


      subroutine DISTRIBUTE (root, globalmtx, descmtx, ishift, jshift, & 
                             globalmtxrow, globalmtxcol, localmtx, istart, & 
                             jstart, rowdim, coldim, info)

      implicit none 

      include "mpif.h"

      integer :: CSRC_, CTXT_, DLEN_, DTYPE_, LLD_, MB_, M_, NB_, N_, &
                 RSRC_
      parameter (CSRC_ = 8, CTXT_ = 2, DLEN_ = 9, DTYPE_ = 1, &
                 LLD_ = 9, MB_ = 5, M_ = 3, NB_ = 6, N_ = 4, RSRC_ = 7)

      integer, external :: numroc, iceil, blacs_pnum, indxg2p
      logical, external :: findnode

      integer :: nbrow, nbcol, myrow, mycol, nprow, npcol, locr, &
                 locc, node, irow, icol, istatus(MPI_STATUS_SIZE), & 
                 localmtxcol, localmtxrow, nbrows, nbcolumns, ibcolumns, &
                 col_max, ibrows, row_max, columns, rows, indexi, &
                 indexj, lindexi, lindexj, root, context, rootrow, &
                 rootcol, noderow, nodecol, noderowidx, nodecolidx, &
                 mpierr, rank, numprocs, descmtx(*), istart, jstart, &
                 rowdim, coldim, info, i, j, k, indxg2pactual,  &
                 indxg2pold, il2g, jl2g, ishift, jshift, globalmtxrow, &
                 globalmtxcol

      integer, allocatable :: nodeset(:), nodesetrow(:), nodesetcol(:)
      
      TIPO :: localmtx(descmtx(LLD_ ), *), globalmtx(globalmtxrow, globalmtxcol)

      TIPO, allocatable :: tmpmtx(:,:)


      if (.not.((((descmtx(M_ ) - globalmtxrow) .eq. 0) .and. &
                ((descmtx(N_ ) - globalmtxcol) .eq. 0)) .or. &
               (((globalmtxrow - rowdim) .eq. 0) .and. &
                ((globalmtxcol - coldim) .eq. 0)))) then
        info = -1
        return 
      endif
   
      info = 0

      if ((istart.gt.0).and.(jstart.gt.0).and. &
          (rowdim.gt.0).and.(coldim.gt.0)) then
      if (((istart + rowdim - 1).gt.descmtx(M_ )).or. &
            ((jstart + coldim - 1).gt.descmtx(N_ ))) then
        info = 1
        return
      endif
        else
        info = 1
        return 
      endif

      context = descmtx(CTXT_ )
      nbrow = descmtx(MB_ )
      nbcol = descmtx(NB_ )
      rootrow = descmtx(RSRC_ )
      rootcol = descmtx(CSRC_ )

      call blacs_pinfo (rank, numprocs )
      call blacs_gridinfo( context, nprow, npcol, myrow, mycol )

      allocate (nodeset(numprocs), nodesetrow(nprow), &
                nodesetcol(npcol))
      nodeset = -1
      nodesetrow = -1
      nodesetcol = -1

      ! compute the node that must be involved in the communication
      j = 1
      indxg2pold = -1
      do i=istart, istart+rowdim-1
        indxg2pactual = mod (rootrow + (i - 1) / nbrow, nprow)
        if ((indxg2pactual.ne.indxg2pold) .and. &
            (.not.findnode(indxg2pactual, nodesetrow, nprow))) then
          nodesetrow(j) = indxg2pactual
          indxg2pold = indxg2pactual
          j = j + 1
        endif
      enddo

      i = 1
      indxg2pold = -1
      do j=jstart, jstart+coldim-1
        indxg2pactual = mod (rootcol + (j - 1) / nbcol, npcol)
        if ((indxg2pactual.ne.indxg2pold) .and. &
            (.not.findnode(indxg2pactual, nodesetcol, npcol))) then
          nodesetcol(i) = indxg2pactual
          indxg2pold = indxg2pactual
          i = i + 1
        endif
      enddo

      k = 1
      do i=1, nprow
        if (nodesetrow(i) .eq. -1) exit
        do j=1, npcol
          if (nodesetcol(j) .eq. -1) exit
          nodeset(k) = blacs_pnum (context, nodesetrow(i),  &
                                   nodesetcol(j))
          k = k + 1
        enddo
      enddo

      localmtxrow = numroc( descmtx(M_ ), nbrow, myrow, rootrow, nprow)
      localmtxcol = numroc( descmtx(N_ ), nbcol, mycol, rootcol, npcol)

      if (rank.eq.root) then
        !node is the rank of the processor that sent the data
        !numprocs, is the number of processors in the communicator
        do noderowidx = 0, nprow-1
          noderow = mod (noderowidx + rootrow, nprow)
          do nodecolidx = 0, npcol-1
            nodecol = mod (nodecolidx + rootcol, npcol)

            node = blacs_pnum (context, noderow, nodecol)
            
            if (findnode(node, nodeset, numprocs)) then
              if (node.eq.root) then
                 irow = myrow
                 icol = mycol
                 !obtain dimensions of the local array on master node
                 locr = localmtxrow
                 locc = localmtxcol
              else
                 !status -- message status array(of type integer)
                 call MPI_Recv(irow, 1, MPI_INTEGER, node, 10, & 
                      MPI_COMM_WORLD, istatus, mpierr)
                 call MPI_Recv(icol, 1, MPI_INTEGER, node, 20, &
                      MPI_COMM_WORLD, istatus, mpierr)
              
                 CALL MPI_Recv(locr, 1, MPI_INTEGER, node, 40, &
                      MPI_COMM_WORLD, istatus, mpierr)
                 CALL MPI_Recv(locc, 1, MPI_INTEGER, node, 50, &
                      MPI_COMM_WORLD, istatus, mpierr)
                 
                 !create a local matrix with the size of the matrix passed
                 allocate(tmpmtx(locr,locc))
              endif
              
              !compute the number of blocks in each local array
              nbrows = ceiling(real(locr)/real(nbrow))    !number of blocks in the row
              nbcolumns = ceiling(real(locc)/real(nbcol))    !number of blocks in the columns
              
              !loop over each block in the column
              do ibcolumns = 1, nbcolumns  
                 !special case - number of columns is less than NB     
                 if (ibcolumns.eq.nbcolumns) then
                    col_max = locc-(nbcolumns-1)*nbcol 
                 else
                    col_max = nbcol !number of columns in a block        
                 endif
              
                 !loop over each block in the row
                 do ibrows = 1, nbrows
                    !special case - number of columns is less than NBCO
                    if (ibrows.eq.nbrows) then
                       row_max = locr-(nbrows-1)*nbrow
                    else
                       row_max = nbrow !number of rows in a block
                    endif
              
                    !for each column in the block, loop over each row in the block
                    do columns = 1, col_max
                      indexj = (nodecolidx+(ibcolumns-1)*npcol)* &
                                nbcol + columns
                      lindexj = nbcol*((indexj-1)/(nbcol*npcol))+ &
                                MOD(indexj-1,nbcol)+1

                      if ((indexj .ge. jstart) .and. (indexj .le. &
                          (jstart+coldim-1)))  then
                        do rows = 1, row_max
                          indexi = (noderowidx+(ibrows-1)*nprow)* &
                                    nbrow + rows
                          if ((indexi .ge. istart) .and. (indexi .le.  &
                               (istart+rowdim-1))) then
                        
                            lindexi = nbrow*((indexi-1)/(nbrow*nprow))+ &
                                      MOD(indexi-1,nbrow)+1

                            if (((descmtx(M_ ) - globalmtxrow) .eq. 0) .and.  &
                                ((descmtx(N_ ) - globalmtxcol) .eq. 0)) then
                              if (node.eq.root) then
                                localmtx (lindexi, lindexj) =  &
                                           globalmtx(indexi + ishift, & 
                                           indexj + jshift)
                              else
                                tmpmtx (lindexi, lindexj) =  &
                                            globalmtx(indexi + ishift, & 
                                            indexj + jshift)
                              endif
                            else if (((globalmtxrow - rowdim) .eq. 0) .and. &
                                     ((globalmtxcol - coldim) .eq. 0)) then
                              if (node.eq.root) then
                                localmtx (lindexi, lindexj) =  &
                                           globalmtx(indexi + ishift - istart, & 
                                           indexj + jshift - jstart)
                              else
                                tmpmtx (lindexi, lindexj) =  &
                                            globalmtx(indexi + ishift - istart, & 
                                            indexj + jshift - jstart)
                              endif
                            endif

                          endif
                        enddo
                      endif
                    enddo
                 enddo
              enddo
              
              if (node.ne.root) then 
                 call MPI_Send(tmpmtx, locr*locc, TYPE_MPI, &
                            node, 30, MPI_COMM_WORLD, mpierr)
                 deallocate (tmpmtx)
              endif
            endif
          enddo
        enddo
      else 
        if (findnode(rank, nodeset, numprocs)) then
          !send the processor coordinates in grid
          call MPI_Send(myrow, 1, MPI_INTEGER, root, 10, MPI_COMM_WORLD,  &
                    mpierr)
          call MPI_Send(mycol, 1, MPI_INTEGER, root, 20, MPI_COMM_WORLD,  &
                    mpierr)
          !send number rows and columns
          call MPI_Send(localmtxrow, 1, MPI_INTEGER, root, 40, &
                    MPI_COMM_WORLD, mpierr)
          call MPI_Send(localmtxcol, 1, MPI_INTEGER, root, 50, &
                    MPI_COMM_WORLD, mpierr)

          if ((globalmtxrow .eq. rowdim) .and. &
              (globalmtxcol .eq. coldim) .and. &
              (istart .eq. 1) .and. (jstart .eq. 1)) then 
            call MPI_Recv(localmtx, localmtxrow*localmtxcol, &
                          TYPE_MPI, root, 30, &
                          MPI_COMM_WORLD, istatus, mpierr)
          else 
            allocate(tmpmtx(localmtxrow, localmtxcol))
            ! receive the local matrix
            call MPI_Recv(tmpmtx, localmtxrow*localmtxcol, &
                   TYPE_MPI, root, 30, &
                   MPI_COMM_WORLD, istatus, mpierr)
            
            do i=1, localmtxrow
              il2g = nprow * nbrow * ((i - 1)/nbrow) + &
                     mod(i-1, nbrow) + mod(nprow+myrow-rootrow,  &
                     nprow) * nbrow + 1
              if ((il2g .ge. istart) .and. (il2g .le. (istart+rowdim-1))) &
                 then
                do j=1, localmtxcol
                  jl2g = npcol * nbcol * ((j - 1)/nbcol) + &
                         mod(j-1, nbcol) + mod(npcol+mycol-rootcol, &
                         npcol) * nbcol + 1
                  if ((jl2g .ge. jstart) .and. (jl2g .le. (jstart+ &
                      coldim-1))) then
                    localmtx(i, j) = tmpmtx(i, j)
                  endif
                enddo
              endif
            enddo
            deallocate(tmpmtx)
          endif

        endif
      endif

      deallocate (nodeset, nodesetcol, nodesetrow)

      call MPI_Barrier (MPI_COMM_WORLD, mpierr)

      return 
      end subroutine DISTRIBUTE

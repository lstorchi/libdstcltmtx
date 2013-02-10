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

      function findnode (node, nodeset, dim)

      ! migliorala, magari usando la find C++

      implicit none

      integer :: node, dim, i, nodeset(*)
      logical :: findnode

      findnode = .false.

      do i=1, dim
        if (nodeset(i) .eq. node) then
          findnode = .true.
          return
        else if (nodeset(i) .eq. -1) then
          return
        endif
      enddo

      return
      end function findnode

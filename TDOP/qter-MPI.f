      program qtetAtypAtyp_mpi
      use mpi
      implicit none

c     MPI variables
      integer ierr, rank, nprocs
      integer my_start, my_end, my_atoms
      integer status(MPI_STATUS_SIZE)

c     Original variables
      character *80 title,Atyp
      character CN*8,O*7
      integer nconf,nAtyp,iconf
      integer ii,j,k,l,num,n
      double precision fshell
      parameter(nconf = 5001)
      parameter(nAtyp = 87909)
      parameter(fshell = 0.175)
      double precision xbox,massAtyp,chrAtyp,ybox,zbox
      double precision Atypx(nAtyp)
      double precision Atypy(nAtyp),Atypz(nAtyp)
      double precision xi,yi,zi,swap
      double precision sumcos,ro_val
      double precision cosval,a,b,c,qtet(nAtyp)
      double precision mAtyp,chAtyp
      double precision atm(nAtyp),atid(nAtyp),ax(nAtyp)
      double precision by(nAtyp),cz(nAtyp)
      integer numAtyp(nAtyp),Ntot,ip,cAtyp
      integer Atypa

c     Optimized data structures
      double precision, allocatable :: neighbor_dist(:), neighbor_pos(:)
      integer, allocatable :: sorted_indices(:)
      double precision fshell_cutoff
      integer max_neighbors

c     Local arrays for MPI
      double precision, allocatable :: local_qtet(:)
      integer, allocatable :: local_cAtyp(:)
      integer, allocatable :: local_atom_ids(:)

c     Temporary arrays for receiving data
      double precision, allocatable :: recv_qtet(:)
      integer, allocatable :: recv_cAtyp(:)
      integer, allocatable :: recv_atom_ids(:)

c     File name variables
      character(len=50) :: output_filename

c     Initialize MPI
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

c     Pre-calculate cutoff
      fshell_cutoff = fshell * 2.0d0

c     Distribute atoms across processors
      my_atoms = nAtyp / nprocs
      my_start = rank * my_atoms + 1
      if (rank .eq. nprocs - 1) then
          my_end = nAtyp
      else
          my_end = (rank + 1) * my_atoms
      endif
      my_atoms = my_end - my_start + 1

      if (rank .eq. 0) then
          write(*,*) 'Running with', nprocs, 'MPI processes'
          write(*,*) 'Processing', nAtyp, 'atoms total'
          write(*,*) 'Processing', nconf, 'configurations'
          write(*,*) 'Process 0 handling atoms:', my_start, 'to', my_end
      endif

c     Allocate local arrays
      allocate(local_qtet(my_atoms))
      allocate(local_cAtyp(my_atoms))
      allocate(local_atom_ids(my_atoms))

c     Store local atom IDs
      do ii = 1, my_atoms
          local_atom_ids(ii) = my_start + ii - 1
      enddo

c     Estimate maximum possible neighbors for memory allocation
      max_neighbors = min(nAtyp, 2000)
      allocate(neighbor_dist(max_neighbors))
      allocate(neighbor_pos(max_neighbors))
      allocate(sorted_indices(max_neighbors))

c*************************************************
c******* MAIN LOOP OVER ALL FRAMES **************
c*************************************************
      if (rank .eq. 0) then
          open(unit=10,status='unknown',
     &     file='../OP-plumed/mutINP_ONC_OW1.gro')
      endif

      do 5000 iconf = 1, nconf
          if (rank .eq. 0) then
              write(*,*) 'Processing frame', iconf, 'of', nconf
          endif

c         File reading - only rank 0 reads
          if (rank .eq. 0) then
              read(10,*)
              read(10,*)
              do j=1,nAtyp
                  read(10,'(A8,A7,1i5,3f8.3)')CN,O,num,Atypx(j),
     &                 Atypy(j),Atypz(j)
              enddo
              read(10,*)xbox,ybox,zbox
          endif

c         Broadcast coordinate data to all processes
          call MPI_BCAST(Atypx, nAtyp, MPI_DOUBLE_PRECISION, 0, 
     &                   MPI_COMM_WORLD, ierr)
          call MPI_BCAST(Atypy, nAtyp, MPI_DOUBLE_PRECISION, 0, 
     &                   MPI_COMM_WORLD, ierr)
          call MPI_BCAST(Atypz, nAtyp, MPI_DOUBLE_PRECISION, 0, 
     &                   MPI_COMM_WORLD, ierr)
          call MPI_BCAST(xbox, 1, MPI_DOUBLE_PRECISION, 0, 
     &                   MPI_COMM_WORLD, ierr)
          call MPI_BCAST(ybox, 1, MPI_DOUBLE_PRECISION, 0, 
     &                   MPI_COMM_WORLD, ierr)
          call MPI_BCAST(zbox, 1, MPI_DOUBLE_PRECISION, 0, 
     &                   MPI_COMM_WORLD, ierr)

c         Synchronize all processes
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)

c*************************************************
c******* Main calculation loop (parallelized) ***
c*************************************************
          do 200 ii = 1, my_atoms
              j = local_atom_ids(ii)
              
c             Progress reporting
              if (rank .eq. 0 .and. mod(ii, 5000) .eq. 0) then
                  write(*,*) 'Frame', iconf, ': Processed', ii, 
     &                      'atoms on rank 0'
              endif
              
c             Initialize neighbor counter
              cAtyp = 0
              
c             Find all neighbors within cutoff distance
              do k=1,nAtyp
                  if(k .eq. j) cycle  ! Skip self
                  
                  xi=Atypx(j)-Atypx(k)
                  yi=Atypy(j)-Atypy(k)
                  zi=Atypz(j)-Atypz(k)
                  
c                 Apply periodic boundary conditions
                  xi=xi-xbox*dnint(xi/xbox)
                  yi=yi-ybox*dnint(yi/ybox)
                  zi=zi-zbox*dnint(zi/zbox)
                  
                  ro_val=dsqrt(xi**2 +yi**2 +zi**2)
                  
c                 Check if within cutoff
                  if(ro_val.lt.fshell_cutoff) then
                      cAtyp = cAtyp + 1
                      if(cAtyp .le. max_neighbors) then
                          neighbor_dist(cAtyp) = ro_val
                          neighbor_pos(cAtyp) = dble(k)
                      endif
                  endif
              enddo
              
c             Store neighbor count
              local_cAtyp(ii) = cAtyp
              
c             Ensure we have enough neighbors for tetrahedral calculation
              if(cAtyp .lt. 4) then
                  local_qtet(ii) = 0.0d0
                  goto 200
              endif
              
c             Limit to available memory
              if(cAtyp .gt. max_neighbors) then
                  cAtyp = max_neighbors
              endif

c*************************************************
c**** Optimized sorting using selection sort ****
c*************************************************
c             Initialize indices
              do k=1,cAtyp
                  sorted_indices(k) = k
              enddo
              
c             Selection sort
              do k=1,cAtyp-1
c                 Find minimum element
                  n = k
                  do l=k+1,cAtyp
                      if(neighbor_dist(sorted_indices(l)).lt.
     &                   neighbor_dist(sorted_indices(n))) then
                          n = l
                      endif
                  enddo
                  
c                 Swap indices if needed
                  if(n .ne. k) then
                      num = sorted_indices(k)
                      sorted_indices(k) = sorted_indices(n)
                      sorted_indices(n) = num
                  endif
              enddo

c*************************************************
c******* Calculating angles and qtet ************
c*************************************************
              sumcos=0.0d0
              
c             Use the 4 nearest neighbors for tetrahedral calculation
              do k=1,3
                  do l=k+1,4
c                     Get positions of neighbors
                      xi=Atypx(int(neighbor_pos(sorted_indices(k))))-
     &                   Atypx(int(neighbor_pos(sorted_indices(l))))
                      yi=Atypy(int(neighbor_pos(sorted_indices(k))))-
     &                   Atypy(int(neighbor_pos(sorted_indices(l))))
                      zi=Atypz(int(neighbor_pos(sorted_indices(k))))-
     &                   Atypz(int(neighbor_pos(sorted_indices(l))))
                      
c                     Apply periodic boundary conditions
                      xi=xi-xbox*dnint(xi/xbox)
                      yi=yi-ybox*dnint(yi/ybox)
                      zi=zi-zbox*dnint(zi/zbox)
                      
                      c = dsqrt(xi**2 +yi**2 +zi**2)
                      a = neighbor_dist(sorted_indices(l))
                      b = neighbor_dist(sorted_indices(k))
                      
c                     Calculate cosine with bounds checking
                      if(a.gt.1.0d-10 .and. b.gt.1.0d-10) then
                          cosval=(a*a + b*b - c*c)/(2.0d0*a*b)
c                         Clamp cosval to valid range [-1, 1]
                          cosval = max(-1.0d0, min(1.0d0, cosval))
                          sumcos=sumcos+(cosval+(1.0d0/3.0d0))**2.0d0
                      endif
                  enddo
              enddo
              
              local_qtet(ii)=(1.0d0-(3.0d0/8.0d0)*sumcos)
              
200       continue

c         Synchronize all processes
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)

c*************************************************
c******* Gather results for this frame **********
c*************************************************
          if (rank .eq. 0) then
              write(*,*) 'Gathering results for frame', iconf
          endif

c         Create output filename for this frame
          if (rank .eq. 0) then
              write(output_filename, '(A,I4.4,A)') 
     &              '0-20ps-TDOP-frame', iconf, '.dat'
              open(unit=20,file=output_filename,status="unknown")
              
c             Write own results first
              do ii = 1, my_atoms
                  j = local_atom_ids(ii)
                  write(20,'(I5,1x,F12.8,1x,I3)') j, local_qtet(ii), 
     &                  local_cAtyp(ii)
              enddo
              
c             Receive and write results from other processes
              do k = 1, nprocs - 1
c                 Calculate how many atoms this process handled
                  if (k .eq. nprocs - 1) then
                      n = nAtyp - k * (nAtyp / nprocs)
                  else
                      n = nAtyp / nprocs
                  endif
                  
c                 Allocate temporary arrays for receiving
                  allocate(recv_qtet(n))
                  allocate(recv_cAtyp(n))
                  allocate(recv_atom_ids(n))
                  
c                 Receive data
                  call MPI_RECV(recv_qtet, n, MPI_DOUBLE_PRECISION, k, 
     &                         0, MPI_COMM_WORLD, status, ierr)
                  call MPI_RECV(recv_cAtyp, n, MPI_INTEGER, k, 1,
     &                         MPI_COMM_WORLD, status, ierr)
                  call MPI_RECV(recv_atom_ids, n, MPI_INTEGER, k, 2,
     &                         MPI_COMM_WORLD, status, ierr)
                  
c                 Write received data
                  do ii = 1, n
                      write(20,'(I5,1x,F12.8,1x,I3)') recv_atom_ids(ii),
     &                      recv_qtet(ii), recv_cAtyp(ii)
                  enddo
                  
c                 Clean up temporary arrays
                  deallocate(recv_qtet)
                  deallocate(recv_cAtyp)
                  deallocate(recv_atom_ids)
              enddo
              
              close(20)
              write(*,*) 'Frame', iconf, 'results written to', 
     &                  output_filename
              
          else
c             Send results to rank 0
              call MPI_SEND(local_qtet, my_atoms, MPI_DOUBLE_PRECISION,
     &                     0, 0, MPI_COMM_WORLD, ierr)
              call MPI_SEND(local_cAtyp, my_atoms, MPI_INTEGER, 0, 1,
     &                     MPI_COMM_WORLD, ierr)
              call MPI_SEND(local_atom_ids, my_atoms, MPI_INTEGER, 0, 2,
     &                     MPI_COMM_WORLD, ierr)
          endif

c         Synchronize before next frame
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)

5000  continue

c     Close input file
      if (rank .eq. 0) then
          close(10)
      endif

c     Synchronize before cleanup
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

c     Cleanup - only deallocate what was allocated
      if (allocated(neighbor_dist)) deallocate(neighbor_dist)
      if (allocated(neighbor_pos)) deallocate(neighbor_pos)
      if (allocated(sorted_indices)) deallocate(sorted_indices)
      if (allocated(local_qtet)) deallocate(local_qtet)
      if (allocated(local_cAtyp)) deallocate(local_cAtyp)
      if (allocated(local_atom_ids)) deallocate(local_atom_ids)

c     Finalize MPI
      call MPI_FINALIZE(ierr)

      if (rank .eq. 0) then
       write(*,*) 'All', nconf, 'configurations processed successfully!'
          write(*,*) 'Output files: 0-20ps-TDOP-frame001.dat to',
     &              '0-20ps-TDOP-frame', nconf, '.dat'
      endif

      stop
      end

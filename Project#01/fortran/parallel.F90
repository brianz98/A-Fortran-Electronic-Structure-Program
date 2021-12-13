module parallel
    #ifdef PARALLEL
    use mpi_f08
    #endif

    use const

    implicit none

    ! MPI-related constants
    integer :: iproc, nprocs
    integer, parameter :: root = 0
    logical :: parent

    ! MPI communicators
    integer :: inter_node_comm, intra_node_comm

    ! OpenMP thread count
    integer :: nthreads

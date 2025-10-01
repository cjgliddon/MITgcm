subroutine floorice_init(myThid)

    implicit none
    integer, intent(in) :: myThid

    if (FLOORICEfreezingRateFile .NE. ' ' ) THEN
        call READ_FLD_XY_RL( FLOORICEfreezingRateFile, ' ', ice_dmdt, 0, myThid )

end subroutine floorice_init
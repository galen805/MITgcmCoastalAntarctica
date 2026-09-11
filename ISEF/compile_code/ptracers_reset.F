#include "PTRACERS_OPTIONS.h"
#include "GAD_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: PTRACERS_RESET

C     !INTERFACE:
      SUBROUTINE PTRACERS_RESET( myTime, myIter, myThid )

C     !DESCRIPTION:
C     Re-initialize PTRACERS if it is the correct time to do so

C     !USES:
#include "PTRACERS_MOD.h"
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"
#include "GAD.h"

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C added by Yilang Xu on April 8, 2020
C to link surface source of ptracers to sea ice area
#ifdef ALLOW_SEAICE
#include "SEAICE_SIZE.h"
#include "SEAICE_PARAMS.h"
#include "SEAICE.h"
#endif /* ALLOW_SEAICE */
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C     !INPUT PARAMETERS:
C     myThid               :: thread number
      INTEGER myIter
      _RL myTime
      INTEGER myThid

#ifdef ALLOW_PTRACERS

      LOGICAL  DIFFERENT_MULTIPLE
      EXTERNAL DIFFERENT_MULTIPLE

C     !LOCAL VARIABLES:
C     i,j,k,bi,bj,iTracer  :: loop indices
C     msgBuf      :: Informational/error message buffer
      INTEGER i,j,k,bi,bj,iTracer
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      CHARACTER*(MAX_LEN_FNAM) tmpInitialFile

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C added by Yilang Xu on April 8, 2020
C to link surface source of ptracers to sea ice area
#ifdef ALLOW_SEAICE
      INTEGER ig,jg
#endif /* ALLOW_SEAICE */
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef PTRACERS_ALLOW_DYN_STATE
      INTEGER n
#endif
CEOP

C     Loop over tracers
      DO iTracer = 1, PTRACERS_num

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C added by Yilang Xu on April 8, 2020
C to link surface source of ptracers to sea ice area
C Reset the surface tracer at every time step
C For passive tracer #1 only
#ifdef ALLOW_SEAICE
        DO bj = myByLo(myThid), myByHi(myThid)
          DO bi = myBxLo(myThid), myBxHi(myThid)
            DO j=1-OLy,sNy+OLy
              DO i=1-OLx,sNx+OLx
                ig = myXGlobalLo-1+(bi-1)*sNx+i
                jg = myYGlobalLo-1+(bj-1)*sNy+j
                IF ((ig .GE. 121) .AND. (ig .LE. 226)
     &           .AND. (jg .GE. 101) .AND. (jg .LE. 150)
     &           .AND. (AREA(i,j,bi,bj) .LE. 0.70 _d 0)) THEN
C                 pTracer(i,j,1,bi,bj,iTracer)=1. _d 0
                 pTracer(i,j,1,bi,bj,1)=1. _d 0
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
#endif /* ALLOW_SEAICE */
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C     Check if it is time to reset this tracer
       IF ( PTRACERS_resetFreq(iTracer).GT.0. .AND. myIter.GT.0 .AND.
     &  DIFFERENT_MULTIPLE( PTRACERS_resetFreq(iTracer),
     &  myTime + PTRACERS_resetPhase(iTracer), deltaTClock ) ) THEN

C     message
        _BEGIN_MASTER( myThid )
        WRITE(msgBuf,'(A,I2,I10)')
     &           '// PTRACER Resetting, (iTracer,t-step) = ',
     &            iTracer, myIter
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
        _END_MASTER( myThid )

C     Initialise again this tracer arrays
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
          DO k=1,Nr
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
             pTracer(i,j,k,bi,bj,iTracer) = PTRACERS_ref(k,iTracer)
             gpTrNm1(i,j,k,bi,bj,iTracer) = 0. _d 0
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO

C     Read initial conditions and exchange
        tmpInitialFile = ' '
        tmpInitialFile = PTRACERS_initialFile(iTracer)
        IF ( tmpInitialFile .NE. ' ' ) THEN
         CALL READ_FLD_XYZ_RL(tmpInitialFile,' ',
     &        pTracer(1-OLx,1-OLy,1,1,1,iTracer),0,myThid)
         _EXCH_XYZ_RL(pTracer(1-OLx,1-OLy,1,1,1,iTracer),myThid)
        ENDIF

C     Apply mask and reset tendencies
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
          DO k=1,Nr
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
             IF (maskC(i,j,k,bi,bj).EQ.0.)
     &           pTracer(i,j,k,bi,bj,iTracer)=0. _d 0
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO

#ifdef PTRACERS_ALLOW_DYN_STATE
C     Initialize SOM array :
        IF ( PTRACERS_SOM_Advection(iTracer) ) THEN
        _BEGIN_MASTER( myThid )
         WRITE(msgBuf,'(A,I3,A)')'PTRACERS_RESET: iTracer = ',
     &               iTracer, ' : resetting 2nd-order moments '
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT, myThid )
        _END_MASTER( myThid )
         DO bj = myByLo(myThid), myByHi(myThid)
          DO bi = myBxLo(myThid), myBxHi(myThid)
           DO n = 1,nSOM
            DO k=1,Nr
             DO j=1-OLy,sNy+OLy
              DO i=1-OLx,sNx+OLx
               _Ptracers_som(i,j,k,bi,bj,n,iTracer) = 0. _d 0
              ENDDO
             ENDDO
            ENDDO
           ENDDO
          ENDDO
         ENDDO
c        CALL GAD_EXCH_SOM( _Ptracers_som(:,:,:,:,:,:,iTracer),
c    &                       Nr, myThid )
        ENDIF
#endif /* PTRACERS_ALLOW_DYN_STATE */

        _BEGIN_MASTER( myThid )
        WRITE(msgBuf,'(A)')  ' '
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
        _END_MASTER( myThid )

C     end of reset if block
       ENDIF
C     end of Tracer loop
      ENDDO

#endif /* ALLOW_PTRACERS */

      RETURN
      END

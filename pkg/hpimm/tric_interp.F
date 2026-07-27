C POLY_INTERP: subroutines for polynomial interpolation and evaluation

C=====================================================================
C  VAND4SOLVE  --  Solve the 4x4 Vandermonde system  V * C = B
C
C  V(m,n) = T(m)^(n-1),  m,n = 1,2,3,4.
C  Method: Gaussian elimination with partial pivoting.
C  B is not modified on return (a local copy R is used).
C
C  IN  : T(4) -- nodes
C        B(4) -- right-hand side
C  OUT : C(4) -- solution vector
C=====================================================================
      SUBROUTINE VAND4SOLVE(T, B, C)
      IMPLICIT NONE
      DOUBLE PRECISION T(4), B(4), C(4)
      DOUBLE PRECISION V(4,4), R(4), FAC, TMP, PMAX
      INTEGER          I, J, K, IP

C --- Build the Vandermonde matrix; copy RHS to local array R.
      DO 20 I = 1, 4
        V(I,1) = 1.0D0
        DO 10 J = 2, 4
          V(I,J) = V(I,J-1) * T(I)
   10   CONTINUE
        R(I) = B(I)
   20 CONTINUE

C --- Forward elimination with partial pivoting.
      DO 70 K = 1, 3

C       Find the pivot row in column K.
        IP   = K
        PMAX = ABS(V(K,K))
        DO 30 I = K+1, 4
          IF (ABS(V(I,K)) .GT. PMAX) THEN
            PMAX = ABS(V(I,K))
            IP   = I
          END IF
   30   CONTINUE

C       Swap rows K and IP if required.
        IF (IP .NE. K) THEN
          DO 40 J = K, 4
            TMP      = V(K,J)
            V(K,J)   = V(IP,J)
            V(IP,J)  = TMP
   40     CONTINUE
          TMP    = R(K)
          R(K)   = R(IP)
          R(IP)  = TMP
        END IF

C       Eliminate entries below the pivot.
        DO 60 I = K+1, 4
          FAC    = V(I,K) / V(K,K)
          DO 50 J = K+1, 4
            V(I,J) = V(I,J) - FAC * V(K,J)
   50     CONTINUE
          R(I)   = R(I) - FAC * R(K)
          V(I,K) = 0.0D0
   60   CONTINUE

   70 CONTINUE

C --- Back substitution.
      C(4) = R(4) / V(4,4)
      DO 90 I = 3, 1, -1
        TMP = R(I)
        DO 80 J = I+1, 4
          TMP = TMP - V(I,J) * C(J)
   80   CONTINUE
        C(I) = TMP / V(I,I)
   90 CONTINUE

      RETURN
      END

C=====================================================================
C  TRIC_SETUP  --  Precompute monomial coefficients of the tricubic
C                  interpolant on a 4x4x4 rectangular grid.
C
C  The polynomial is stored as:
C    p(x,y,z) = SUM_{i,j,k=1}^{4} A(i,j,k)*x^(i-1)*y^(j-1)*z^(k-1)
C
C  Method: three successive passes of 16 Vandermonde solves (one per
C  axis), exploiting the Kronecker product structure of the grid.
C  Total cost: 48 tiny (4x4) systems.
C
C  IN  : XG(4), YG(4), ZG(4)  -- axis node coordinates
C        F(4,4,4)              -- F(i,j,k) = f(XG(i), YG(j), ZG(k))
C  OUT : A(4,4,4)              -- monomial coefficient tensor
C=====================================================================
      SUBROUTINE TRIC_SETUP(XG, YG, ZG, F, A)
      IMPLICIT NONE
      DOUBLE PRECISION XG(4), YG(4), ZG(4)
      DOUBLE PRECISION F(4,4,4), A(4,4,4)
      DOUBLE PRECISION W(4,4,4), COL(4), COEF(4)
      INTEGER          I, J, K

C --- Pass 1 (x-direction): 16 Vandermonde solves.
C     W(i,j,k) = coefficient of x^(i-1) in the 1-D cubic
C     interpolant through {(XG(m), F(m,j,k)), m=1..4}.
      DO 130 K = 1, 4
        DO 120 J = 1, 4
          DO 100 I = 1, 4
            COL(I) = F(I,J,K)
  100     CONTINUE
          CALL VAND4SOLVE(XG, COL, COEF)
          DO 110 I = 1, 4
            W(I,J,K) = COEF(I)
  110     CONTINUE
  120   CONTINUE
  130 CONTINUE

C --- Pass 2 (y-direction): 16 Vandermonde solves.
C     A(i,j,k) [temp] = coefficient of x^(i-1)*y^(j-1).
      DO 230 K = 1, 4
        DO 220 I = 1, 4
          DO 200 J = 1, 4
            COL(J) = W(I,J,K)
  200     CONTINUE
          CALL VAND4SOLVE(YG, COL, COEF)
          DO 210 J = 1, 4
            A(I,J,K) = COEF(J)
  210     CONTINUE
  220   CONTINUE
  230 CONTINUE

C --- Pass 3 (z-direction): 16 Vandermonde solves.
C     A(i,j,k) [final] = coefficient of x^(i-1)*y^(j-1)*z^(k-1).
C     In-place overwrite of A is safe: the entire z-column is copied
C     into COL before any element of A is modified.
      DO 330 I = 1, 4
        DO 320 J = 1, 4
          DO 300 K = 1, 4
            COL(K) = A(I,J,K)
  300     CONTINUE
          CALL VAND4SOLVE(ZG, COL, COEF)
          DO 310 K = 1, 4
            A(I,J,K) = COEF(K)
  310     CONTINUE
  320   CONTINUE
  330 CONTINUE

      RETURN
      END

C=====================================================================
C  TRIC_EVAL  --  Evaluate the tricubic polynomial or a mixed partial
C                 derivative at a single point (X, Y, Z).
C
C  Computes  VAL = d^(NX+NY+NZ) p / (dx^NX dy^NY dz^NZ) at (X,Y,Z)
C  where  p(x,y,z) = SUM_{i,j,k=1}^{4} A(i,j,k)*x^(i-1)*y^(j-1)*z^(k-1).
C
C  Algorithm: three nested Horner passes (x innermost, z outermost).
C  Falling-factorial prefactors [arising from differentiation] are
C  precomputed once per call before the Horner sweeps.
C
C     FFD(I) = (I-1)*(I-2)*...*(I-ND)   [ND factors; = 1 when ND=0]
C
C  Terms with power < ND vanish upon differentiation; the Horner
C  loops start at index ND+1, so they are never visited.
C
C  Cost: O(1) falling-factorial precomputation + O(64) multiply-adds.
C  VAL = 0 automatically if any derivative order exceeds 3.
C
C  IN  : A(4,4,4)   -- monomial coefficients (output of TRIC_SETUP)
C        X, Y, Z    -- evaluation point
C        NX, NY, NZ -- derivative orders (integers >= 0)
C  OUT : VAL        -- value of the requested (mixed) derivative
C=====================================================================
      SUBROUTINE TRIC_EVAL(A, X, Y, Z, NX, NY, NZ, VAL)
      IMPLICIT NONE
      DOUBLE PRECISION A(4,4,4), X, Y, Z, VAL
      INTEGER          NX, NY, NZ
      DOUBLE PRECISION CX(4,4), CY(4)
      DOUBLE PRECISION FFX(4), FFY(4), FFZ(4)
      DOUBLE PRECISION ACCX, ACCY, ACCZ
      INTEGER          I, J, K, M

C --- Precompute falling-factorial prefactors for each axis.
C     Only entries that will actually be accessed in the Horner
C     sweeps (index >= derivative order + 1) are initialised.
      DO 15 I = NX+1, 4
        FFX(I) = 1.0D0
        DO 10 M = 0, NX-1
          FFX(I) = FFX(I) * DBLE(I - 1 - M)
   10   CONTINUE
   15 CONTINUE

      DO 25 J = NY+1, 4
        FFY(J) = 1.0D0
        DO 20 M = 0, NY-1
          FFY(J) = FFY(J) * DBLE(J - 1 - M)
   20   CONTINUE
   25 CONTINUE

      DO 35 K = NZ+1, 4
        FFZ(K) = 1.0D0
        DO 30 M = 0, NZ-1
          FFZ(K) = FFZ(K) * DBLE(K - 1 - M)
   30   CONTINUE
   35 CONTINUE

C --- Pass 1 (x-direction): for each (J,K), sweep I from 4 down to
C     NX+1, evaluating the NX-th x-derivative of the row polynomial
C     at X by Horner's rule.  Result in CX(J,K).
      DO 60 K = 1, 4
        DO 55 J = 1, 4
          ACCX = 0.0D0
          DO 50 I = 4, NX+1, -1
            ACCX = ACCX * X + A(I,J,K) * FFX(I)
   50     CONTINUE
          CX(J,K) = ACCX
   55   CONTINUE
   60 CONTINUE

C --- Pass 2 (y-direction): for each K, sweep J from 4 down to NY+1,
C     evaluating the NY-th y-derivative at Y using CX(J,K) as the
C     j-indexed coefficient sequence.  Result in CY(K).
      DO 90 K = 1, 4
        ACCY = 0.0D0
        DO 80 J = 4, NY+1, -1
          ACCY = ACCY * Y + CX(J,K) * FFY(J)
   80   CONTINUE
        CY(K) = ACCY
   90 CONTINUE

C --- Pass 3 (z-direction): sweep K from 4 down to NZ+1, evaluating
C     the NZ-th z-derivative at Z using CY(K) as coefficients.
      ACCZ = 0.0D0
      DO 110 K = 4, NZ+1, -1
        ACCZ = ACCZ * Z + CY(K) * FFZ(K)
  110 CONTINUE
      VAL = ACCZ

      RETURN
      END

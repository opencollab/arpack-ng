c     Initialisation of the stat common block.
      block data staini
      common /timing/
     &     nopx, nbx, nrorth, nitref, nrstrt,
     &     tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     &     tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     &     tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     &     tmvopx, tmvbx, tgetv0, titref, trvec
      data nopx, nbx, nrorth, nitref, nrstrt,
     &     tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     &     tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     &     tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     &     tmvopx, tmvbx, tgetv0, titref, trvec
     &     / 0, 0, 0, 0, 0,
     &     0., 0., 0., 0., 0., 0., 0.,
     &     0., 0., 0., 0., 0., 0., 0.,
     &     0., 0., 0., 0., 0., 0., 0.,
     &     0., 0., 0., 0., 0. /
      end

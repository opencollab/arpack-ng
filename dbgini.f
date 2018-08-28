c Initialisation of the debug common block to "no debug".
      block data dbgini
      common /debug/ logfil, ndigit, mgetv0,
     &     msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     &     mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     &     mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
      data logfil, ndigit, mgetv0,
     &     msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     &     mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     &     mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
     &     / 6, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     &     0, 0, 0, 0, 0, 0 /
      end

#!/bin/bash -eu

for eigPb in "--A As.mtx" "--nonSymPb --A An.mtx" "--nonSymPb --cpxPb --A Az.mtx --B Bz.mtx"
do
  for genPb in "" "--genPb"
  do
    for smallMag in "" "--mag SM --noCheck" # SM is known to be difficult to converge.
    do
      export shiftOpt=""
      if [[ "$eigPb" == *nonSymPb* ]]; then
        if [[ "$genPb" == *genPb* ]]; then
          export shiftOpt="--shiftReal 2.5 --shiftImag 2.5 --tol 0.5" # Relax tolerance, tricky to converge.
        else
          export shiftOpt="--shiftReal 100.0 --shiftImag 100.0"
        fi
      else
        if [[ "$genPb" == *genPb* ]]; then
          export shiftOpt="--shiftReal 50.0"
        else
          export shiftOpt="--shiftReal 100.0"
        fi
      fi

      for shiftRI in "" "$shiftOpt"
      do
        for invert in "" "--invert"
        do
          for tol in "" "--tol 1.e-5"
          do
            for slv in "" "--slv CG --slvItrTol 1.e-06 --slvItrMaxIt 100" "--slv LU --slvDrtPvtThd 1.e-06" "--slv QR"
            do
              for rs in "" "--schur"
              do
                export extraGenPb=""
                if [[ "$genPb" == *genPb* ]]; then
                  export extraGenPb="$shiftOpt" # Force shift if genPb.
                fi

                if [[ "$slv" == *CG* ]]; then
                  if [[ "$eigPb" == *nonSymPb* ]]; then
                    continue # Skip CG that could fail (CG is meant to deal with sym matrices).
                  fi
                fi

                # Run arpackmm: use --nbCV 6 to ease convergence, and, --verbose 3 for debug.
                export CMD="./arpackmm $eigPb $genPb $smallMag $shiftRI $invert $tol $slv $extraGenPb $rs --nbCV 6 --verbose 3"
                echo "$CMD"
                eval "$CMD"
                echo ""
                echo "========================================================================================"
                echo ""

                # Run arpackmm: re-run with restart.
                export CMD="$CMD --restart"
                echo "$CMD"
                eval "$CMD"
                echo ""
                echo "========================================================================================"
                echo ""
              done
            done
          done
        done
      done
    done
  done
done

echo "OK"

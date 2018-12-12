#!/bin/bash -eu

for symPb in "--A As.mtx" "--nonSymPb --A An.mtx"
do
  for genPb in "" "--genPb"
  do
    for smallMag in "" "--mag SM --noCheck" # SM is known to be difficult to converge.
    do
      export shiftOpt=""
      if [[ "$symPb" == *nonSymPb* ]]; then
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
            for slv in "" "--slv CG" "--slv LU" "--slv QR"
            do
              export extraGenPb=""
              if [[ "$genPb" == *genPb* ]]; then
                export extraGenPb="$shiftOpt" # Force shift if genPb.
              fi

              # Run arpackmm: use --nbCV 6 to ease convergence, and, --verbose 3 for debug.
              export CMD="./arpackmm $symPb $genPb $smallMag $shiftRI $invert $tol $slv $extraGenPb --nbCV 6 --verbose 3"
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

echo "OK"

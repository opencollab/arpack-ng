#!/bin/bash -eu

trap 'catch' ERR
catch() {
  grep OPT arpackmm.run.log
  grep Error arpackmm.run.log
  exit 1
}

# For all these eigen problems, the first eigen value is about 382 or (382, 0.).

for stdPb in "--A As.mtx" "--nonSymPb --A An.mtx" "--nonSymPb --cpxPb --A Az.mtx"
do
  # Choose B matrix according to A.
  export fileB=""
  if [[ "$stdPb" == *cpxPb* ]]; then
    export fileB="--B Bz.mtx"
  else
    export fileB="--B B.mtx"
  fi

  for genPb in "" "$fileB"
  do
    # SM may not converge so we can not apply checks.
    # LM + invert is equivalent to SM and does converge.
    # Note: this is expected as "power-like" methods are designed to find largest eigen values.
    for magOpt in "--mag LM" "--mag SM --noCheck"
    do
      export mrn="--maxResNorm 1.e-1" # Relax residual check to get stable tests.

      # Shift slightly to avoid the zero-vector starting problem.
      export shiftZV=""
      if [[ "$stdPb" == *cpxPb* ]]; then
        export shiftZV="--shiftReal 1.0 --shiftImag 1.0"
      else
        export shiftZV="--shiftReal 1.0"
      fi

      # Shift according to the estimation of the eigen value we may have.
      export shiftEV=""
      if [[ "$stdPb" == *cpxPb* ]]; then
        export shiftEV="--shiftReal 380.0 --shiftImag 1.0"
      else
        export shiftEV="--shiftReal 380.0"
      fi

      for shiftRI in "$shiftZV" "$shiftEV"
      do
        for invert in "" "--invert"
        do
          # Choose solver according to the sym/non-sym type of the problem.
          export cgSlv=""
          if [[ "$stdPb" == *nonSymPb* ]]; then
            export cgSlv="BiCG"
          else
            export cgSlv="CG"
          fi

          for slv in "--slv $cgSlv" "--slv LU"
          do
            for rs in "" "--schur"
            do
              for dsPrec in "" "--simplePrec"
              do
                for dsMat in "" "--dense true"
                do
                  # Skip not supported cases.
                  if [[ "$dsMat" == *dense* ]]; then
                    if [[ "$slv" == *CG* ]]; then
                      continue # Iterative solvers are not allowed when using dense matrices.
                    fi
                  fi

                  export easeCV="--nbCV 6 --maxIt 200" # Use --nbCV 6 and --maxIt 200 to ease convergence.
                  echo "CLI: ./arpackmm $stdPb $genPb $magOpt $mrn $shiftRI $invert $slv $rs $dsPrec $dsMat $easeCV"
                  echo "----------------------------------------------------------------------------------------"

                  # Run arpackmm.
                  ./arpackmm "$stdPb" "$genPb" "$magOpt" "$mrn" "$shiftRI" "$invert" "$slv" "$rs" "$dsPrec" "$dsMat" \
                             "$easeCV" --verbose 3 &> arpackmm.run.log
                  grep OPT arpackmm.run.log
                  grep OUT arpackmm.run.log
                  echo "----------------------------------------------------------------------------------------"

                  # Run arpackmm: re-run with restart always with small shift to avoid the zero-starting vector problem.
                  ./arpackmm "$stdPb" "$genPb" "$magOpt" "$mrn" "$shiftZV" "$invert" "$slv" "$rs" "$dsPrec" "$dsMat" \
                             --restart                                                                               \
                             "$easeCV" --verbose 3 &> arpackmm.run.log
                  grep OPT arpackmm.run.log
                  grep OUT arpackmm.run.log
                  echo "========================================================================================"
                done
              done
            done
          done
        done
      done
    done
  done
done

echo "arpackmm: OK"

(* ================================================================
   vmmc_2d_field.wl
   2D VMMC with user-defined field and coupling functions on a
   fully-periodic (torus) square lattice.
   ================================================================

   The field f(x,y) and coupling J(a,b,d²) are defined as ordinary
   Mathematica functions in Section 0 below.  The checker proves
   detailed balance symbolically for ALL values of any free parameters
   in these functions simultaneously.

   This example defines:
     fieldF     – sinusoidal field in the x-direction
     couplingJ  – exponential decay with a per-type-pair amplitude

   To use a different model, edit Section 0 only.

   STATE / ENCODING
      Identical to vmmc_2d.wl: flat L²-array, state[[s]]=0 (empty) or
      k∈{1,...,N} (labeled particle).  BitsToState filters to perfect-
      square lengths.

   REFERENCES
      S. Whitelam & P. L. Geissler, J. Chem. Phys. 127, 154101 (2007)
   ================================================================ *)


(* ================================================================
   SECTION 0 — USER INTERFACE: field and coupling functions
   ================================================================

   Edit this section to define your physical model.

   fieldF[x, y, L]
     Field energy felt by any particle at lattice column x, row y
     (both 0-indexed integers, x,y ∈ {0,...,L-1}) on an L×L lattice.
     Evaluated with integer arguments; use _Integer patterns.

   couplingJ[a, b, d2]
     Pair interaction between particle types a and b at squared
     distance d2 (always a positive integer on an integer lattice).
     MUST satisfy couplingJ[a,b,d2] = couplingJ[b,a,d2].
     The checker verifies this before running; an asymmetric coupling
     immediately implies a detailed-balance violation.
     Use $jPairSym[a,b] for a per-type-pair coupling symbol
     (it is symmetric by construction: $jPairSym[a,b] = $jPairSym[b,a]).

   $fieldFreeParams
     List every free symbolic parameter that appears in fieldF
     (do not include x, y, L — those are always concrete integers).

   $couplingFreeParams
     List every free symbolic parameter that appears in couplingJ
     beyond the per-pair $jPairSym symbols (which are auto-discovered).

   $maxD2
     Maximum squared distance included in bonds and neighbour shells.
     1 = nearest-neighbour only (fastest checker, d = 1).
     5 = includes d = 1, √2, 2, √5 (more realistic, slower).
     The coupling function is evaluated at all active bond distances,
     so setting $maxD2 = 1 while couplingJ is nonzero at larger d2 is
     valid — it simply restricts the model to nearest-neighbour bonds.
   ================================================================ *)

(* ---- EXAMPLE FIELD: sinusoidal in the x (column) direction ---- *)
(* fieldAmp is the field amplitude — a single free symbolic parameter.
   For an L×L lattice with x ∈ {0,...,L-1}:
     x=0 → fieldAmp·Sin[π/L]
     x=1 → fieldAmp·Sin[2π/L]
     ...
   Mathematica evaluates Sin at rational multiples of π exactly. *)
fieldF[x_Integer, y_Integer, L_Integer] :=
  fieldAmp * Sin[\[Pi] * (x + 1) / L]

(* ---- EXAMPLE COUPLING: exponential decay with per-pair amplitude ---- *)
(* lambdaJ is the decay rate (one shared free parameter).
   $jPairSym[a,b] is an independent amplitude symbol per type pair.
   couplingJ[a,b,d2] = Jpair<lo><hi> · exp(−λ d²)
   This is symmetric: couplingJ[a,b,d2] = couplingJ[b,a,d2]. ✓ *)
couplingJ[a_Integer, b_Integer, d2_Integer] :=
  $jPairSym[a, b] * Exp[-lambdaJ * d2]

(* Free parameters for the checker
   (per-pair $jPairSym symbols are auto-added by DynamicSymParams) *)
$fieldFreeParams    = {fieldAmp};
$couplingFreeParams = {lambdaJ};

(* Nearest-neighbour only by default for checker speed *)
$maxD2 = 1;


(* ================================================================
   SECTION 1 — Bijective integer encoding   (identical to vmmc_2d.wl)
   ================================================================ *)

$cL[L_]        := $cL[L]       = Sum[Binomial[L, k] * k!, {k, 0, L}]
$cLPre[L_]     := $cLPre[L]    = Sum[$cL[l], {l, 0, L - 1}]
$cLNPre[L_,N_] := $cLNPre[L,N] = Sum[Binomial[L, k] * k!, {k, 0, N - 1}]

$rankCombo[pos_List] := Sum[Binomial[pos[[i]], i], {i, Length[pos]}]

$unrankCombo[rank_, L_, N_] :=
  Module[{pos = ConstantArray[0, N], x = L - 1, r = rank},
    Do[While[Binomial[x, i] > r, x--]; pos[[i]] = x; r -= Binomial[x, i]; x--,
       {i, N, 1, -1}]; pos]

$rankPerm[perm_List] :=
  Module[{n = Length[perm], elems = Range[Length[perm]], rank = 0, idx},
    Do[idx = FirstPosition[elems, perm[[i]]][[1]] - 1;
       rank += idx * Factorial[n - i]; elems = Delete[elems, idx + 1],
       {i, n}]; rank]

$unrankPerm[k_, n_] :=
  Module[{elems = Range[n], perm = {}, r = k, idx},
    Do[idx = Quotient[r, Factorial[i - 1]]; r = Mod[r, Factorial[i - 1]];
       AppendTo[perm, elems[[idx + 1]]]; elems = Delete[elems, idx + 1],
       {i, n, 1, -1}]; perm]

$decode[id_Integer] :=
  Module[{L = 0, N = 0, r, rpos, rperm, pos, perm, arr},
    While[$cLPre[L + 1] <= id, L++];
    r = id - $cLPre[L];
    While[$cLNPre[L, N + 1] <= r, N++];
    r -= $cLNPre[L, N];
    rpos = Quotient[r, Factorial[N]]; rperm = Mod[r, Factorial[N]];
    pos  = $unrankCombo[rpos, L, N]; perm  = $unrankPerm[rperm, N];
    arr  = ConstantArray[0, L];
    Do[arr[[pos[[i]] + 1]] = perm[[i]], {i, N}]; arr]


(* ================================================================
   SECTION 2 — Periodic lattice helpers   (identical to vmmc_2d.wl)
   ================================================================ *)

(* Translate site s by direction vector {dRow, dCol} on the L×L torus *)
$applyDir[s_, {dr_, dc_}, L_] :=
  Mod[Ceiling[s/L] - 1 + dr, L]*L + Mod[Mod[s-1, L] + dc, L] + 1

(* Row and column of site s (1-indexed) *)
$row[s_, L_] := Ceiling[s/L]
$col[s_, L_] := Mod[s-1, L] + 1

(* 0-indexed coordinates for fieldF (x = column index, y = row index) *)
$siteX[s_, L_] := $col[s, L] - 1   (* x ∈ {0,...,L-1} *)
$siteY[s_, L_] := $row[s, L] - 1   (* y ∈ {0,...,L-1} *)


(* ================================================================
   SECTION 3 — Minimum-image squared distance on the L×L torus
   ================================================================ *)

$torusD2[s1_, s2_, L_] :=
  With[{dr0 = Abs[$row[s1,L] - $row[s2,L]],
         dc0 = Abs[$col[s1,L] - $col[s2,L]]},
    With[{dr = Min[dr0, L - dr0],
           dc = Min[dc0, L - dc0]},
      dr^2 + dc^2]]


(* ================================================================
   SECTION 4 — Neighbour shells (memoised)
   ================================================================

   Includes all sites within squared distance $maxD2.
   Restricting this shell to active bond distances is critical for
   checker speed: a site with zero coupling would still trigger a
   RandomReal[] draw in the cluster builder (consuming a BFS bit)
   even though the link weight is always 0. *)

$neighborsD2[s_, L_] := $neighborsD2[s, L] =
  Select[Range[L^2],
    Function[q, Module[{d2 = $torusD2[s, q, L]}, d2 > 0 && d2 <= $maxD2]]]


(* ================================================================
   SECTION 5 — Unique undirected bonds (memoised)
   ================================================================

   Each bond {s1, s2, d²} with s1 < s2 at squared distance ≤ $maxD2.
   Used by energy[] to sum all pair contributions. *)

$uniqueBondsExt[L_] := $uniqueBondsExt[L] =
  Flatten[
    Table[
      With[{d2 = $torusD2[s1, s2, L]},
        If[d2 > 0 && d2 <= $maxD2, {{s1, s2, d2}}, {}]],
      {s1, L^2}, {s2, s1+1, L^2}],
    2]


(* ================================================================
   SECTION 6 — Energy function
   ================================================================

   Total energy = pair energy + field energy.
   Pair energy: sum over all undirected bonds within $maxD2, using
     couplingJ[type_a, type_b, d²].
   Field energy: each occupied site s contributes fieldF[x, y, L]
     where x = $siteX[s,L] and y = $siteY[s,L] are 0-indexed integers. *)

energy[state_List] :=
  With[{L = Round[Sqrt[Length[state]]]},
    (* Pair energy: sum over bonds, skip holes *)
    Total[Map[Function[bond,
      With[{t1 = state[[bond[[1]]]], t2 = state[[bond[[2]]]], d2 = bond[[3]]},
        If[t1 != 0 && t2 != 0, couplingJ[t1, t2, d2], 0]]],
      $uniqueBondsExt[L]]] +
    (* Field energy: each occupied site contributes fieldF *)
    Total @ Table[
      If[state[[s]] != 0,
         fieldF[$siteX[s, L], $siteY[s, L], L],
         0],
      {s, L^2}]]


(* ================================================================
   SECTION 7 — Virtual pair energy for VMMC link weights
   ================================================================

   Pair energy between a particle of type typeI at virtual site vI
   and a particle of type typeJ at site qSite.
   Uses couplingJ at the torus distance; Infinity at hard-core overlap. *)

$virtualPairEnergy[typeI_, typeJ_, vI_, qSite_, L_] :=
  Which[
    vI === qSite, Infinity,
    True,
      With[{d2 = $torusD2[vI, qSite, L]},
        If[d2 > 0 && d2 <= $maxD2, couplingJ[typeI, typeJ, d2], 0]]]


(* ================================================================
   SECTION 8 — VMMC cluster builder
   ================================================================

   Whitelam-Geissler link-weight logic.  Identical structure to
   vmmc_2d.wl; uses $neighborsD2 (respects $maxD2) and
   $virtualPairEnergy (calls couplingJ). *)

$vmmcBuildCluster[state_, L_, seed_, dir_] :=
  Module[{
    cluster   = {seed},
    inCluster = <|seed -> True|>,
    queue     = {seed},
    frustrated = False,
    p, pType, pPost, pRev, nbrs, q, qType,
    eInit, eFwd, eRev, wFwd, wRev, r1, r2
  },

    While[queue =!= {} && !frustrated,
      p     = First[queue]; queue = Rest[queue];
      pType = state[[p]];
      pPost = $applyDir[p,  dir, L];
      pRev  = $applyDir[p, {-dir[[1]], -dir[[2]]}, L];

      (* Union of neighbour shells of p, pPost, and pRev *)
      nbrs = DeleteDuplicates @ Join[
               $neighborsD2[p,     L],
               $neighborsD2[pPost, L],
               $neighborsD2[pRev,  L]];

      Do[
        q = nbrs[[k]];
        If[state[[q]] =!= 0 && !KeyExistsQ[inCluster, q],
          qType = state[[q]];

          eInit = $virtualPairEnergy[pType, qType, p,     q, L];
          eFwd  = $virtualPairEnergy[pType, qType, pPost, q, L];
          eRev  = $virtualPairEnergy[pType, qType, pRev,  q, L];

          wFwd = Piecewise[{
              {1,                                eFwd === Infinity},
              {1 - Exp[\[Beta] (eInit - eFwd)],  eInit < eFwd}},
            0];
          wRev = Piecewise[{
              {1,                                eRev === Infinity},
              {1 - Exp[\[Beta] (eInit - eRev)],  eInit < eRev}},
            0];

          r1 = RandomReal[];
          If[r1 <= wFwd,
            r2 = RandomReal[];
            If[r2 > Piecewise[{
                  {1,
                      eFwd === Infinity && eRev === Infinity},
                  {1 - Exp[\[Beta] (eInit - eRev)],
                      eFwd === Infinity && eInit < eRev},
                  {0,
                      eFwd === Infinity},
                  {1,
                      eRev === Infinity && eInit < eFwd},
                  {Min[(1 - Exp[\[Beta] (eInit - eRev)]) /
                        (1 - Exp[\[Beta] (eInit - eFwd)]), 1],
                   eInit < eFwd && eInit < eRev},
                  {0,
                      eInit < eFwd}},
                0],
              frustrated = True; Break[],
              AppendTo[cluster, q];
              inCluster[q] = True;
              AppendTo[queue, q]
            ]
          ]
        ],
        {k, Length[nbrs]}
      ]
    ];

    If[frustrated, None, cluster]
  ]


(* ================================================================
   SECTION 9 — Algorithm
   ================================================================ *)

Algorithm[state_List] :=
  Module[{L, occupied, seed, dirs, dir, cluster, newState, dEField, dE},
    L = Round[Sqrt[Length[state]]];

    occupied = Flatten[Position[state, _?(# > 0 &)]];
    If[occupied === {}, Return[state]];
    seed = RandomChoice[occupied];

    dirs = {{0,1},{0,-1},{1,0},{-1,0}};
    dir  = RandomChoice[dirs];

    cluster = $vmmcBuildCluster[state, L, seed, dir];
    If[cluster === None, Return[state]];

    newState = state;
    Do[newState[[cluster[[i]]]] = 0, {i, Length[cluster]}];
    Do[With[{dest = $applyDir[cluster[[i]], dir, L]},
         If[newState[[dest]] =!= 0, Return[state, Module]];
         newState[[dest]] = state[[cluster[[i]]]]
       ], {i, Length[cluster]}];

    (* Post-cluster Metropolis for field energy change.
       ΔE_field = Σ_{p in cluster} [fieldF(dest(p)) − fieldF(p)]
       fieldF and $siteX/$siteY always produce integer arguments. *)
    dEField = Total @ Table[
      fieldF[$siteX[$applyDir[cluster[[i]], dir, L], L],
             $siteY[$applyDir[cluster[[i]], dir, L], L], L] -
      fieldF[$siteX[cluster[[i]], L],
             $siteY[cluster[[i]], L], L],
      {i, Length[cluster]}];
    If[RandomReal[] >= MetropolisProb[dEField], Return[state]];

    dE = energy[newState] - energy[state];
    MetropolisProb[dE];

    newState
  ]


(* ================================================================
   SECTION 10 — Checker interface   (identical to vmmc_2d.wl)
   ================================================================ *)

BitsToState[bits_List] :=
  Module[{id = FromDigits[bits, 2], state, sqrtM},
    If[id == 0, Return[None]];
    state = $decode[id];
    sqrtM = Sqrt[Length[state]];
    If[!IntegerQ[sqrtM], Return[None]];
    state]

DisplayState[state_List] :=
  With[{L = Round[Sqrt[Length[state]]]},
    StringJoin @ Riffle[
      Table["{" <> StringRiffle[ToString /@ state[[(r-1)*L+1 ;; r*L]], ","] <> "}",
            {r, 1, L}],
      "|"]]

ValidStateIDs[maxId_Integer] :=
  Module[{L = 1, ids = {}},
    While[$cLPre[L^2] <= maxId,
      ids = Join[ids, Range[$cLPre[L^2], Min[$cLPre[L^2+1]-1, maxId]]];
      L++];
    ids]

(* DynamicSymParams: returns the free symbolic parameters for this
   component.  The per-pair $jPairSym symbols are auto-generated from
   the particle types present; the user's $fieldFreeParams and
   $couplingFreeParams supply any additional free parameters. *)
DynamicSymParams[states_List] :=
  Module[{types = Sort[DeleteCases[Union @@ states, 0]], jPairSyms},
    jPairSyms = Flatten @ Table[
      If[a < b, {$jPairSym[a, b]}, Nothing],
      {a, types}, {b, types}];
    <|"couplings" -> Join[$fieldFreeParams, $couplingFreeParams, jPairSyms]|>]

numBeta = 1
